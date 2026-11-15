mod db;
mod vcf;

use axum::{
    extract::{Query, State},
    http::{header, StatusCode},
    response::{Html, IntoResponse, Response},
    routing::{get, post},
    Json, Router,
};
use rusqlite::params;
use serde::Deserialize;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::Instant;

struct AppState {
    conn: Mutex<rusqlite::Connection>,
    reference: vcf::Reference,
}

const INDEX_HTML: &str = include_str!("../static/index.html");

#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |k: &str, d: &str| -> String {
        args.windows(2)
            .find(|w| w[0] == k)
            .map(|w| w[1].clone())
            .unwrap_or_else(|| d.to_string())
    };
    let db_path = arg("--db", "variant.db");
    let ref_path = arg("--ref", "ref.fa");
    let genes_path = arg("--genes", "");
    let port = arg("--port", "8000");

    let reference = vcf::Reference::from_fasta(&ref_path).expect("load reference failed");
    let conn = rusqlite::Connection::open(&db_path).expect("open db failed");
    db::init_schema(&conn).expect("init schema failed");
    if !genes_path.is_empty() {
        let n = db::load_genes(&conn, &genes_path).expect("load genes failed");
        println!("loaded {n} genes");
    }
    let state = Arc::new(AppState {
        conn: Mutex::new(conn),
        reference,
    });

    let app = Router::new()
        .route("/", get(|| async { Html(INDEX_HTML) }))
        .route("/api/stats", get(api_stats))
        .route("/api/samples", get(api_samples))
        .route("/api/duplicates", get(api_duplicates))
        .route("/api/batches", get(api_batches))
        .route("/api/import", post(api_import))
        .route("/api/variants", get(api_variants))
        .route("/api/frequencies", get(api_frequencies))
        .route("/api/qc", get(api_qc))
        .route("/api/export.vcf", get(api_export_vcf))
        .route("/api/bench", get(api_bench))
        .with_state(state);

    let listener = tokio::net::TcpListener::bind(format!("0.0.0.0:{port}"))
        .await
        .unwrap();
    println!("variantdb listening on http://0.0.0.0:{port}");
    axum::serve(listener, app).await.unwrap();
}

// ---------- 通用过滤条件 ----------

#[derive(Deserialize, Default)]
struct Filter {
    chrom: Option<String>,
    start: Option<i64>,
    end: Option<i64>,
    gene: Option<String>,
    sample: Option<String>,
    person: Option<String>,
    population: Option<String>,
    limit: Option<i64>,
}

/// 解析区间(显式坐标优先,其次基因名),返回 (chrom, start, end)
fn resolve_region(conn: &rusqlite::Connection, f: &Filter) -> Result<Option<(String, i64, i64)>, String> {
    if let Some(g) = f.gene.as_ref().filter(|g| !g.is_empty()) {
        let row = conn.query_row(
            "SELECT chrom,start,end FROM genes WHERE name=?1",
            params![g],
            |r| Ok((r.get::<_, String>(0)?, r.get::<_, i64>(1)?, r.get::<_, i64>(2)?)),
        );
        return match row {
            Ok(t) => Ok(Some(t)),
            Err(_) => Err(format!("基因 {g} 不存在")),
        };
    }
    if let Some(c) = f.chrom.as_ref().filter(|c| !c.is_empty()) {
        return Ok(Some((
            c.clone(),
            f.start.unwrap_or(1),
            f.end.unwrap_or(i64::MAX / 2),
        )));
    }
    Ok(None)
}

/// 构造 WHERE 片段与参数(区间 + 样本携带条件)
fn build_variant_query(
    conn: &rusqlite::Connection,
    f: &Filter,
) -> Result<(String, Vec<Box<dyn rusqlite::types::ToSql>>), String> {
    let mut wh = String::new();
    let mut args: Vec<Box<dyn rusqlite::types::ToSql>> = Vec::new();
    if let Some((chrom, start, end)) = resolve_region(conn, f)? {
        wh.push_str(" AND v.chrom=? AND v.pos BETWEEN ? AND ?");
        args.push(Box::new(chrom));
        args.push(Box::new(start));
        args.push(Box::new(end));
    }
    let mut sample_conds: Vec<String> = Vec::new();
    if let Some(s) = f.sample.as_ref().filter(|s| !s.is_empty()) {
        sample_conds.push("s.name=?".into());
        args.push(Box::new(s.clone()));
    }
    if let Some(p) = f.person.as_ref().filter(|p| !p.is_empty()) {
        sample_conds.push("s.person_id=?".into());
        args.push(Box::new(p.clone()));
    }
    if let Some(p) = f.population.as_ref().filter(|p| !p.is_empty()) {
        sample_conds.push("s.population=?".into());
        args.push(Box::new(p.clone()));
    }
    if !sample_conds.is_empty() {
        // 位点需被满足条件的样本携带(allele_count>0)
        wh.push_str(&format!(
            " AND EXISTS (SELECT 1 FROM genotypes g JOIN samples s ON s.id=g.sample_id \
             WHERE g.variant_id=v.id AND g.allele_count>0 AND ({}))",
            sample_conds.join(" OR ")
        ));
    }
    Ok((wh, args))
}

// ---------- API ----------

async fn api_stats(State(st): State<Arc<AppState>>) -> Json<Value> {
    let conn = st.conn.lock().unwrap();
    let q = |sql: &str| -> i64 {
        conn.query_row(sql, [], |r| r.get(0)).unwrap_or(0)
    };
    Json(json!({
        "variants": q("SELECT COUNT(*) FROM variants"),
        "samples": q("SELECT COUNT(*) FROM samples"),
        "genotypes": q("SELECT COUNT(*) FROM genotypes"),
        "batches": q("SELECT COUNT(*) FROM batches WHERE status='done'"),
        "qc_failures": q("SELECT COUNT(*) FROM qc_failures"),
        "genes": q("SELECT COUNT(*) FROM genes"),
    }))
}

async fn api_samples(State(st): State<Arc<AppState>>, Query(f): Query<Filter>) -> Json<Value> {
    let conn = st.conn.lock().unwrap();
    let mut sql = "SELECT s.id,s.name,s.person_id,s.population,s.batch_id,b.name,
        (SELECT COUNT(*) FROM genotypes g WHERE g.sample_id=s.id),
        (SELECT COUNT(*) FROM genotypes g WHERE g.sample_id=s.id AND g.allele_count>0)
        FROM samples s JOIN batches b ON b.id=s.batch_id WHERE 1=1".to_string();
    let mut args: Vec<Box<dyn rusqlite::types::ToSql>> = Vec::new();
    if let Some(p) = f.person.as_ref().filter(|p| !p.is_empty()) {
        sql.push_str(" AND s.person_id=?");
        args.push(Box::new(p.clone()));
    }
    if let Some(p) = f.population.as_ref().filter(|p| !p.is_empty()) {
        sql.push_str(" AND s.population=?");
        args.push(Box::new(p.clone()));
    }
    if let Some(s) = f.sample.as_ref().filter(|s| !s.is_empty()) {
        sql.push_str(" AND s.name LIKE ?");
        args.push(Box::new(format!("%{s}%")));
    }
    sql.push_str(" ORDER BY s.id LIMIT 500");
    let mut stmt = conn.prepare(&sql).unwrap();
    let rows = stmt
        .query_map(rusqlite::params_from_iter(args.iter().map(|a| a.as_ref())), |r| {
            Ok(json!({
                "id": r.get::<_, i64>(0)?, "name": r.get::<_, String>(1)?,
                "person_id": r.get::<_, String>(2)?, "population": r.get::<_, String>(3)?,
                "batch_id": r.get::<_, i64>(4)?, "batch": r.get::<_, String>(5)?,
                "n_called": r.get::<_, i64>(6)?, "n_carrier": r.get::<_, i64>(7)?,
            }))
        })
        .unwrap();
    let out: Vec<Value> = rows.map(|r| r.unwrap()).collect();
    let pops: Vec<String> = conn
        .prepare("SELECT DISTINCT population FROM samples ORDER BY 1")
        .unwrap()
        .query_map([], |r| r.get(0))
        .unwrap()
        .map(|r| r.unwrap())
        .collect();
    Json(json!({"samples": out, "populations": pops}))
}

async fn api_duplicates(State(st): State<Arc<AppState>>) -> Json<Value> {
    let conn = st.conn.lock().unwrap();
    // 1) 样本表声明的同人重复(person_id 相同)
    let mut stmt = conn
        .prepare("SELECT person_id, COUNT(*), GROUP_CONCAT(name), GROUP_CONCAT(population)
                  FROM samples GROUP BY person_id HAVING COUNT(*)>1 ORDER BY 2 DESC")
        .unwrap();
    let declared: Vec<Value> = stmt
        .query_map([], |r| {
            Ok(json!({"person_id": r.get::<_,String>(0)?, "n": r.get::<_,i64>(1)?,
                      "samples": r.get::<_,String>(2)?, "populations": r.get::<_,String>(3)?}))
        })
        .unwrap()
        .map(|r| r.unwrap())
        .collect();
    // 2) 基因型指纹检测:所有携带基因型完全一致的样本对
    let mut stmt = conn
        .prepare("SELECT g.sample_id, s.name, s.person_id, g.variant_id, g.allele_count
                  FROM genotypes g JOIN samples s ON s.id=g.sample_id
                  WHERE g.allele_count>0 ORDER BY g.sample_id, g.variant_id")
        .unwrap();
    let rows = stmt
        .query_map([], |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, String>(2)?,
                r.get::<_, i64>(3)?,
                r.get::<_, i64>(4)?,
            ))
        })
        .unwrap();
    let mut fps: HashMap<i64, (String, String, String)> = HashMap::new(); // sample_id -> (name, person, fp)
    for r in rows {
        let (sid, name, person, vid, ac) = r.unwrap();
        let e = fps
            .entry(sid)
            .or_insert_with(|| (name.clone(), person.clone(), String::new()));
        e.2.push_str(&format!("{vid}:{ac};"));
    }
    let mut by_fp: HashMap<&str, Vec<(i64, &str, &str)>> = HashMap::new();
    for (sid, (name, person, fp)) in &fps {
        by_fp.entry(fp.as_str()).or_default().push((*sid, name, person));
    }
    let detected: Vec<Value> = by_fp
        .values()
        .filter(|v| v.len() > 1)
        .map(|v| {
            json!({
                "samples": v.iter().map(|(_,n,_)| n).collect::<Vec<_>>(),
                "persons": v.iter().map(|(_,_,p)| p).collect::<Vec<_>>(),
            })
        })
        .collect();
    // 3) 相似度检测:稀有变异(携带者<=200)上的携带集合重叠系数 >= 0.5
    //    同一人重复送检:两次实验的 QC 过滤是独立噪声,携带集合高度重叠;
    //    无关个体在稀有变异上几乎不重叠。
    let mut stmt = conn
        .prepare("SELECT g.variant_id, g.sample_id FROM genotypes g WHERE g.allele_count>0
                  AND g.variant_id IN (SELECT variant_id FROM genotypes WHERE allele_count>0
                                       GROUP BY variant_id HAVING COUNT(*) BETWEEN 2 AND 200)
                  ORDER BY g.variant_id")
        .unwrap();
    let rows = stmt
        .query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?)))
        .unwrap();
    let mut by_variant: HashMap<i64, Vec<i64>> = HashMap::new();
    let mut totals: HashMap<i64, u32> = HashMap::new();
    for r in rows {
        let (vid, sid) = r.unwrap();
        by_variant.entry(vid).or_default().push(sid);
        *totals.entry(sid).or_default() += 1;
    }
    let mut pair_shared: HashMap<(i64, i64), u32> = HashMap::new();
    for sids in by_variant.values() {
        for i in 0..sids.len() {
            for j in (i + 1)..sids.len() {
                let (a, b) = (sids[i].min(sids[j]), sids[i].max(sids[j]));
                *pair_shared.entry((a, b)).or_default() += 1;
            }
        }
    }
    let mut name_of: HashMap<i64, (String, String)> = HashMap::new();
    {
        let mut stmt = conn.prepare("SELECT id,name,person_id FROM samples").unwrap();
        let rows = stmt
            .query_map([], |r| {
                Ok((r.get::<_, i64>(0)?, r.get::<_, String>(1)?, r.get::<_, String>(2)?))
            })
            .unwrap();
        for r in rows {
            let (id, n, p) = r.unwrap();
            name_of.insert(id, (n, p));
        }
    }
    let mut similar: Vec<Value> = pair_shared
        .iter()
        .filter_map(|((a, b), shared)| {
            let ta = *totals.get(a).unwrap_or(&0);
            let tb = *totals.get(b).unwrap_or(&0);
            let denom = ta.min(tb);
            if *shared < 4 || denom == 0 {
                return None;
            }
            let overlap = *shared as f64 / denom as f64;
            if overlap < 0.4 {
                return None;
            }
            let (na, pa) = name_of.get(a).cloned().unwrap_or_default();
            let (nb, pb) = name_of.get(b).cloned().unwrap_or_default();
            Some(json!({
                "sample_a": na, "sample_b": nb, "person_a": pa, "person_b": pb,
                "shared": shared, "overlap": (overlap * 1000.0).round() / 1000.0,
            }))
        })
        .collect();
    similar.sort_by(|x, y| {
        y["overlap"]
            .as_f64()
            .partial_cmp(&x["overlap"].as_f64())
            .unwrap()
    });
    Json(json!({"declared": declared, "detected": detected, "similar": similar}))
}

async fn api_batches(State(st): State<Arc<AppState>>) -> Json<Value> {
    let conn = st.conn.lock().unwrap();
    let mut stmt = conn
        .prepare("SELECT id,name,created,n_samples,n_records,n_genotypes,n_qc_fail,duration_ms,status FROM batches ORDER BY id")
        .unwrap();
    let rows: Vec<Value> = stmt
        .query_map([], |r| {
            Ok(json!({
                "id": r.get::<_,i64>(0)?, "name": r.get::<_,String>(1)?,
                "created": r.get::<_,String>(2)?, "n_samples": r.get::<_,Option<i64>>(3)?,
                "n_records": r.get::<_,Option<i64>>(4)?, "n_genotypes": r.get::<_,Option<i64>>(5)?,
                "n_qc_fail": r.get::<_,Option<i64>>(6)?, "duration_ms": r.get::<_,Option<i64>>(7)?,
                "status": r.get::<_,String>(8)?,
            }))
        })
        .unwrap()
        .map(|r| r.unwrap())
        .collect();
    Json(json!({"batches": rows}))
}

#[derive(Deserialize)]
struct ImportReq {
    name: String,
    vcf_dir: String,
    sample_sheet: Option<String>,
    min_qual: Option<f64>,
    min_dp: Option<i64>,
    min_gq: Option<f64>,
}

async fn api_import(
    State(st): State<Arc<AppState>>,
    Json(req): Json<ImportReq>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let mut conn = st.conn.lock().unwrap();
    let qc = db::QcThresholds {
        min_qual: req.min_qual.unwrap_or(30.0),
        min_dp: req.min_dp.unwrap_or(10),
        min_gq: req.min_gq.unwrap_or(20.0),
    };
    match db::import_batch(
        &mut conn,
        &st.reference,
        &req.name,
        &req.vcf_dir,
        req.sample_sheet.as_deref(),
        &qc,
    ) {
        Ok(r) => Ok(Json(json!({
            "ok": true, "batch_id": r.batch_id, "n_samples": r.n_samples,
            "n_records": r.n_records, "n_new_variants": r.n_new_variants,
            "n_genotypes": r.n_genotypes, "n_qc_fail": r.n_qc_fail,
            "duration_ms": r.duration_ms,
            "records_per_sec": (r.n_records as f64 / (r.duration_ms as f64 / 1000.0)) as i64,
        }))),
        Err(e) => Err((StatusCode::BAD_REQUEST, Json(json!({"ok": false, "error": e})))),
    }
}

async fn api_variants(
    State(st): State<Arc<AppState>>,
    Query(f): Query<Filter>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let t0 = Instant::now();
    let conn = st.conn.lock().unwrap();
    let (wh, args) = build_variant_query(&conn, &f)
        .map_err(|e| (StatusCode::BAD_REQUEST, Json(json!({"error": e}))))?;
    let limit = f.limit.unwrap_or(500).min(5000);
    let sql = format!(
        "SELECT v.id,v.chrom,v.pos,v.ref,v.alt,
                COALESCE(SUM(g.allele_count),0), COUNT(g.sample_id)*2
         FROM variants v LEFT JOIN genotypes g ON g.variant_id=v.id
         WHERE 1=1 {wh} GROUP BY v.id ORDER BY v.chrom,v.pos LIMIT {limit}"
    );
    let mut stmt = conn.prepare(&sql).unwrap();
    let rows = stmt
        .query_map(rusqlite::params_from_iter(args.iter().map(|a| a.as_ref())), |r| {
            let ac: i64 = r.get(5)?;
            let an: i64 = r.get(6)?;
            Ok(json!({
                "id": r.get::<_,i64>(0)?, "chrom": r.get::<_,String>(1)?,
                "pos": r.get::<_,i64>(2)?, "ref": r.get::<_,String>(3)?,
                "alt": r.get::<_,String>(4)?, "ac": ac, "an": an,
                "af": if an>0 { ac as f64 / an as f64 } else { 0.0 },
            }))
        })
        .unwrap();
    let out: Vec<Value> = rows.map(|r| r.unwrap()).collect();
    let elapsed = t0.elapsed().as_secs_f64() * 1000.0;
    Ok(Json(json!({"variants": out, "count": out.len(), "query_ms": elapsed})))
}

async fn api_frequencies(
    State(st): State<Arc<AppState>>,
    Query(f): Query<Filter>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let t0 = Instant::now();
    let conn = st.conn.lock().unwrap();
    let (wh, args) = build_variant_query(&conn, &f)
        .map_err(|e| (StatusCode::BAD_REQUEST, Json(json!({"error": e}))))?;
    let limit = f.limit.unwrap_or(200).min(2000);
    // 总频率
    let sql = format!(
        "SELECT v.id,v.chrom,v.pos,v.ref,v.alt,COALESCE(SUM(g.allele_count),0),COUNT(g.sample_id)*2
         FROM variants v LEFT JOIN genotypes g ON g.variant_id=v.id
         WHERE 1=1 {wh} GROUP BY v.id ORDER BY v.chrom,v.pos LIMIT {limit}"
    );
    let mut stmt = conn.prepare(&sql).unwrap();
    let pref = rusqlite::params_from_iter(args.iter().map(|a| a.as_ref()));
    let rows = stmt
        .query_map(pref, |r| {
            Ok((
                r.get::<_, i64>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, i64>(2)?,
                r.get::<_, String>(3)?,
                r.get::<_, String>(4)?,
                r.get::<_, i64>(5)?,
                r.get::<_, i64>(6)?,
            ))
        })
        .unwrap();
    let mut variants: Vec<Value> = Vec::new();
    let mut ids: Vec<i64> = Vec::new();
    for r in rows {
        let (id, chrom, pos, rf, alt, ac, an) = r.unwrap();
        ids.push(id);
        variants.push(json!({
            "id": id, "chrom": chrom, "pos": pos, "ref": rf, "alt": alt,
            "ac": ac, "an": an, "af": if an>0 {ac as f64/an as f64} else {0.0},
            "populations": {},
        }));
    }
    // 分人群频率
    if !ids.is_empty() {
        let placeholders = ids.iter().map(|_| "?").collect::<Vec<_>>().join(",");
        let sql = format!(
            "SELECT g.variant_id, s.population, SUM(g.allele_count), COUNT(*)*2
             FROM genotypes g JOIN samples s ON s.id=g.sample_id
             WHERE g.variant_id IN ({placeholders}) GROUP BY g.variant_id, s.population"
        );
        let mut stmt = conn.prepare(&sql).unwrap();
        let rows = stmt
            .query_map(rusqlite::params_from_iter(ids.iter()), |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, i64>(2)?,
                    r.get::<_, i64>(3)?,
                ))
            })
            .unwrap();
        let mut pmap: HashMap<i64, HashMap<String, (i64, i64)>> = HashMap::new();
        for r in rows {
            let (vid, pop, ac, an) = r.unwrap();
            pmap.entry(vid).or_default().insert(pop, (ac, an));
        }
        for v in variants.iter_mut() {
            let id = v["id"].as_i64().unwrap();
            if let Some(m) = pmap.get(&id) {
                let mut obj = serde_json::Map::new();
                for (pop, (ac, an)) in m {
                    obj.insert(
                        pop.clone(),
                        json!({"ac": ac, "an": an, "af": if *an>0 {*ac as f64/ *an as f64} else {0.0}}),
                    );
                }
                v["populations"] = Value::Object(obj);
            }
        }
    }
    let elapsed = t0.elapsed().as_secs_f64() * 1000.0;
    Ok(Json(json!({"variants": variants, "count": variants.len(), "query_ms": elapsed})))
}

async fn api_qc(State(st): State<Arc<AppState>>, Query(q): Query<HashMap<String, String>>) -> Json<Value> {
    let conn = st.conn.lock().unwrap();
    let batch_id = q.get("batch_id").and_then(|s| s.parse::<i64>().ok());
    let (wh, args): (String, Vec<i64>) = match batch_id {
        Some(b) => ("WHERE batch_id=?".into(), vec![b]),
        None => (String::new(), vec![]),
    };
    let mut stmt = conn
        .prepare(&format!(
            "SELECT reason, COUNT(*) FROM qc_failures {wh} GROUP BY reason ORDER BY 2 DESC"
        ))
        .unwrap();
    let summary: Vec<Value> = stmt
        .query_map(rusqlite::params_from_iter(args.iter()), |r| {
            Ok(json!({"reason": r.get::<_,String>(0)?, "count": r.get::<_,i64>(1)?}))
        })
        .unwrap()
        .map(|r| r.unwrap())
        .collect();
    let mut stmt = conn
        .prepare(&format!(
            "SELECT batch_id,sample_name,chrom,pos,ref,alt,reason,detail FROM qc_failures {wh} ORDER BY id DESC LIMIT 500"
        ))
        .unwrap();
    let rows: Vec<Value> = stmt
        .query_map(rusqlite::params_from_iter(args.iter()), |r| {
            Ok(json!({
                "batch_id": r.get::<_,i64>(0)?, "sample": r.get::<_,String>(1)?,
                "chrom": r.get::<_,String>(2)?, "pos": r.get::<_,i64>(3)?,
                "ref": r.get::<_,String>(4)?, "alt": r.get::<_,String>(5)?,
                "reason": r.get::<_,String>(6)?, "detail": r.get::<_,String>(7)?,
            }))
        })
        .unwrap()
        .map(|r| r.unwrap())
        .collect();
    Json(json!({"summary": summary, "records": rows}))
}

async fn api_export_vcf(
    State(st): State<Arc<AppState>>,
    Query(f): Query<Filter>,
) -> Result<Response, (StatusCode, Json<Value>)> {
    let conn = st.conn.lock().unwrap();
    let (wh, args) = build_variant_query(&conn, &f)
        .map_err(|e| (StatusCode::BAD_REQUEST, Json(json!({"error": e}))))?;
    let limit = f.limit.unwrap_or(10000).min(100000);
    // 样本集合(按条件过滤,否则全部)
    let mut ssql = "SELECT id,name FROM samples".to_string();
    let mut sargs: Vec<Box<dyn rusqlite::types::ToSql>> = Vec::new();
    let mut conds: Vec<String> = Vec::new();
    if let Some(s) = f.sample.as_ref().filter(|s| !s.is_empty()) {
        conds.push("name=?".into());
        sargs.push(Box::new(s.clone()));
    }
    if let Some(p) = f.person.as_ref().filter(|p| !p.is_empty()) {
        conds.push("person_id=?".into());
        sargs.push(Box::new(p.clone()));
    }
    if let Some(p) = f.population.as_ref().filter(|p| !p.is_empty()) {
        conds.push("population=?".into());
        sargs.push(Box::new(p.clone()));
    }
    if !conds.is_empty() {
        ssql.push_str(&format!(" WHERE {}", conds.join(" OR ")));
    }
    ssql.push_str(" ORDER BY id");
    let mut stmt = conn.prepare(&ssql).unwrap();
    let samples: Vec<(i64, String)> = stmt
        .query_map(rusqlite::params_from_iter(sargs.iter().map(|a| a.as_ref())), |r| {
            Ok((r.get::<_, i64>(0)?, r.get::<_, String>(1)?))
        })
        .unwrap()
        .map(|r| r.unwrap())
        .collect();

    let sql = format!(
        "SELECT v.id,v.chrom,v.pos,v.ref,v.alt FROM variants v WHERE 1=1 {wh}
         ORDER BY v.chrom,v.pos LIMIT {limit}"
    );
    let mut stmt = conn.prepare(&sql).unwrap();
    let variants: Vec<(i64, String, i64, String, String)> = stmt
        .query_map(rusqlite::params_from_iter(args.iter().map(|a| a.as_ref())), |r| {
            Ok((
                r.get(0)?,
                r.get(1)?,
                r.get(2)?,
                r.get(3)?,
                r.get(4)?,
            ))
        })
        .unwrap()
        .map(|r| r.unwrap())
        .collect();

    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.2\n##source=variantdb\n");
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    out.push_str("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (_, name) in &samples {
        out.push('\t');
        out.push_str(name);
    }
    out.push('\n');
    let mut gstmt = conn
        .prepare("SELECT sample_id,gt,dp,gq FROM genotypes WHERE variant_id=?1")
        .unwrap();
    for (vid, chrom, pos, rf, alt) in &variants {
        out.push_str(&format!("{chrom}\t{pos}\t.\t{rf}\t{alt}\t.\tPASS\t.\tGT:DP:GQ"));
        let gmap: HashMap<i64, (String, Option<i64>, Option<f64>)> = gstmt
            .query_map(params![vid], |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, Option<i64>>(2)?,
                    r.get::<_, Option<f64>>(3)?,
                ))
            })
            .unwrap()
            .map(|r| r.unwrap())
            .map(|(sid, gt, dp, gq)| (sid, (gt, dp, gq)))
            .collect();
        for (sid, _) in &samples {
            match gmap.get(sid) {
                Some((gt, dp, gq)) => out.push_str(&format!(
                    "\t{}:{}:{}",
                    gt,
                    dp.map(|d| d.to_string()).unwrap_or(".".into()),
                    gq.map(|g| g.to_string()).unwrap_or(".".into())
                )),
                None => out.push_str("\t./.:.:."),
            }
        }
        out.push('\n');
    }
    Ok((
        [(header::CONTENT_TYPE, "text/plain; charset=utf-8")],
        out,
    )
        .into_response())
}

/// 基准测试:随机区间查询延迟
async fn api_bench(State(st): State<Arc<AppState>>, Query(q): Query<HashMap<String, String>>) -> Json<Value> {
    let n: usize = q.get("n").and_then(|s| s.parse().ok()).unwrap_or(200);
    let conn = st.conn.lock().unwrap();
    let max_pos: i64 = conn
        .query_row("SELECT COALESCE(MAX(pos),1) FROM variants", [], |r| r.get(0))
        .unwrap_or(1);
    let chrom: String = conn
        .query_row("SELECT chrom FROM variants LIMIT 1", [], |r| r.get(0))
        .unwrap_or_else(|_| "chr1".into());
    let mut lat: Vec<f64> = Vec::with_capacity(n);
    // 简单 LCG 随机数,避免额外依赖
    let mut seed: u64 = 0x9e3779b97f4a7c15;
    let mut next = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) as i64
    };
    let mut stmt = conn
        .prepare("SELECT COUNT(*) FROM variants v WHERE v.chrom=?1 AND v.pos BETWEEN ?2 AND ?3")
        .unwrap();
    let mut stmt2 = conn
        .prepare("SELECT v.id, COALESCE(SUM(g.allele_count),0), COUNT(g.sample_id)*2
                  FROM variants v LEFT JOIN genotypes g ON g.variant_id=v.id
                  WHERE v.chrom=?1 AND v.pos BETWEEN ?2 AND ?3 GROUP BY v.id")
        .unwrap();
    for _ in 0..n {
        let width = 1 + next() % 5000;
        let start = 1 + next() % (max_pos.max(2) - 1);
        let end = start + width;
        let t0 = Instant::now();
        let _: i64 = stmt
            .query_row(params![chrom, start, end], |r| r.get(0))
            .unwrap_or(0);
        let rows = stmt2
            .query_map(params![chrom, start, end], |r| {
                Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?, r.get::<_, i64>(2)?))
            })
            .unwrap();
        let _: Vec<_> = rows.collect();
        lat.push(t0.elapsed().as_secs_f64() * 1000.0);
    }
    lat.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let avg = lat.iter().sum::<f64>() / lat.len() as f64;
    let p = |q: f64| lat[((lat.len() as f64 - 1.0) * q) as usize];
    Json(json!({
        "n_queries": n, "avg_ms": avg, "p50_ms": p(0.5), "p95_ms": p(0.95),
        "p99_ms": p(0.99), "max_ms": lat.last(),
    }))
}
