// SQLite 存储层:schema、事务化批量导入、查询
use rusqlite::{params, Connection};
use std::time::Instant;

use crate::vcf::{self, Reference};

pub struct QcThresholds {
    pub min_qual: f64,
    pub min_dp: i64,
    pub min_gq: f64,
}

pub struct ImportResult {
    pub batch_id: i64,
    pub n_samples: i64,
    pub n_records: i64,
    pub n_new_variants: i64,
    pub n_genotypes: i64,
    pub n_qc_fail: i64,
    pub duration_ms: i64,
}

pub fn init_schema(conn: &Connection) -> rusqlite::Result<()> {
    conn.execute_batch(
        "
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE IF NOT EXISTS batches(
            id INTEGER PRIMARY KEY,
            name TEXT UNIQUE NOT NULL,
            created TEXT DEFAULT (datetime('now')),
            n_samples INTEGER, n_records INTEGER, n_genotypes INTEGER,
            n_qc_fail INTEGER, duration_ms INTEGER,
            status TEXT NOT NULL DEFAULT 'importing'
        );
        CREATE TABLE IF NOT EXISTS samples(
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL,
            person_id TEXT NOT NULL,
            population TEXT NOT NULL,
            batch_id INTEGER NOT NULL REFERENCES batches(id),
            UNIQUE(name, batch_id)
        );
        CREATE TABLE IF NOT EXISTS variants(
            id INTEGER PRIMARY KEY,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            UNIQUE(chrom, pos, ref, alt)
        );
        CREATE TABLE IF NOT EXISTS genotypes(
            variant_id INTEGER NOT NULL REFERENCES variants(id),
            sample_id INTEGER NOT NULL REFERENCES samples(id),
            allele_count INTEGER NOT NULL,
            gt TEXT NOT NULL,
            dp INTEGER, gq REAL,
            PRIMARY KEY(variant_id, sample_id)
        );
        CREATE TABLE IF NOT EXISTS qc_failures(
            id INTEGER PRIMARY KEY,
            batch_id INTEGER NOT NULL,
            sample_name TEXT,
            chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
            reason TEXT NOT NULL,
            detail TEXT
        );
        CREATE TABLE IF NOT EXISTS genes(
            name TEXT PRIMARY KEY,
            chrom TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_variants_region ON variants(chrom, pos);
        CREATE INDEX IF NOT EXISTS idx_gt_sample ON genotypes(sample_id);
        CREATE INDEX IF NOT EXISTS idx_gt_variant ON genotypes(variant_id);
        CREATE INDEX IF NOT EXISTS idx_qc_batch ON qc_failures(batch_id);
        CREATE INDEX IF NOT EXISTS idx_samples_person ON samples(person_id);
        ",
    )?;
    Ok(())
}

pub fn load_genes(conn: &Connection, path: &str) -> Result<usize, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("read genes: {e}"))?;
    let mut n = 0;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 4 {
            continue;
        }
        conn.execute(
            "INSERT OR REPLACE INTO genes(name,chrom,start,end) VALUES(?1,?2,?3,?4)",
            params![c[0], c[1], c[2].parse::<i64>().unwrap_or(0), c[3].parse::<i64>().unwrap_or(0)],
        )
        .map_err(|e| e.to_string())?;
        n += 1;
    }
    Ok(n)
}

/// 样本表 CSV: name,person_id,population (带表头)
fn load_sample_sheet(path: &str) -> Result<std::collections::HashMap<String, (String, String)>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("read sample sheet: {e}"))?;
    let mut map = std::collections::HashMap::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || (i == 0 && line.to_ascii_lowercase().starts_with("name")) {
            continue;
        }
        let c: Vec<&str> = line.split(',').collect();
        if c.len() >= 3 {
            map.insert(c[0].to_string(), (c[1].to_string(), c[2].to_string()));
        }
    }
    Ok(map)
}

/// 事务化导入一个批次(目录下全部 .vcf)。任何错误 -> 整体回滚。
pub fn import_batch(
    conn: &mut Connection,
    reference: &Reference,
    name: &str,
    vcf_dir: &str,
    sample_sheet: Option<&str>,
    qc: &QcThresholds,
) -> Result<ImportResult, String> {
    let t0 = Instant::now();
    let sheet = match sample_sheet {
        Some(p) if !p.is_empty() => load_sample_sheet(p)?,
        _ => std::collections::HashMap::new(),
    };
    let mut files: Vec<_> = std::fs::read_dir(vcf_dir)
        .map_err(|e| format!("read dir {vcf_dir}: {e}"))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|x| x == "vcf").unwrap_or(false))
        .collect();
    files.sort();
    if files.is_empty() {
        return Err(format!("no .vcf files in {vcf_dir}"));
    }

    let tx = conn.transaction().map_err(|e| e.to_string())?;
    let result = (|| -> Result<ImportResult, String> {
        tx.execute(
            "INSERT INTO batches(name,status) VALUES(?1,'importing')",
            params![name],
        )
        .map_err(|e| format!("insert batch (name 重复?): {e}"))?;
        let batch_id = tx.last_insert_rowid();

        let mut n_samples = 0i64;
        let mut n_records = 0i64;
        let mut n_genotypes = 0i64;
        let mut n_qc_fail = 0i64;
        let mut new_variants = 0i64;

        {
            let mut ins_sample = tx
                .prepare("INSERT INTO samples(name,person_id,population,batch_id) VALUES(?1,?2,?3,?4)")
                .map_err(|e| e.to_string())?;
            let mut get_variant = tx
                .prepare("SELECT id FROM variants WHERE chrom=?1 AND pos=?2 AND ref=?3 AND alt=?4")
                .map_err(|e| e.to_string())?;
            let mut ins_variant = tx
                .prepare("INSERT INTO variants(chrom,pos,ref,alt) VALUES(?1,?2,?3,?4)")
                .map_err(|e| e.to_string())?;
            let mut ins_gt = tx
                .prepare("INSERT OR IGNORE INTO genotypes(variant_id,sample_id,allele_count,gt,dp,gq) VALUES(?1,?2,?3,?4,?5,?6)")
                .map_err(|e| e.to_string())?;
            let mut ins_qc = tx
                .prepare("INSERT INTO qc_failures(batch_id,sample_name,chrom,pos,ref,alt,reason,detail) VALUES(?1,?2,?3,?4,?5,?6,?7,?8)")
                .map_err(|e| e.to_string())?;

            for f in &files {
                let fpath = f.to_string_lossy().to_string();
                let (vcf_samples, records) = vcf::parse_vcf(&fpath, reference)
                    .map_err(|e| format!("{fpath}: {e}"))?;
                // 建样本
                let mut sample_ids: Vec<i64> = Vec::with_capacity(vcf_samples.len());
                for s in &vcf_samples {
                    let (person, pop) = sheet
                        .get(s)
                        .cloned()
                        .unwrap_or_else(|| (s.clone(), "UNKNOWN".to_string()));
                    ins_sample
                        .execute(params![s, person, pop, batch_id])
                        .map_err(|e| format!("insert sample {s}: {e}"))?;
                    sample_ids.push(tx.last_insert_rowid());
                    n_samples += 1;
                }
                for rec in &records {
                    n_records += 1;
                    let v = &rec.variant;
                    // 位点质控:QUAL 不达标 -> 整行记录失败原因,不入库
                    if let Some(q) = rec.qual {
                        if q < qc.min_qual {
                            for s in &vcf_samples {
                                ins_qc
                                    .execute(params![
                                        batch_id, s, v.chrom, v.pos as i64, v.ref_, v.alt,
                                        "LOW_QUAL", format!("QUAL={q} < {}", qc.min_qual)
                                    ])
                                    .map_err(|e| e.to_string())?;
                                n_qc_fail += 1;
                            }
                            continue;
                        }
                    }
                    // 变异合并:标准化后 (chrom,pos,ref,alt) 相同即同一位点
                    let vid: i64 = match get_variant.query_row(
                        params![v.chrom, v.pos as i64, v.ref_, v.alt],
                        |r| r.get(0),
                    ) {
                        Ok(id) => id,
                        Err(rusqlite::Error::QueryReturnedNoRows) => {
                            ins_variant
                                .execute(params![v.chrom, v.pos as i64, v.ref_, v.alt])
                                .map_err(|e| e.to_string())?;
                            new_variants += 1;
                            tx.last_insert_rowid()
                        }
                        Err(e) => return Err(e.to_string()),
                    };
                    for (si, g) in rec.genotypes.iter().enumerate() {
                        if !g.called {
                            continue; // 无呼叫,不计失败也不算 AN
                        }
                        if let Some(dp) = g.dp {
                            if dp < qc.min_dp {
                                ins_qc
                                    .execute(params![
                                        batch_id, g.sample, v.chrom, v.pos as i64, v.ref_, v.alt,
                                        "LOW_DP", format!("DP={dp} < {}", qc.min_dp)
                                    ])
                                    .map_err(|e| e.to_string())?;
                                n_qc_fail += 1;
                                continue;
                            }
                        }
                        if let Some(gq) = g.gq {
                            if gq < qc.min_gq {
                                ins_qc
                                    .execute(params![
                                        batch_id, g.sample, v.chrom, v.pos as i64, v.ref_, v.alt,
                                        "LOW_GQ", format!("GQ={gq} < {}", qc.min_gq)
                                    ])
                                    .map_err(|e| e.to_string())?;
                                n_qc_fail += 1;
                                continue;
                            }
                        }
                        ins_gt
                            .execute(params![
                                vid,
                                sample_ids[si],
                                g.allele_count,
                                g.gt_norm,
                                g.dp,
                                g.gq
                            ])
                            .map_err(|e| e.to_string())?;
                        n_genotypes += 1;
                    }
                }
            }
        }

        let duration_ms = t0.elapsed().as_millis() as i64;
        tx.execute(
            "UPDATE batches SET n_samples=?1,n_records=?2,n_genotypes=?3,n_qc_fail=?4,duration_ms=?5,status='done' WHERE id=?6",
            params![n_samples, n_records, n_genotypes, n_qc_fail, duration_ms, batch_id],
        )
        .map_err(|e| e.to_string())?;
        Ok(ImportResult {
            batch_id,
            n_samples,
            n_records,
            n_new_variants: new_variants,
            n_genotypes,
            n_qc_fail,
            duration_ms,
        })
    })();
    match result {
        Ok(r) => {
            tx.commit().map_err(|e| format!("commit: {e}"))?;
            Ok(r)
        }
        Err(e) => {
            let _ = tx.rollback(); // 整体回滚,不留半截数据
            Err(format!("批次 {name} 导入失败,已整体回滚: {e}"))
        }
    }
}
