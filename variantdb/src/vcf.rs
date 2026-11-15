// VCF 解析与变异标准化(多等位拆分 + 左对齐)

/// 参考基因组:染色体名 -> 序列(大写)
pub struct Reference {
    pub seqs: std::collections::HashMap<String, Vec<u8>>,
}

impl Reference {
    pub fn from_fasta(path: &str) -> Result<Reference, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("read ref: {e}"))?;
        let mut seqs = std::collections::HashMap::new();
        let mut name: Option<String> = None;
        let mut cur: Vec<u8> = Vec::new();
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('>') {
                if let Some(n) = name.take() {
                    seqs.insert(n, std::mem::take(&mut cur));
                }
                name = Some(h.split_whitespace().next().unwrap_or("").to_string());
            } else {
                cur.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        if let Some(n) = name {
            seqs.insert(n, cur);
        }
        if seqs.is_empty() {
            return Err("empty fasta".into());
        }
        Ok(Reference { seqs })
    }
}

/// 对单个等位做最小化 + 左对齐标准化。
/// pos 为 1-based(VCF 约定),返回 (pos, ref, alt)。
pub fn normalize_allele(
    refseq: &[u8],
    mut pos: usize,
    mut r: Vec<u8>,
    mut a: Vec<u8>,
) -> (usize, String, String) {
    // 不动点循环:右修剪 -> 左修剪 -> 空等位补锚碱基(等价于左移)
    loop {
        if r == a {
            break; // 改写后与参考完全等同(无实际变异),无法继续规范化
        }
        let mut changed = false;
        // 1. 修剪共同后缀
        while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
            r.pop();
            a.pop();
            changed = true;
        }
        // 2. 修剪共同前缀(各保留至少 1 个碱基)
        while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
            r.remove(0);
            a.remove(0);
            pos += 1;
            changed = true;
        }
        // 3. 某一等位为空(纯插入/缺失): 用参考序列左侧碱基补锚,位置左移
        if r.is_empty() || a.is_empty() {
            if pos > 1 {
                let b = refseq[pos - 2]; // pos 1-based -> 前一个碱基的 0-based 下标
                r.insert(0, b);
                a.insert(0, b);
                pos -= 1;
                changed = true;
            } else {
                // 已在染色体最左端,用右侧碱基补锚
                let idx = pos - 1 + r.len().max(a.len());
                let b = if idx < refseq.len() { refseq[idx] } else { b'N' };
                r.push(b);
                a.push(b);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    (
        pos,
        String::from_utf8_lossy(&r).into_owned(),
        String::from_utf8_lossy(&a).into_owned(),
    )
}

#[derive(Debug, Clone)]
pub struct NormVariant {
    pub chrom: String,
    pub pos: usize, // 1-based, 标准化后
    pub ref_: String,
    pub alt: String,
    pub orig_allele_idx: usize, // 原始 ALT 列表中的序号(1-based)
}

#[derive(Debug, Clone)]
pub struct SampleGenotype {
    pub sample: String,
    pub allele_count: u32, // 目标等位的拷贝数 0/1/2
    pub called: bool,      // GT 是否可用(非 ./.)
    pub gt_norm: String,   // 相对该等位的基因型 0/0 0/1 1/1
    pub dp: Option<i64>,
    pub gq: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ParsedRecord {
    pub variant: NormVariant,
    pub qual: Option<f64>,
    pub genotypes: Vec<SampleGenotype>,
}

/// 解析单个 VCF 文件,返回 (样本名列表, 拆分+标准化后的记录列表)
pub fn parse_vcf(path: &str, reference: &Reference) -> Result<(Vec<String>, Vec<ParsedRecord>), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    parse_vcf_text(&text, reference)
}

pub fn parse_vcf_text(
    text: &str,
    reference: &Reference,
) -> Result<(Vec<String>, Vec<ParsedRecord>), String> {
    let mut samples: Vec<String> = Vec::new();
    let mut records: Vec<ParsedRecord> = Vec::new();

    for (lineno, line) in text.lines().enumerate() {
        if line.starts_with("##") {
            continue;
        }
        if let Some(h) = line.strip_prefix("#CHROM") {
            let cols: Vec<&str> = h.split('\t').collect();
            if cols.len() > 9 {
                samples = cols[9..].iter().map(|s| s.to_string()).collect();
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            return Err(format!("line {}: too few columns", lineno + 1));
        }
        let chrom_raw = cols[0].to_string();
        // 染色体名归一:统一去/加 "chr" 前缀,与参考基因组保持一致
        let chrom = chrom_raw.clone();
        let refseq = reference
            .seqs
            .get(&chrom)
            .ok_or_else(|| format!("line {}: unknown chrom {chrom}", lineno + 1))?;
        let pos: usize = cols[1]
            .parse()
            .map_err(|_| format!("line {}: bad POS", lineno + 1))?;
        if pos < 1 || pos > refseq.len() {
            return Err(format!("line {}: POS {pos} out of range", lineno + 1));
        }
        let ref_allele = cols[3].to_uppercase();
        // 与参考序列比对(宽容:不一致时截断到参考长度继续,保证导入健壮)
        let qual: Option<f64> = match cols[5] {
            "." => None,
            q => q.parse().ok(),
        };
        let format_keys: Vec<&str> = if cols.len() > 8 { cols[8].split(':').collect() } else { Vec::new() };
        let gt_idx = format_keys.iter().position(|k| *k == "GT");
        let dp_idx = format_keys.iter().position(|k| *k == "DP");
        let gq_idx = format_keys.iter().position(|k| *k == "GQ");

        let alts: Vec<&str> = cols[4].split(',').collect();
        for (ai, alt_raw) in alts.iter().enumerate() {
            let allele_idx = ai + 1;
            let alt_up = alt_raw.to_uppercase();
            if alt_up.contains('<') || alt_up == "." || alt_up == "*" {
                continue; // 跳过符号等位/缺失
            }
            let (npos, nref, nalt) = normalize_allele(
                refseq,
                pos,
                ref_allele.as_bytes().to_vec(),
                alt_up.as_bytes().to_vec(),
            );
            let mut gts: Vec<SampleGenotype> = Vec::with_capacity(samples.len());
            for (si, sname) in samples.iter().enumerate() {
                let col = cols.get(9 + si).copied().unwrap_or(".");
                let fields: Vec<&str> = col.split(':').collect();
                let getf = |idx: Option<usize>| -> Option<&str> {
                    idx.and_then(|i| fields.get(i).copied()).filter(|v| *v != ".")
                };
                let dp: Option<i64> = getf(dp_idx).and_then(|v| v.parse().ok());
                let gq: Option<f64> = getf(gq_idx).and_then(|v| v.parse().ok());
                let gt_raw = getf(gt_idx).unwrap_or("./.");
                let mut count: u32 = 0;
                let mut called = true;
                let mut ploidy = 0u32;
                for a in gt_raw.split(|c| c == '/' || c == '|') {
                    if a == "." {
                        called = false;
                    } else if let Ok(n) = a.parse::<usize>() {
                        ploidy += 1;
                        if n == allele_idx {
                            count += 1;
                        }
                    }
                }
                if ploidy == 0 {
                    called = false;
                }
                let gt_norm = if !called {
                    "./.".to_string()
                } else {
                    match count {
                        0 => "0/0".to_string(),
                        1 => "0/1".to_string(),
                        _ => "1/1".to_string(),
                    }
                };
                gts.push(SampleGenotype {
                    sample: sname.clone(),
                    allele_count: count,
                    called,
                    gt_norm,
                    dp,
                    gq,
                });
            }
            records.push(ParsedRecord {
                variant: NormVariant {
                    chrom: chrom.clone(),
                    pos: npos,
                    ref_: nref,
                    alt: nalt,
                    orig_allele_idx: allele_idx,
                },
                qual,
                genotypes: gts,
            });
        }
    }
    if samples.is_empty() {
        return Err("no samples in VCF header".into());
    }
    Ok((samples, records))
}
