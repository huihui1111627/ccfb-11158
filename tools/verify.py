#!/usr/bin/env python3
"""暴力校验:独立实现的 VCF 规范化 + 全文件扫描,与 variantdb API 查询结果对比。
用法: python3 -I verify.py [data_dir] [api_base]
"""
import os, sys, csv, json, time, random, urllib.request, urllib.parse, glob
from collections import defaultdict

DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
API = sys.argv[2] if len(sys.argv) > 2 else "http://localhost:8000"
MIN_QUAL, MIN_DP, MIN_GQ = 30.0, 10, 20.0
random.seed(7)

# ---------- 参考基因组 ----------
ref = {}
name = None
for line in open(f"{DATA}/ref.fa"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        ref[name] = []
    else:
        ref[name].extend(line)
ref = {k: "".join(v) for k, v in ref.items()}

genes = {}
for line in open(f"{DATA}/genes.tsv"):
    c = line.strip().split("\t")
    if len(c) >= 4:
        genes[c[0]] = (c[1], int(c[2]), int(c[3]))

# 样本表
sample_meta = {}  # name -> (person, pop)
for f in glob.glob(f"{DATA}/samples_batch*.csv"):
    for row in csv.DictReader(open(f)):
        sample_meta[row["name"]] = (row["person_id"], row["population"])

# ---------- 独立实现的左对齐标准化 ----------
def normalize(seq, pos, r, a):
    """pos 1-based;不动点:右修剪/左修剪/补锚"""
    r, a = list(r), list(a)
    while True:
        if r == a:
            break  # 与参考等同,无法继续规范化
        changed = False
        while r and a and r[-1] == a[-1]:
            r.pop(); a.pop(); changed = True
        while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
            r.pop(0); a.pop(0); pos += 1; changed = True
        if not r or not a:
            if pos > 1:
                b = seq[pos - 2]
                r.insert(0, b); a.insert(0, b); pos -= 1; changed = True
            else:
                idx = pos - 1 + max(len(r), len(a))
                b = seq[idx] if idx < len(seq) else "N"
                r.append(b); a.append(b); changed = True
        if not changed:
            break
    return pos, "".join(r), "".join(a)

# ---------- 暴力扫描全部 VCF ----------
ac = defaultdict(int)          # variant -> allele count
an = defaultdict(int)          # variant -> allele number (2*called)
seen = set()                   # 通过 QUAL 过滤的变异(对应库中 variants 表)
carriers = defaultdict(set)    # variant -> {sample}
pop_ac = defaultdict(lambda: defaultdict(int))   # variant -> pop -> ac
pop_an = defaultdict(lambda: defaultdict(int))
sample_vars = defaultdict(set) # sample -> {variant carried}
t0 = time.time()
n_files = 0
for f in sorted(glob.glob(f"{DATA}/batch*/*.vcf")):
    n_files += 1
    samples = []
    for line in open(f):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            continue
        c = line.rstrip("\n").split("\t")
        chrom, pos, r, alts, qual = c[0], int(c[1]), c[3].upper(), c[4].split(","), c[5]
        qual = None if qual == "." else float(qual)
        if qual is not None and qual < MIN_QUAL:
            continue
        fmt = c[8].split(":")
        gi, di, qi = fmt.index("GT"), fmt.index("DP"), fmt.index("GQ")
        seq = ref[chrom]
        for ai, alt in enumerate(alts, start=1):
            alt = alt.upper()
            if "<" in alt or alt in (".", "*"):
                continue
            np_, nr, na = normalize(seq, pos, r, alt)
            v = (chrom, np_, nr, na)
            seen.add(v)  # 库语义:QUAL 通过即入 variants 表,DP/GQ 只影响基因型
            for si, sname in enumerate(samples):
                fields = c[9 + si].split(":")
                gt = fields[gi]
                if "." in gt:
                    continue
                dp = int(fields[di]) if di < len(fields) and fields[di] != "." else None
                gq = float(fields[qi]) if qi < len(fields) and fields[qi] != "." else None
                if dp is not None and dp < MIN_DP:
                    continue
                if gq is not None and gq < MIN_GQ:
                    continue
                cnt = sum(1 for x in gt.replace("|", "/").split("/") if x == str(ai))
                ac[v] += cnt
                an[v] += 2
                pop = sample_meta.get(sname, (sname, "UNKNOWN"))[1]
                pop_ac[v][pop] += cnt
                pop_an[v][pop] += 2
                if cnt > 0:
                    carriers[v].add(sname)
                    sample_vars[sname].add(v)
scan_sec = time.time() - t0
print(f"[brute] scanned {n_files} VCFs in {scan_sec:.1f}s -> {len(ac)} distinct normalized variants")

def api(path, params=None):
    url = API + path + ("?" + urllib.parse.urlencode(params) if params else "")
    t0 = time.time()
    with urllib.request.urlopen(url) as r:
        d = json.load(r)
    return d, (time.time() - t0) * 1000

fails = 0
def check(label, ok, extra=""):
    global fails
    if not ok:
        fails += 1
        print(f"  FAIL {label} {extra}")

# ---------- 1. 随机区间查询对比 ----------
print("[1] interval queries")
lat = []
all_vars = sorted(seen)
for i in range(60):
    center = random.choice(all_vars)[1]
    start = max(1, center - random.randint(0, 3000))
    end = center + random.randint(0, 3000)
    d, ms = api("/api/variants", {"chrom": "chr1", "start": start, "end": end, "limit": 5000})
    lat.append(ms)
    got = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v for v in d["variants"]}
    want = {v for v in seen if v[0] == "chr1" and start <= v[1] <= end}
    check(f"interval {start}-{end} set", set(got) == want,
          f"got {len(got)} want {len(want)} diff {list(set(got)^want)[:3]}")
    for v in list(want)[:20]:
        g = got.get(v)
        check(f"interval {start}-{end} ac/an {v}", g and g["ac"] == ac[v] and g["an"] == an[v],
              f"got {g and (g['ac'],g['an'])} want {(ac[v],an[v])}")
lat.sort()
print(f"  interval query latency: avg {sum(lat)/len(lat):.1f}ms p50 {lat[len(lat)//2]:.1f}ms p95 {lat[int(len(lat)*.95)]:.1f}ms max {lat[-1]:.1f}ms")

# ---------- 2. 基因名查询 ----------
print("[2] gene queries")
for g in random.sample(sorted(genes), 5):
    chrom, s, e = genes[g]
    d, ms = api("/api/variants", {"gene": g, "limit": 5000})
    got = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in d["variants"]}
    want = {v for v in seen if v[0] == chrom and s <= v[1] <= e}
    check(f"gene {g}", got == want, f"got {len(got)} want {len(want)}")

# ---------- 3. 样本/人群条件组合查询 ----------
print("[3] sample/population filter queries")
for s in random.sample(sorted(sample_vars), 10):
    d, _ = api("/api/variants", {"sample": s, "limit": 5000})
    got = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in d["variants"]}
    check(f"sample {s}", got == sample_vars[s], f"got {len(got)} want {len(sample_vars[s])}")
for pop in ["EAS", "AFR"]:
    want = {v for v, ss in carriers.items() if any(sample_meta.get(s, ("", ""))[1] == pop for s in ss)}
    d, _ = api("/api/variants", {"population": pop, "limit": 5000})
    got = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in d["variants"]}
    check(f"pop {pop}", got == want, f"got {len(got)} want {len(want)}")

# ---------- 4. 频率统计对比(总体 + 分人群) ----------
print("[4] allele frequencies")
vs = random.sample(all_vars, 15)
for v in vs:
    d, _ = api("/api/frequencies", {"chrom": v[0], "start": v[1], "end": v[1]})
    rows = [x for x in d["variants"] if (x["chrom"], x["pos"], x["ref"], x["alt"]) == v]
    check(f"freq found {v}", len(rows) == 1)
    if rows:
        row = rows[0]
        check(f"freq total {v}", row["ac"] == ac[v] and row["an"] == an[v],
              f"got {row['ac']}/{row['an']} want {ac[v]}/{an[v]}")
        for pop, (pac, pan) in ((p, (pop_ac[v][p], pop_an[v][p])) for p in pop_ac[v]):
            got = row["populations"].get(pop)
            check(f"freq pop {pop} {v}", got and got["ac"] == pac and got["an"] == pan,
                  f"got {got} want {pac}/{pan}")

# ---------- 5. 重复样本识别 ----------
print("[5] duplicate detection")
d, _ = api("/api/duplicates")
declared = {x["person_id"] for x in d["declared"]}
want_declared = {p for p, _ in sample_meta.values()
                 if sum(1 for q, _ in sample_meta.values() if q == p) > 1}
check("declared dup persons", declared == want_declared,
      f"got {len(declared)} want {len(want_declared)}")
# 指纹:暴力侧同样计算
fp = defaultdict(set)
for s, vs in sample_vars.items():
    sig = tuple(sorted(vs))
    fp[sig].add(s)
want_fp_groups = {tuple(sorted(g)) for g in fp.values() if len(g) > 1}
got_fp_groups = {tuple(sorted(x["samples"])) for x in d["detected"]}
check("fingerprint dup groups", got_fp_groups == want_fp_groups,
      f"got {len(got_fp_groups)} want {len(want_fp_groups)} diff {list(got_fp_groups ^ want_fp_groups)[:3]}")
# 相似度检测:40 组重复送检应全部被找回
sim_pairs = {tuple(sorted((x["sample_a"], x["sample_b"]))) for x in d.get("similar", [])}
person_samples = defaultdict(list)
for s, (p, _) in sample_meta.items():
    person_samples[p].append(s)
dup_person_pairs = set()
for p, ss in person_samples.items():
    if len(ss) > 1:
        for i in range(len(ss)):
            for j in range(i + 1, len(ss)):
                dup_person_pairs.add(tuple(sorted((ss[i], ss[j]))))
missing = dup_person_pairs - sim_pairs
check("similarity recovers all re-submitted pairs", not missing,
      f"missing {len(missing)}/{len(dup_person_pairs)}: {list(missing)[:3]}")
# 抽查 API 报出的相似对,暴力侧重算重叠系数(同样只用携带者 2..200 的变异)
carrier_n = {v: len(ss) for v, ss in carriers.items()}
rare = {v for v, n in carrier_n.items() if 2 <= n <= 200}
for x in d.get("similar", [])[:20]:
    a, b = x["sample_a"], x["sample_b"]
    sa = sample_vars[a] & rare
    sb = sample_vars[b] & rare
    ov = len(sa & sb) / max(1, min(len(sa), len(sb)))
    check(f"similar pair {a}-{b} overlap", ov >= 0.4, f"brute overlap {ov:.3f}")

# ---------- 6. 导出 VCF 合法性 ----------
print("[6] VCF export")
center = random.choice(all_vars)[1]
params = {"chrom": "chr1", "start": center - 500, "end": center + 500, "population": "EAS"}
url = API + "/api/export.vcf?" + urllib.parse.urlencode(params)
text = urllib.request.urlopen(url).read().decode()
lines = text.strip().split("\n")
header = [l for l in lines if l.startswith("#")]
data = [l for l in lines if not l.startswith("#")]
ok = header and header[-1].startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
ok = ok and all(len(l.split("\t")) >= 10 for l in data)
ok = ok and all(l.split("\t")[3] and l.split("\t")[4] for l in data)
check("export vcf format", ok, f"{len(data)} rows")
# 导出内容与查询一致
d, _ = api("/api/variants", {**params, "limit": 5000})
want = {(v["chrom"], str(v["pos"]), v["ref"], v["alt"]) for v in d["variants"]}
got = {(l.split("\t")[0], l.split("\t")[1], l.split("\t")[3], l.split("\t")[4]) for l in data}
check("export vcf content", got == want, f"got {len(got)} want {len(want)}")

# ---------- 7. 服务端基准 ----------
d, _ = api("/api/bench", {"n": 300})
print(f"[7] server-side interval bench: n={d['n_queries']} avg {d['avg_ms']:.2f}ms "
      f"p50 {d['p50_ms']:.2f}ms p95 {d['p95_ms']:.2f}ms p99 {d['p99_ms']:.2f}ms max {d['max_ms']:.2f}ms")

print(f"\n{'ALL PASS' if fails == 0 else str(fails) + ' FAILURES'}")
sys.exit(1 if fails else 0)
