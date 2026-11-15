#!/usr/bin/env python3
"""生成合成参考基因组 + 几千个样本的 VCF(故意混入未标准化写法和重复样本)"""
import random, os, sys, csv

random.seed(42)
OUT = sys.argv[1] if len(sys.argv) > 1 else "data"
os.makedirs(OUT, exist_ok=True)

CHROM = "chr1"
REFLEN = 100_000

# ---------- 1. 合成参考基因组:随机序列 + 串联重复/均聚物区域 ----------
bases = list("ACGT")
ref = [random.choice(bases) for _ in range(REFLEN)]
repeat_regions = []  # (start0, end0) 0-based,用于放置 indel
pos = 5000
while pos < REFLEN - 2000:
    kind = random.choice(["homo", "tandem2", "tandem3"])
    if kind == "homo":
        unit, n = random.choice(bases), random.randint(8, 20)
    elif kind == "tandem2":
        unit, n = random.choice(["AT", "GC", "CA", "TG"]), random.randint(5, 12)
    else:
        unit, n = random.choice(["GAT", "CAG", "TTA"]), random.randint(4, 10)
    seq = unit * n
    ref[pos:pos+len(seq)] = list(seq)
    repeat_regions.append((pos, pos + len(seq)))
    pos += random.randint(3000, 6000)

with open(f"{OUT}/ref.fa", "w") as f:
    f.write(f">{CHROM}\n")
    for i in range(0, REFLEN, 60):
        f.write("".join(ref[i:i+60]) + "\n")

# ---------- 2. 合成基因 ----------
genes = []
gpos = 2000
for i in range(1, 13):
    start = gpos
    end = start + random.randint(4000, 8000)
    genes.append((f"GENE{i}", CHROM, start, end))
    gpos = end + random.randint(1000, 3000)
with open(f"{OUT}/genes.tsv", "w") as f:
    for g in genes:
        f.write("\t".join(map(str, g)) + "\n")

# ---------- 3. 变异池:每个变异存标准(左对齐)形式 ----------
def normalize(seq, pos, r, a):
    """与库内一致的左对齐标准化,用于把池内变异规范为唯一形式"""
    r, a = list(r), list(a)
    while True:
        if r == a:
            break
        changed = False
        while r and a and r[-1] == a[-1]:
            r.pop(); a.pop(); changed = True
        while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
            r.pop(0); a.pop(0); pos += 1; changed = True
        if not r or not a:
            if pos > 1:
                b = seq[pos - 2]; r.insert(0, b); a.insert(0, b); pos -= 1; changed = True
            else:
                idx = pos - 1 + max(len(r), len(a))
                b = seq[idx] if idx < len(seq) else "N"
                r.append(b); a.append(b); changed = True
        if not changed:
            break
    return pos, "".join(r), "".join(a)

refseq = "".join(ref)

def norm_snp():
    p = random.randint(100, REFLEN - 100)
    r = ref[p-1]
    a = random.choice([b for b in bases if b != r])
    return (p, r, a)

def norm_indel():
    s, e = random.choice(repeat_regions)
    p = random.randint(s, e - 2)
    if random.random() < 0.5:  # 缺失:删 1~4 bp
        ln = random.randint(1, 4)
        r = "".join(ref[p-1:p+ln])  # 锚定 1bp + 删除序列
        return (p, r, r[0])
    else:  # 插入 1~4 bp
        ins = "".join(random.choice(bases) for _ in range(random.randint(1, 4)))
        a = ref[p-1]
        return (p, a, a + ins)

pool = []
seen = set()
while len(pool) < 2000:
    v = norm_snp() if random.random() < 0.7 else norm_indel()
    # 规范为左对齐唯一形式后去重(重复区里不同原始写法可能是同一变异)
    cp, cr, ca = normalize(refseq, v[0], v[1], v[2])
    if cr == ca or len(cr) == 0 or len(ca) == 0:
        continue
    if (cp, cr, ca) not in seen:
        seen.add((cp, cr, ca))
        pool.append((cp, cr, ca))
pool.sort()

# 一部分位点人为设成多等位(同一位置两个 ALT)
multiallelic = set()
for i, (p, r, a) in enumerate(pool):
    if len(r) == 1 and len(a) == 1 and random.random() < 0.08:
        multiallelic.add(p)

# ---------- 4. 非标准化写法生成器(全部保持等价变换) ----------
def denormalize(p, r, a):
    """把标准写法改写成等价但非标准的写法"""
    choice = random.random()
    if choice < 0.4 and len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        # 右移:去掉共同首碱基,双等位补同一个参考下游碱基(左对齐的逆操作)
        k = random.randint(1, 4)
        for _ in range(k):
            if len(r) > 1 and len(a) > 1 and r[0] == a[0] and p - 1 + len(r) < REFLEN:
                b = ref[p - 1 + len(r)]
                r, a = r[1:] + b, a[1:] + b
                p += 1
        return p, r, a
    elif choice < 0.75:
        # 右侧补 padding(共同后缀)
        pad = "".join(ref[p-1+len(r):p-1+len(r)+random.randint(1, 3)])
        if not pad:
            return p, r, a
        return p, r + pad, a + pad
    else:
        # 左侧补 padding(共同前缀)
        k = random.randint(1, 3)
        if p - 1 - k < 0:
            return p, r, a
        pad = "".join(ref[p-1-k:p-1])
        return p - k, pad + r, pad + a

# ---------- 5. 样本与批次 ----------
POPS = ["EAS", "EUR", "AFR", "SAS", "AMR"]
N_BATCHES = 6
SAMPLES_PER_BATCH = 500

persons = []
samples = []  # (sample_name, person_id, pop, batch)
n = 0
for b in range(1, N_BATCHES + 1):
    for i in range(SAMPLES_PER_BATCH):
        n += 1
        person = f"P{n:05d}"
        pop = random.choice(POPS)
        persons.append((person, pop))
        samples.append((f"S{n:05d}", person, pop, b))
# 重复送检:挑 40 个人,在后面的批次里再送一次(不同样本名)
dup_people = random.sample(persons, 40)
extra = []
for j, (person, pop) in enumerate(dup_people):
    b = random.randint(1, N_BATCHES)
    extra.append((f"S{90000+j}", person, pop, b))
samples += extra

# 每个人的基因型:对池子里每个位点按群体频率抽样
# 位点频率:不同人群略有差异
var_af = {}
for (p, r, a) in pool:
    base = random.choice([0.001, 0.005, 0.02, 0.05, 0.15, 0.3])
    var_af[(p, r, a)] = {pop: min(0.99, max(0.0005, base * random.uniform(0.5, 1.5))) for pop in POPS}

# 每个人一套基因型(重复样本共用同一人基因型 -> 指纹一致)
person_gt = {}
for person, pop in persons:
    gts = {}
    for v in pool:
        af = var_af[v][pop]
        gts[v] = 2 if random.random() < af*af else (1 if random.random() < 2*af*(1-af) else 0)
    person_gt[person] = gts

# ---------- 6. 写 VCF ----------
MIN_QUAL, MIN_DP, MIN_GQ = 30, 10, 20
by_batch = {}
for s in samples:
    by_batch.setdefault(s[3], []).append(s)

for b, slist in sorted(by_batch.items()):
    d = f"{OUT}/batch{b}"
    os.makedirs(d, exist_ok=True)
    with open(f"{OUT}/samples_batch{b}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "person_id", "population"])
        for s in slist:
            w.writerow([s[0], s[1], s[2]])
    for (sname, person, pop, _) in slist:
        gts = person_gt[person]
        lines = []
        # 每个样本随机携带池子的子集进行记录(模拟检出集合)
        chosen = [v for v in pool if gts[v] > 0 or random.random() < 0.05]
        # 按位置合并多等位
        by_pos = {}
        for v in chosen:
            by_pos.setdefault(v[0], []).append(v)
        for p in sorted(by_pos):
            vs = by_pos[p]
            # 基因型与质控字段(大部分达标,少量触发过滤)
            qual = round(random.uniform(28, 99), 1)
            dp = random.randint(6, 60)
            gq = random.randint(15, 60)
            if len(vs) > 1 and p in multiallelic:
                # 多等位合并写法:一行两个 ALT
                (p1, r1, a1), (p2, r2, a2) = vs[0], vs[1]
                gt1, gt2 = gts[vs[0]], gts[vs[1]]
                gt = f"{1 if gt1 else 0}/{2 if gt2 else 0}" if (gt1 or gt2) else "0/0"
                lines.append(f"{CHROM}\t{p1}\t.\t{r1}\t{a1},{a2}\t{qual}\tPASS\t.\tGT:DP:GQ\t{gt}:{dp}:{gq}")
            else:
                for v in vs:
                    (vp, vr, va) = v
                    ac = gts[v]
                    gt = {0: "0/0", 1: "0/1", 2: "1/1"}[ac]
                    # 30% 的位点故意写成非标准形式
                    if random.random() < 0.3:
                        vp, vr, va = denormalize(vp, vr, va)
                    lines.append(f"{CHROM}\t{vp}\t.\t{vr}\t{va}\t{qual}\tPASS\t.\tGT:DP:GQ\t{gt}:{dp}:{gq}")
        with open(f"{d}/{sname}.vcf", "w") as f:
            f.write("##fileformat=VCFv4.2\n##source=synthetic\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sname + "\n")
            f.write("\n".join(lines) + "\n")

print(f"done: {len(samples)} samples ({len(extra)} 重复送检), {len(pool)} pooled variants, {N_BATCHES} batches -> {OUT}")
