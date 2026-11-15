# variantdb — 变异位点库

面向临床基因检测场景的变异位点库:VCF 批量入库、位点标准化合并、区间/基因/样本组合查询、
人群频率统计、质控过滤留痕、重复送检识别、子集导出合法 VCF,自带检索页面。

## 架构

```
variantdb/            Rust 后端 (axum + rusqlite/SQLite, 单文件数据库)
  src/vcf.rs          VCF 解析 + 多等位拆分 + 左对齐标准化
  src/db.rs           schema、事务化批量导入(失败整体回滚)
  src/main.rs         REST API + 基准测试端点
  static/index.html   检索页面(样本管理/批次导入/位点检索/频率统计/质控报告)
tools/
  gen_data.py         合成参考基因组 + 几千个样本 VCF(混入未标准化写法、重复样本)
  verify.py           独立实现的规范化 + 暴力全扫描,与 API 结果逐项对比
data/                 生成的数据与 variant.db
```

## 构建与运行

```bash
cd variantdb && cargo build --release
python3 -I tools/gen_data.py data          # 生成合成数据
./variantdb/target/release/variantdb \
    --db data/variant.db --ref data/ref.fa --genes data/genes.tsv --port 8000
# 打开 http://localhost:8000
```

导入批次(API 或页面"批次导入"):

```bash
curl -X POST localhost:8000/api/import -H 'Content-Type: application/json' -d '{
  "name": "batch1", "vcf_dir": "/abs/path/data/batch1",
  "sample_sheet": "/abs/path/data/samples_batch1.csv",
  "min_qual": 30, "min_dp": 10, "min_gq": 20 }'
```

## 关键设计

1. **标准化**:每个 ALT 等位单独拆分;不动点循环做"去共同后缀 → 去共同前缀 → 空等位借参考
   碱基补锚左移",等价于 bcftools norm 的最小化+左对齐。不同批次的不同写法归一到同一
   `(chrom,pos,ref,alt)`,靠该列上的 UNIQUE 约束合并成一条位点记录。
2. **事务导入**:一个批次一个事务,任何一条记录解析/插入失败即 `ROLLBACK`,批次行本身也在
   事务内,库里不留半截数据。
3. **质控**:QUAL(位点级)、DP/GQ(基因型级)阈值在导入时指定;不达标记录写入
   `qc_failures`(含原因和实测值),不进入频率统计。
4. **频率**:AC/AN 按通过质控的基因型实时汇总,可按人群分组。
5. **重复送检**:三条路径——样本表 `person_id` 声明、携带基因型指纹完全一致、
   稀有变异(携带者 2..200)携带集合重叠系数 ≥0.5 的相似对。
6. **导出**:按检索条件过滤后重建 VCFv4.2(头部 + GT:DP:GQ),位点已是标准化写法。

## 主要 API

| 端点 | 说明 |
|---|---|
| `POST /api/import` | 导入批次(事务,失败回滚),返回耗时与吞吐 |
| `GET /api/variants` | 区间/基因/样本/个体/人群组合查询 |
| `GET /api/frequencies` | 位点频率(总体 + 分人群) |
| `GET /api/samples` `/api/duplicates` | 样本管理、重复送检识别 |
| `GET /api/qc` | 质控失败原因汇总与明细 |
| `GET /api/export.vcf` | 按条件导出合法 VCF |
| `GET /api/bench` | 服务端随机区间查询基准 |

## 验证

```bash
python3 -I tools/verify.py data http://localhost:8000
```

暴力扫描全部 VCF(独立实现的规范化),与 API 对比:随机区间查询的位点集合与 AC/AN、
基因查询、样本/人群条件查询、总体及分人群频率、重复样本识别、导出 VCF 内容,
并报告区间查询延迟。
