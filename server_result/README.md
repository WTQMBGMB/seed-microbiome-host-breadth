# Precomputed results / 现成重计算结果

本目录只读。默认 Knit 直接读取这些文件；缺失时停止并提示重算入口，不会静默重算。`MANIFEST.tsv` 记录原路径、大小和 SHA256。

## Contents

| 目录 | 内容 | 是否必须 |
|---|---|---|
| `beta/` | 六个 Bray-Curtis 距离 RDS、三个 adonis2 999 置换 checkpoint、beta source 分解表 | 是 |
| `phylo/` | λ bootstrap 样本、999 次系统发育置换 null、汇总表 | 是 |
| `bnti/` | 自测序 V5–V7 与 V4 Genus 的 pair 表、prep/pool checkpoint、输入树与 metadata | 是 |
| `gos/` | presence 矩阵、n_hosts/MPD 零模型、阈值、全 ASV 树与 iTOL 注释、GOS 下游结果表与桑基配色表 | 是 |
| `asv_beta_phylo/` | mapped-bin ASV 系统发育候选树、profiles、999 置换 null、paper_main 结果 | 是 |
| `sensitivity/` | Table S4 四个源表 | 是 |
| `function_prediction/` | 16S 性状表、trait/KO two-part 结果、KO/KEGG/PGPG 桑基节点与连线、source bootstrap | 是 |

## Recompute policy

- 默认 `ALLOW_SERVER_RECOMPUTE=false`。
- 需要本地重算时，设置 `ALLOW_SERVER_RECOMPUTE=true` 并运行 `scripts/server/` 中对应入口；所有线程硬上限 16。
- bNTI 历史计算使用修改版 iCAMP 引擎；当前归档 pair 表是稿件结果，本地重算只作为可选敏感性，不声称逐位等价。
- GOS 全 ASV 树由服务器 MAFFT → trimAl → FastTreeMP 构建；默认复用现成树，不重建。

The manuscript uses pooled source-bootstrap stability and full-data S/P comparisons. Historical S/P bootstrap rows are retained for provenance, not as additional manuscript evidence.
