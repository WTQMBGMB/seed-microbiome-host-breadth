# Archived server results / 服务器结果存档

Download the inputs from Zenodo by running `python scripts/tools/download_data.py`
from the repository root. The archive restores `data/` and `server_result/`.
/ 在仓库根目录运行下载脚本，解压恢复输入数据和服务器存档。

The default analysis copies required archived results into `result/`.
See `MANIFEST.tsv` for file size, SHA-256, and provenance.

## Contents / 内容

| Path | Content | Used by |
|---|---|---|
| `ch1/host_overview_tree.pdf` | Fig1C; manual/iTOL export, no generating code | Ch1 |
| `ch2/beta/bray_dist_{v4,v5_v7,its1}.rds` | per-region Bray-Curtis distances | Ch2, Ch3 |
| `ch2/adonis/{partA_units,partB_pools}.rds` | within-study and per-pool adonis/alpha objects | Ch2 |
| `ch2/mixed_pool_adonis.rds` | mixed-pool adonis2, 1999 permutations, 200 threads | Ch2 |
| `ch2/mixed_pool_raw/` | raw mixed-pool design and adonis outputs | provenance |
| `ch3/asv_ref_silva.nwk` | SILVA 138.2 all-ASV tree (22,845 tips) | Ch3 |
| `ch3/boot_lambda_samples.rds` | 500 species-bootstrap λ samples | Ch3 |
| `ch3/asv_host_phylo/` | mapped-bin Mantel/MRM/sensitivity and 999-permutation null | Ch3 |
| `ch4/gos_hostn_pos_16s_midpoint.nwk` | microbial tree (host breadth ≥ 1) | Ch4 |
| `ch4/Tab_4_tree_{phylosig,generalist_D}.tsv` | λ/K and Fritz–Purvis D statistics | Ch4 |
| `ch4/server_raw/` | build logs and raw phylosig/D tables | provenance |
| `ch4/itol_annotation/` | iTOL annotation for Fig4A | manuscript assembly |

## Policy / 规则

- Do not write here from the pipeline.
- Do not edit archived files; regenerate only via `scripts/server/` and update
  `MANIFEST.tsv`.
- `host_overview_tree.pdf` is the only result image that cannot be regenerated
  from raw data; it is kept as a manual asset and marked in the manifest.
- Fig4A is rendered in iTOL, not by the pipeline; only its annotation and
  statistics are archived.

- 流程不写入本目录；不直接编辑存档文件。
- `host_overview_tree.pdf` 是唯一无法由原始数据重建的结果图片，作为手工资产保留。
- Fig4A 在 iTOL 中渲染，此处仅存档注释与统计表。
