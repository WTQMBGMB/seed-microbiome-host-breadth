# Server-equivalent computations / 服务器等价计算

The default pipeline never runs these steps. It reads the archived outputs in
`server_result/`; missing files stop the run with a message pointing here.
Only set `SEED_ALLOW_SERVER_RECOMPUTE=true` if you intend to recompute locally.

默认流程不运行本目录脚本，只读取 `server_result/` 中的现成结果；文件缺失时
流程停止并提示本目录。只有显式设置 `SEED_ALLOW_SERVER_RECOMPUTE=true` 才本地重算。

All scripts use `<SERVER_BASE>` as a placeholder for the working directory; edit
it or export `BASE_DIR` before running. Paths to Rscript/MAFFT/trimAl/FastTree
are intentionally not hard-coded.

所有脚本以 `<SERVER_BASE>` 作为工作目录占位符；运行前修改或导出 `BASE_DIR`。
Rscript、MAFFT、trimAl、FastTree 的路径未硬编码。

## Contents / 内容

| Step | Script | Output consumed by |
|---|---|---|
| Ch2 mixed-pool adonis, 1999 permutations | `ch2_mixed_pool_adonis.R <work_dir> 1999 <ncore>` | `server_result/ch2/mixed_pool_adonis.rds` |
| Ch2 per-region Bray distances | `ch2_beta_distances.R` (optional; see note) | `server_result/ch2/beta/bray_dist_*.rds` |
| Ch2 Part A/B adonis objects | same Ch2 code as the Rmd (`params$force_recompute`) | `server_result/ch2/adonis/part{A_units,B_pools}.rds` |
| Ch3 ASV tree + mapped-bin backend | `asv_host_phylo/server_build_candidate_trees.sh`, then `00`–`20` scripts | `server_result/ch3/` |
| Ch4 GOS tree (MAFFT→trimAl→FastTreeMP→midpoint) | `ch4/run_ch4_server.sh` | `server_result/ch4/gos_hostn_pos_16s_midpoint.nwk` |
| Ch4 λ/K and Fritz–Purvis D | `ch4/run_ch4_server.sh` | `server_result/ch4/Tab_4_tree_*.tsv` |

## Notes / 说明

- Ch2 β distances and Part A/B adonis objects were produced by the same code as
  `Ch2_host_effect.Rmd` on the server. If the archived objects are absent, the
  chapter recomputes them locally when `params$force_recompute = TRUE`; this is
  the supported fallback. / Ch2 的 β 距离与 Part A/B adonis 对象由同一章代码在服务器
  生成；缺失时可在章内重算。
- Ch3 ASV backend scripts were written for a fixed server layout; they are kept
  for provenance and optional rebuilding. / Ch3 后端脚本保留用于溯源与可选重建。
- Ch4 tree construction requires MAFFT, trimAl and FastTreeMP; the R statistics
  require `phylosig_fast.R`, `phylo_d_fast.R` and the helper scripts in `ch4/`.
  / Ch4 建树需要 MAFFT、trimAl、FastTreeMP；统计需要 `ch4/` 中的 R 脚本。
