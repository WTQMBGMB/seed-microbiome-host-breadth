# Optional server-style recomputation / 可选服务器等价重算

默认流程不运行本目录。所有高耗时结果已放在 `server_result/`，`Analysis_pipeline.Rmd` 只读取并校验。

需要本地重算时：

1. 设置 `ALLOW_SERVER_RECOMPUTE=true`。
2. 设置 `SEED_RECOMPUTE` 为具体步骤（如 `adonis_v4`、`gos`），或 `all`。
3. 保证线程不超过 16（`SEED_THREADS=16`）。

当前本地可重算的步骤在对应章节内实现：

- β 距离 / adonis2：`analysis/02_alpha_beta.Rmd`
- GOS 零模型：`analysis/04_generalists_specialists.Rmd`
- bNTI：默认只读 `server_result/bnti/local` 的 pair 表；历史 iCAMP 引擎未归档，不提供等价重算。
- 全 ASV 树：默认复用 `server_result/gos/tree_16s_allASV`；建树命令见 `server_result/README.md` 和原项目归档。

本目录不包含自动安装或自动下载逻辑。
