# 投稿代码编写规范 / Code style

版本 v2.0，2026-11-12。适用于 `Downstream_analysis/github2`。本目录是从 `filtered_merge`、`filtered_region` 和现成功能预测结果出发的可复现下游流程；不重跑原始数据下载、质控和功能预测。

## 组织 / Organization

- `Analysis_pipeline.Rmd` 是唯一入口；`analysis/` 按 5 个主图章节 + 预处理/表格拆分。
- 一次性分析直接写在章节 Rmd 中；只有跨章节复用、性能关键或服务器可选步骤保留为 `R/` 或 `scripts/server/` 文件。
- 同一章节内逻辑紧密的代码块可以共享对象；跨章节通过 `results/checkpoints/` 和 `server_result/` 显式传递。
- 不维护“最终版”“副本”或多个同义入口。

## 数据 / Data

- `data/` 与 `server_result/` 只读；`results/` 是唯一可写区域。
- `standard_merge`/`standard_region` 是生成结果，写入 `results/derived/css/`；metadata 放在 `data/metadata/`。
- 输入不原位覆盖；派生结果按分析命名目录写入。
- 每个输入和 server 结果都在 `MANIFEST.tsv` 中记录来源、大小和 SHA256。

## 注释 / Comments

- 必要注释中英文并列，每个逻辑步骤一条简短说明。
- 不逐行解释赋值、读表、建目录或绘图函数。
- 科学问题、统计口径和解释边界写在 Rmd 正文。

## 函数 / Functions

- 先写顺序分析；不为一次读表、一次检验或一次保存文件新建函数。
- 保留确有复用或性能需要的函数；`R/` 只放共享工具和性能代码。
- 公共函数使用显式参数，不依赖调用环境中的隐式全局变量。

## 复现 / Reproducibility

- 全部使用项目相对路径；不出现个人盘符、服务器账户或 `setwd()`。
- R 包版本记录在 `environment/required_packages.tsv` 和 `renv.lock`；分析不自动安装或更新包。
- 随机分析明确设置已批准种子；置换/bootstrap 次数和样本范围不得在重构中改变。
- 线程统一上限 16：`SEED_THREADS` 默认 16，R/BLAS/外部命令都遵守。
- 默认读取 `server_result/`；缺失资产时停止并指向重算入口，不静默重算。
- 明确区分“读取审核结果重新出图”和“重新拟合模型”。

## 图件 / Figures

- 每行 `Material_Description2.csv` 的 `Subplot` 输出一个独立文件；单行多面板图（FigS1、FigS2、FigS9、FigS12）按整图输出。
- 文件名 = `Subplot_Short_English_Description.pdf/.jpg`，例如 `Fig1A_Global_sample_distribution_map.pdf`。
- 图内不嵌 A/B/C 面板字母；最终拼版在 Illustrator 完成。
- 统一用 `save_panel()` 输出 PDF + 300 dpi JPG；不另建多套 `save_fig`。
- 非稿件诊断图放 `results/diagnostics/`，用 `diagnostic_` 前缀。

## 表格 / Tables

- 稿件表统一命名 `Table S1.xlsx`–`Table S6.xlsx`；代码内部历史名 `Table_S7` 不再使用。
- Table S1 是人工整理输入；Table S2/S3 由数据集章节生成，Table S4 由宿主系统发育章节生成，Table S5 由 GOS 分类表生成，Table S6 由 16S 性状表生成。

## 修改与验收 / Change control

- 修改前备份；科学方法或结论变更须作者确认。
- 重构后用 `scripts/tools/check_outputs.py` 检查 83 个 Subplot 是否齐全，并对关键结果表做数值对齐。
- 不通过观察 JPG/TIF 判断结论；图件只检查底层数据和文件完整性。
