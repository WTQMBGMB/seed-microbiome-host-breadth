# Code and comment conventions / 代码与注释规范

Applies to the five chapter Rmd files, `00_preprocessing.Rmd` and `R/`.
适用于五个章节 Rmd、`00_preprocessing.Rmd` 与 `R/`。

## Comments / 注释

- Use bilingual comments where an explanation is needed:
  `# Explain the method or parameter / 说明方法或参数`.
- Keep comments only for: scientific method, thresholds/parameters,
  non-obvious handling, special cases, seeds, and input/output contracts.
- Do not keep: change history (`2026-…`, “previously”, “replaced”),
  prototype notes, debug prints, commented-out old code, decorative separators,
  or machine-specific paths.
- Do not change `#` inside strings, colours, regular expressions or code.

- 必要处使用中英双语注释。
- 只保留方法、阈值/参数、非直观处理、特殊情况、随机种子与输入输出约定。
- 删除修改历史、原型说明、调试输出、注释掉的旧代码、装饰分隔线与本机路径。
- 不改变字符串中的 `#`、颜色值、正则与代码行为。

## Paths and inputs / 路径与输入

- All paths resolve from `PROJECT_ROOT` (`R/config.R`); no drive letters,
  usernames or server paths in the pipeline.
- Inputs are grouped in `data/`; server results in `server_result/`.
- Never overwrite an input in place. Derived files go to `result/` (or
  `result/prepared/` when a shipped prepared table is regenerated).
- When a prepared file exists, read it and skip the computation; when it is
  missing, either regenerate (if cheap and deterministic) or stop with a clear
  message.

- 所有路径从 `PROJECT_ROOT` 解析，不出现盘符、用户名或服务器路径。
- 输入集中在 `data/`，服务器结果集中在 `server_result/`。
- 不原位覆盖输入；派生结果写入 `result/`。
- 预计算文件存在即读取并跳过计算；缺失时重算或明确报错。

## Analysis code / 分析代码

- Do not change seeds, permutation/bootstrap counts, thresholds, FDR scopes,
  model formulas or filtering rules when editing paths or comments.
- Shared helpers live in `R/`; chapter-specific logic stays in the chapter.
- Avoid duplicate helper definitions within a chapter; give distinct helpers
  distinct names.
- Server-equivalent steps are guarded by `SEED_ALLOW_SERVER_RECOMPUTE`; the
  default is read-only.

- 修改路径或注释时，不改变种子、置换/bootstrap 次数、阈值、FDR 范围、
  模型公式与筛选规则。
- 公共函数放 `R/`，章节特有逻辑留在章节内。
- 避免同章重复定义同名辅助函数。
- 服务器等价步骤由 `SEED_ALLOW_SERVER_RECOMPUTE` 控制，默认只读。

## Outputs / 输出

- Figures use `save_fig()` (PDF + JPG, millimetres); Chapter 3 keeps its
  original stems and disables the `Fig_<n>_` check.
- Each saved figure is embedded in place with `embed_fig()`; the rendered
  document uses the JPG with a path relative to the HTML output.
- Tables use `write_tab()` and keep complete numeric precision.
- Do not rename existing output stems; manuscript mapping depends on them.
- `manifests/output_manifest.tsv` defines the expected output set.

- 图件用 `save_fig()` 输出 PDF+JPG；第三章沿用原 stem 并关闭命名检查。
- 每张图随后用 `embed_fig()` 原位嵌入渲染文档，使用 JPG 与相对 HTML 的路径。
- 表格用 `write_tab()`，保留完整数值精度。
- 不重命名现有输出 stem；稿件映射依赖这些名称。
- 期望输出集合见 `manifests/output_manifest.tsv`。

## Verification / 验证

- After changing a chapter, render it and compare tables against the reference
  with `scripts/tools/compare_results.py` (numeric tolerance, not file hash).
- A full clean run from an empty `result/` is the acceptance test.
- Missing or mismatching outputs are errors, not warnings to ignore.

- 修改章节后渲染，并用 `compare_results.py` 按数值容差比对参考结果。
- 从空 `result/` 完整运行是验收标准。
- 输出缺失或数值不一致按错误处理。
