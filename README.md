# Seed microbiome meta-analysis: reproducible downstream pipeline
**种子微生物组 meta 分析：可复现下游流程**

This repository rebuilds the downstream analysis from filtered feature
tables, metadata, functional predictions and archived server results. It
generates the analysis figures and tables from the accompanying dataset.
The manuscript figures include separately assembled panels and iTOL exports.

本仓库从 filtered ASV 表、metadata、功能预测结果与服务器存档出发，重建 v2
下游分析。默认运行生成下游分析图表，使用归档的服务器运算结果。稿件主图包含单独拼版与 iTOL 导出。

Companion processed dataset (version 2.0.0): https://doi.org/10.5281/zenodo.24523342

## 1. Quick start / 快速开始

1. Install R 4.5.2 and packages listed in `environment/required_packages.tsv`
   (see `environment/README.md`).
2. Run `python scripts/tools/download_data.py` from the repository root.
   The script downloads and checks the input archive listed in
   `manifests/data_download.tsv` and extracts `data/` and `server_result/`.
3. Run the whole pipeline:

```bash
Rscript 00_run_all.R
```

or render chapters individually:

```bash
Rscript scripts/render_one.R 00_preprocessing.Rmd result/html
Rscript scripts/render_one.R Ch1_dataset.Rmd        result/html
# ... Ch2 to Ch5
```

4. Check outputs:

```bash
python scripts/tools/check_outputs.py
python scripts/tools/compare_results.py --ref <reference_result> --new result
```

1. 安装 R 4.5.2 与 `environment/required_packages.tsv` 中的包。
2. 在仓库根目录运行 `python scripts/tools/download_data.py`，下载并校验归档，
   解压得到 `data/` 和 `server_result/`。
3. `Rscript 00_run_all.R` 一次跑完；也可用 `scripts/render_one.R` 分章渲染。
4. 用两个 Python 脚本检查输出完整性与数值一致性。

## 2. Repository layout / 目录

```text
00_preprocessing.Rmd   CSS tables, host trees, alpha/beta, server assets
Ch1–Ch5_dataset...Rmd  Five analysis chapters
00_run_all.R           Preflight, render, output check
R/                     config, shared helpers, preprocessing, legacy plot helpers
scripts/server/        Optional server-equivalent computations
scripts/tools/         Data download, output check, result comparison
data/                  Inputs (not tracked; see data/README.md)
server_result/         Read-only server results from the archive + MANIFEST.tsv
manifests/             Output manifest and data download manifest
environment/           R/package versions and setup notes
result/                Generated figures, tables, trees and caches
```

## 3. Server results / 服务器结果

`server_result/` is read-only. The pipeline never writes there. If a required
file exists it is copied into `result/` and the step is skipped; if it is
missing the run stops and points to `scripts/server/`.

- Ch2: per-region Bray distances, Part A/B unit objects (adonis), mixed-pool
  adonis.
- Ch3: SILVA all-ASV tree, 500-bootstrap λ samples, mapped-bin backend tables
  and null model.
- Ch4: microbial tree, λ/K and Fritz–Purvis D statistics, build logs, iTOL
  annotation for Fig4A.
- Ch1: `host_overview_tree.pdf` (Fig1C), a manual/iTOL export with no generating
  code; it is provided as an asset and flagged in `server_result/MANIFEST.tsv`.

Local recomputation is opt-in only:

```bash
SEED_ALLOW_SERVER_RECOMPUTE=true Rscript 00_run_all.R
```

`server_result/` 只读，流程不会写入。文件存在时复制到 `result/` 并跳过计算；
缺失时停止并指向 `scripts/server/`。本地重算需显式设置环境变量。

## 4. Recompute switches / 重算开关

| Variable | Default | Meaning |
|---|---|---|
| `SEED_RECOMPUTE` | empty | `css`, `alphabeta`, `traits`, `all`; regenerate prepared inputs |
| `SEED_ALLOW_SERVER_RECOMPUTE` | false | allow local execution of server-equivalent steps |
| `SEED_NCORE` | 16 | thread cap for parallel steps |

Shipped CSS tables (`data/standard_*`) and alpha/trait caches
(`data/precomputed/`) are used by default; the switches force regeneration into
`result/prepared/`.

随仓库提供的 CSS 表与 α/性状缓存默认直接使用；上述开关可强制重算到
`result/prepared/`。

## 5. Outputs / 输出

- Figures: `result/ch*/figures/*.pdf` and `*.jpg` (76 files: 39 PDFs, 37 JPGs;
  two PDFs are exported without a JPG companion).
- Tables: `result/ch*/tables/Tab_*.tsv` (36 files).
- Trees and iTOL assets: `result/ch*/tree/`.
- Rendered reports: `result/html/*.html`, with each figure embedded inline.
- Caches and logs: `result/rds/`, `result/intermediate/`, `result/logs/` (not
  required for reading the results).

`manifests/output_manifest.tsv` lists every expected figure and table.
`scripts/tools/check_outputs.py` verifies presence; `compare_results.py`
compares table values against a reference run with a numeric tolerance.

图件 76 个文件、表格 36 个 TSV；清单见 `manifests/output_manifest.tsv`。
渲染报告位于 `result/html/*.html`，图件已内嵌在对应位置。

## 6. Analysis outline / 分析概要

1. **Dataset**: 6,349 samples (16S 4,012; ITS 2,337), 59 papers, 91 sources,
   45 countries, 122 host genera.
2. **Host effect**: within-study ICC/Δ and adonis, meta-analysis, and
   study-design variance partitioning.
3. **Host phylogeny**: Blomberg K / Pagel λ / Mantel / MRM / dbRDA, permutation
   nulls, and mapped-bin ASV phylogeny.
4. **Generalists and specialists**: representative-T host breadth, null-based
   threshold (4), O/E, HBS logistic models, microbial-tree signal.
5. **Traits**: BugBase/FAPROTAX2/metaTraits and PICRUSt2 KO associations with
   host breadth, merged and per-pool, with Seed/Germ consistency.

Methods, thresholds, seeds, permutation counts and FDR scopes are described in
the chapter Rmd files and unchanged from the analysed pipeline. See
`CODE_STYLE.md` for code and comment conventions.

## 7. Data provenance / 数据来源

- `data/MANIFEST.tsv`: size and SHA256 of every input, with raw/prepared role.
- `server_result/MANIFEST.tsv`: source, size, SHA256, generating step and
  parameters for every server result.
- `data/README.md`, `server_result/README.md`: inclusion and exclusion rules.
- Processed inputs and archived server results are distributed on Zenodo.
  Raw data for the authors’ sequencing are held by NGDC/GSA under
  PRJCA066385 (CRA044477). Source-study accessions are listed in the metadata.

## 8. Notes / 说明

- No absolute paths are used in code; all paths resolve from the repository
  root or `SEED_PROJECT_ROOT`.
- Figures are exported as editable PDF plus JPG and embedded inline in the
  rendered HTML through `embed_fig()`; panel letters are added during
  manuscript assembly, not here.
- Fig4A is an iTOL rendering: the repository supplies the annotation folder and
  statistics tables, not a rendered image.
- The trait table fallback keeps only columns used by Chapter 5; the shipped
  table additionally contains legacy GOS columns that Chapter 5 does not use.

- 代码不含绝对路径，均从仓库根目录或 `SEED_PROJECT_ROOT` 解析。
- 图件同时输出可编辑 PDF 与 JPG，并通过 `embed_fig()` 内嵌到渲染后的 HTML；
  分面字母在稿件排版阶段添加。
