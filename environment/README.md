# Environment / 运行环境

Validated on R 4.5.2 (2025-10-31), Windows 11. The analysis reads local data and does not install packages automatically.

在 R 4.5.2（2025-10-31）、Windows 11 上验证。流程只读取本地数据，不在分析过程中自动安装包。

## R

- Required R packages are listed in `required_packages.tsv`; `renv.lock` records the current workstation library.
- `renv.lock` is an environment record, not a clean-room-validated lock. Mixed Bioconductor releases and an lme4/Rcpp declared-version mismatch were present in the source project. Do not describe `renv::restore()` as guaranteed.
- First-time restore (manual, once):

```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore(project = ".", lockfile = "environment/renv.lock", prompt = FALSE)
```

- `Analysis_pipeline.Rmd` starts with a preflight chunk that checks required packages and stops with a bilingual message when something is missing.

## Python

Map generation uses Python 3 with `numpy`, `pandas`, `geopandas`, `matplotlib`, and `shapely`.
Set `SEED_PYTHON` to the interpreter if `python` is not on `PATH`.

```bash
python -m pip install numpy pandas geopandas matplotlib shapely
```

## External tools

Only needed if optional server-style recomputation is enabled (`ALLOW_SERVER_RECOMPUTE=true`):
MAFFT, trimAl, FastTreeMP, IQ-TREE 3, MUSCLE, and GNU parallel. The default run uses `server_result/` and does not require them.

## Threads

Every local or optional step is capped at 16 threads (`SEED_THREADS`, default 16). External commands must be given the same cap explicitly.

## Pandoc

HTML output is self-contained. RStudio provides Pandoc; outside RStudio set `RSTUDIO_PANDOC` or put `pandoc` on `PATH`.
