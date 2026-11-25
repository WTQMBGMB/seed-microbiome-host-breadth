# Environment / 运行环境

- R 4.5.2 (Windows/Linux); exact package versions in `required_packages.tsv`.
- Run R with a UTF-8 locale supported by your operating system. On Windows,
  an inherited `LC_ALL=C.UTF-8` can make R escape Unicode in table labels.
- pandoc 3.1.11.1 or newer (bundled with RStudio is also fine).
- `metagenomeSeq` is a Bioconductor package; install once with
  `BiocManager::install("metagenomeSeq")`. Other packages come from CRAN,
  except `ggtree`/`ggtreeExtra` (Bioconductor) and `itol.toolkit`.
- Python 3.12 or later is only needed for `scripts/tools/*.py` (standard library only).
- Fonts: figures use Arial; on Linux install a metric-compatible Arial
  (`ttf-mscorefonts-installer` or Liberation Sans as fallback).
- MAFFT, trimAl and FastTreeMP are only needed for optional server rebuilds;
  the default run reads `server_result/`.

R 4.5.2；包版本见 `required_packages.tsv`。`metagenomeSeq`、`ggtree` 等来自
Bioconductor。默认运行不需要 MAFFT/trimAl/FastTreeMP，仅服务器重建时需要。
Python 仅用于检查脚本，依赖标准库。
