# Data / 输入数据

Download the inputs from Zenodo by running `python scripts/tools/download_data.py`
from the repository root. The archive restores `data/` and `server_result/`.
/ 在仓库根目录运行下载脚本，解压恢复输入数据和服务器存档。

## Layout / 目录

```text
data/
├─ filtered_merge/{16S,ITS}/        filtered ASV tables, rep-seqs, taxonomy
├─ filtered_region/{6 regions}/     filtered per-region tables
├─ standard_merge/                  CSS-normalised merged tables (prepared)
├─ standard_region/{6 regions}/     CSS-normalised region tables (prepared)
├─ metadata/                        metadata_16S.txt, metadata_ITS.txt
├─ function_prediction/             BugBase, FAPROTAX2, metaTraits, PICRUSt2 KO
├─ db/                              KO level and PGPG annotation tables
├─ plant_tree/                      TimeTree export and derived host trees
├─ gis/                             Natural Earth 1:110m boundaries
├─ manual/                          fixed phylum colour map
└─ precomputed/                     alpha-diversity caches, integrated trait table
```

## Raw vs prepared / 原始与预计算

- **Raw**: `filtered_*`, `metadata/`, `function_prediction/`, `db/`,
  `plant_tree/species_all.nwk`, `gis/`, `manual/`.
- **Prepared**: `standard_*` (CSS tables) and `precomputed/` (alpha caches and
  `taxonomy_16s_traits.txt`). They are shipped for speed and exact reproduction;
  `00_preprocessing.Rmd` and Chapter 5 can regenerate them when absent.
- β distances and Chapter 2 adonis objects are server results and live in
  `server_result/ch2/`, not here.

`MANIFEST.tsv` lists every file with size and SHA256. Files not listed are not
used by the pipeline.
