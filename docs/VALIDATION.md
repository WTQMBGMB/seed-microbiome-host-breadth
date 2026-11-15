# Public release preparation checks on 2026-11-15

- All three data archives extracted successfully using the shipped downloader.
- SHA-256 and sizes matched for all 537 listed files.
- Output check passed: 84 mappings and 73 figure PDFs, with corresponding JPGs.
- Metadata: 3,652 16S and 2,310 ITS profiles.
- Host-breadth classification counts matched the manuscript.
- Robust17 / discordant9 sets and Table S1–S7 existence passed.
- All Python scripts, R helper files and R chapter chunks parsed successfully.
- CITATION.cff parsed and includes seven manuscript authors.

These are packaging and archived-result checks. The complete pipeline was not rerun
on a clean machine during release preparation. Historical checks follow below.

# Release preparation update

The historical validation below used the original 83-row map and Tables S1–S6. The public release adds Table S7 (84 rows), refreshes the 52-publication manual list, and provides portable validation. The numerical analysis inputs/results are unchanged. This is a record of the prior workstation check, not a clean-machine certification.

# Validation / 验证记录

验证只比较结果表、计数和文件完整性，不观察或解读 JPG/TIF 图片。

## 已通过

| 检查 | 结果 |
|---|---|
| metadata 行数 | 16S 3,652 + ITS 2,310 = 5,962 |
| adonis2 数值 | 18 行，与归档结果最大绝对差 0 |
| phylosig K / λ 估计 | 与归档结果最大绝对差 0 |
| Mantel / MRM / dbRDA 估计 | 最大绝对差 9.99e-16 |
| GOS 16S 分类 | NRH 12,190；generalist 1,095；intermediate 3,813；specialist 5,747 |
| GOS ITS 分类 | NRH 2,860；generalist 234；intermediate 522；specialist 1,302 |
| robust17 / discordant9 | 17 / 9 |
| Subplot 映射 | 83 行 |
| 输出文件 | 73 个图件 PDF/JPG 成对，Table S1–S6 齐全 |
| 桑基图几何 | 直接使用归档 nodes/links 的 `x`、`ymin`、`ymax`；16S/ITS/KO 柱宽分别为 0.4048、0.529、0.36 |
| 桑基图配色 | 16S 使用 `phylum_color_map_16s.tsv`；ITS 使用归档 PAL1 的全局频次升序；KO 使用归档 `colorlist2` 前 10/15 色 |
| Genus 3×3 | 16S 182 genera / 39 labels，11.25×9 in；ITS 45 genera / 8 labels，4×3 in；任一 part 显著即标注 |
| Boolean trait 3×3 | 119 traits / 45 labels，10.5×8.2 in；x(-0.75, 1)，y(-1, 1) |
| KO 3×3 | 9,110 个 both_finite KO，5×4 in；九宫格 `n=` 标注；0.5–99.5% 分位裁剪，169 个点越界 |

## 蒙特卡洛差异

Mantel、dbRDA、Blomberg's K 的置换 P 值与 ARCHIVE 结果可能在小范围内不同；点估计、R²、λ/K 和分类计数一致。原因是旧结果未固定统一的置换种子，本流程显式设置 `set.seed(123L)`。这是蒙特卡洛波动，不是方法变更。

## 复现命令

```bash
Rscript -e 'rmarkdown::render("Analysis_pipeline.Rmd")'
python scripts/tools/check_outputs.py
python scripts/tools/validate_results.py
```

验证报告写入 `results/diagnostics/validation_report.tsv`。

## 未覆盖

- 未在全新机器上验证 `renv::restore()`；见 `environment/README.md`。
- 未重跑原始数据下载、质控、分类和功能预测。
- 未重建 bNTI 历史修改版 iCAMP 引擎；归档 pair 表为稿件结果。
- 未以像素/视觉方式比较图件；只核对了底层数据和文件完整性。
