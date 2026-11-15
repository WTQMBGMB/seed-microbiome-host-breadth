# Output map / 输出映射

`manifests/figure_map.tsv` 是 `Paper_writing2/Material_Description2.csv` 的冻结版映射，共 83 行：

- `Subplot`：CSV 的 Subplot token（如 `Fig1A`、`FigS6C`、`TableS5A`）。
- `Order`：逻辑顺序 1–22。
- `Chapter`：对应的 `analysis/*.Rmd`。
- `Output_kind`：`figure` 或 `table`。
- `Expected_file`：`results/figures/` 或 `results/tables/` 下的实际文件名。

## 命名规则

```text
Fig{主图号}{子图字母}_{Short_English_Description}.pdf
Fig{主图号}{子图字母}_{Short_English_Description}.jpg
```

- 单行多面板的补充图（FigS1、FigS2、FigS9、FigS12）按 CSV 一行输出整图。
- 组合源图按 CSV Subplot 拆分为独立文件。
- 图内不嵌 A/B/C 字母；最终拼版在 Illustrator 完成。
- 表格文件固定为 `Table S1.xlsx`–`Table S6.xlsx`。

## 检查

```bash
python scripts/tools/check_outputs.py
```

检查 83 个 Subplot 是否全部有文件、PDF/JPG 是否成对、表格是否齐全，并报告未映射的多余图件。
