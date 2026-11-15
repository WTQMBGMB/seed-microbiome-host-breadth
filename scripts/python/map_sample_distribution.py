#!/usr/bin/env python3
"""
map_sample_distribution.py — Global amplicon sample distribution map

参考 `Map heatmap_blue.py`，将其参数化：
  - 样本计数不再硬编码，从输入 TSV 读取（列：country, n_16s, n_its）
  - 国家名做同义词归一（metadata 中的 Country 写法 -> Natural Earth NAME）
  - 台湾几何合并入中国
  - 对数（log10）绿色渐变配色 + 共享色条
  - 输出 16S+ITS 合计地图，作为 Fig. 1A

用法：
  python map_sample_distribution.py \
      --input country_counts.tsv \
      --shapefile ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp \
      --outdir ../../result/dataset_overview/result_figure \
      --prefix "Fig. 1A_Sample distribution map"

输出：<outdir>/<prefix>.jpg
依赖：geopandas, pandas, numpy, matplotlib
"""
import argparse
import os
import sys

# Windows 控制台默认 GBK，无法打印部分国家名（如 Côte d'Ivoire 的 ô），强制 UTF-8
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# 国家名同义词（metadata Country -> Natural Earth 110m NAME）
# ---------------------------------------------------------------------------
SYNONYMS = {
    "USA": "United States of America",
    "United States": "United States of America",
    "U.S.A.": "United States of America",
    "U.S.": "United States of America",
    "US": "United States of America",
    "UK": "United Kingdom",
    "United Kingdom of Great Britain and Northern Ireland": "United Kingdom",
    "England": "United Kingdom",
    "Korea": "South Korea",
    "Republic of Korea": "South Korea",
    "Russia": "Russia",
    "Democratic Republic of the Congo": "Dem. Rep. Congo",
    "DRC": "Dem. Rep. Congo",
    "Congo (DRC)": "Dem. Rep. Congo",
    "Republic of Congo": "Congo",
    "Congo-Brazzaville": "Congo",
    "Ivory Coast": "Côte d'Ivoire",
    "Czech Republic": "Czechia",
    "Czechia": "Czechia",
    "United Republic of Tanzania": "Tanzania",
    "Tanzania, United Republic of": "Tanzania",
    "Central African Republic": "Central African Rep.",
    "Dominican Republic": "Dominican Rep.",
    "Equatorial Guinea": "Eq. Guinea",
    "Western Sahara": "W. Sahara",
    "South Sudan": "S. Sudan",
    "São Tomé and Príncipe": "São Tomé and Principe",
    "Sao Tome and Principe": "São Tomé and Principe",
    "Bosnia and Herzegovina": "Bosnia and Herz.",
    "North Macedonia": "Macedonia",
    "Macedonia (FYROM)": "Macedonia",
    "Myanmar": "Myanmar",
    "Turkey": "Turkey",
    "Vietnam": "Viet Nam",
    "Viet Nam": "Viet Nam",
    "Laos": "Laos",
    "Syria": "Syria",
    "Federated States of Micronesia": "Micronesia",
    "Palestinian Territory": "Palestine",
    "Palestine": "Palestine",
    "Brunei": "Brunei",
    "Bolivia": "Bolivia",
    "Iran": "Iran",
    "Netherlands": "Netherlands",
    "Timor-Leste": "Timor-Leste",
    # 以下在 R 数据准备层已做国名归一化，这里保留作兜底
    "Cote d Ivoire": "Côte d'Ivoire",
    "Italia": "Italy",
}


def main():
    ap = argparse.ArgumentParser(description="Fig. 1A 样本全球分布地图")
    ap.add_argument("--input", required=True,
                    help="输入 TSV：列 country / n_16s / n_its")
    ap.add_argument("--shapefile", required=True,
                    help="Natural Earth 110m 国家 shapefile 路径")
    ap.add_argument("--outdir", required=True, help="输出目录")
    ap.add_argument("--prefix", default="Fig. 1A_Sample distribution map",
                    help="JPG 输出文件名前缀")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- 读取样本计数 ----
    df = pd.read_csv(args.input, sep="\t")
    # 列名大小写兼容（Country / country / n_16s / n_its / total 均可）
    df.columns = [c.strip().lower() for c in df.columns]
    for col in ("country", "n_16s", "n_its"):
        if col not in df.columns:
            raise ValueError(f"输入缺少列: {col}")
    df["country"] = df["country"].astype(str).str.strip()
    # 应用同义词
    df["country_ne"] = df["country"].map(lambda c: SYNONYMS.get(c, c))
    df = df.rename(columns={"country_ne": "Country"})

    # ---- 读取世界地图 ----
    world = gpd.read_file(args.shapefile)

    # 台湾合并入中国
    if "Taiwan" in world["NAME"].values:
        china_idx = world[world["NAME"] == "China"].index
        taiwan = world[world["NAME"] == "Taiwan"]
        if not taiwan.empty and len(china_idx) > 0:
            world.loc[china_idx, "geometry"] = world.loc[
                china_idx, "geometry"].union(taiwan.iloc[0]["geometry"])
            world = world[world["NAME"] != "Taiwan"]

    map_names = set(world["NAME"].dropna().str.strip())
    data_names = set(df["Country"])
    unmatched = data_names - map_names
    matched = data_names & map_names
    print("=" * 60)
    print(f"匹配成功 ({len(matched)}):", sorted(matched))
    print(f"未匹配 ({len(unmatched)}):", sorted(unmatched))
    print("=" * 60)

    world = world.merge(df, how="left", left_on="NAME", right_on="Country")

    # ---- 配色 ----
    unselected_color = "#eeeeee"
    green_cmap = LinearSegmentedColormap.from_list(
        "green_gradient", ["#A2CCAD", "#2A4019"])

    world["total"] = (world["n_16s"].fillna(0).astype(int)
                      + world["n_its"].fillna(0).astype(int))
    panels = [("total", "16S + ITS")]

    # 对数变换（各面板共用色标，以最大计数定标）
    all_counts = pd.concat([world[c] for c, _ in panels]).fillna(0).astype(int)
    max_count = int(all_counts.max())
    tick_values = [1, 10, 100, 1000, 10000]
    tick_values = [t for t in tick_values if t <= max_count] or [1]
    tick_logs = np.log10(tick_values)
    vmax_log = np.log10(max_count)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6.5))
    axes = [ax]

    for ax, (col, title) in zip(axes, panels):
        vals = world[col].fillna(0).astype(int)

        unsel = world[vals == 0]
        sel = world[vals > 0]
        if not unsel.empty:
            unsel.plot(color=unselected_color, linewidth=0.6, ax=ax,
                       edgecolor="black")
        if not sel.empty:
            logv = np.log10(sel[col])
            norm = plt.Normalize(vmin=0.0, vmax=vmax_log)
            sel.plot(color=green_cmap(norm(logv)), linewidth=0.6, ax=ax,
                     edgecolor="black")
        ax.set_title(f"{title}  (n = {int(vals.sum())})", fontsize=14)
        ax.set_axis_off()

    # 图例 + 共享色条
    legend_elements = [
        Patch(facecolor=unselected_color, edgecolor="black", label="No samples"),
        Patch(facecolor=green_cmap(0.0), edgecolor="black", label="1"),
        Patch(facecolor=green_cmap(0.5), edgecolor="black", label="10"),
        Patch(facecolor=green_cmap(1.0), edgecolor="black",
              label=str(tick_values[-1])),
    ]
    axes[0].legend(handles=legend_elements, loc="lower left", fontsize=9)

    sm = ScalarMappable(norm=plt.Normalize(0, vmax_log), cmap=green_cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.6, pad=0.01)
    cbar.set_ticks(tick_logs)
    cbar.set_ticklabels([str(t) for t in tick_values])
    cbar.set_label("Number of samples (log10 scale)", fontsize=12)

    fig.suptitle("Global distribution of seed microbiome amplicon samples (16S + ITS combined)",
                 fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    jpg_path = os.path.join(args.outdir, f"{args.prefix}.jpg")
    pdf_path = os.path.join(args.outdir, f"{args.prefix}.pdf")
    fig.savefig(jpg_path, format="jpeg", dpi=args.dpi, bbox_inches="tight")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    print(f"已输出: {jpg_path}; {pdf_path}")


if __name__ == "__main__":
    main()
