# -*- coding: utf-8 -*-
"""Compare key github2 outputs with the reviewed archive.

Images are never interpreted. Only tables and counts are compared.
"""
from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OLD = Path(r"D:\seed_meta_new\Downstream_analysis")
REPORT = ROOT / "results" / "diagnostics" / "validation_report.tsv"
rows: list[dict] = []


def add(check: str, status: str, detail: str = "") -> None:
    rows.append({"check": check, "status": status, "detail": detail})


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return math.nan


# 1. Dataset size
md16 = read_tsv(ROOT / "data" / "metadata" / "metadata_16S.txt")
mdit = read_tsv(ROOT / "data" / "metadata" / "metadata_ITS.txt")
add("metadata rows", "PASS" if (len(md16), len(mdit)) == (3652, 2310) else "FAIL",
    f"16S={len(md16)}, ITS={len(mdit)}")

# 2. adonis2 exact numeric comparison
new = read_tsv(ROOT / "results" / "tables" / "adonis2_results_3regions.tsv")
old = read_tsv(OLD / "result" / "Host_Seed_microbiota" / "result_table" / "adonis2_results_3regions.tsv")
maxdiff = 0.0
for a, b in zip(new, old):
    for col in ("Df", "SumOfSqs", "R2", "F", "P"):
        if a.get(col) in ("", "NA", None) or b.get(col) in ("", "NA", None):
            continue
        maxdiff = max(maxdiff, abs(num(a[col]) - num(b[col])))
add("adonis2 numeric", "PASS" if len(new) == len(old) and maxdiff < 1e-10 else "FAIL",
    f"rows={len(new)}/{len(old)}, max_abs_diff={maxdiff:.3g}")

# 3. phylosig estimates (permutation P only checked loosely)
newp = {(r["index"], r["region"]): r for r in read_tsv(ROOT / "results" / "tables" / "phylosig_results.tsv")}
oldp = {}
with (OLD / "result" / "Host_Seed_microbiota" / "result_table" / "phylosig_results.csv").open(encoding="utf-8-sig", newline="") as fh:
    for r in csv.DictReader(fh):
        oldp[(r["Alpha_Index"], r["Region"])] = r
max_est = 0.0
for k, r in newp.items():
    if k not in oldp:
        continue
    for a, b in (("K", "K"), ("lambda", "Lambda")):
        max_est = max(max_est, abs(num(r[a]) - num(oldp[k][b])))
add("phylosig estimates", "PASS" if max_est < 1e-10 else "FAIL",
    f"max_abs_diff={max_est:.3g}; permutation P values use fixed seed and may differ within Monte Carlo tolerance")

# 4. Mantel/MRM/dbRDA estimates
newm = {(r["Region"]): r for r in read_tsv(ROOT / "results" / "tables" / "mantel_mrm_results.tsv")}
with (OLD / "result" / "Host_Seed_microbiota" / "result_table" / "mantel_mrm_results.csv").open(encoding="utf-8-sig", newline="") as fh:
    oldm = {(r["Region"]): r for r in csv.DictReader(fh)}
max_m = 0.0
for k, r in newm.items():
    if k not in oldm:
        continue
    for col in ("Mantel_r", "MRM_R2", "dbRDA_R2"):
        max_m = max(max_m, abs(num(r[col]) - num(oldm[k][col])))
add("Mantel/MRM/dbRDA estimates", "PASS" if max_m < 1e-10 else "FAIL",
    f"max_abs_diff={max_m:.3g}; permutation P values may vary within Monte Carlo tolerance")

# 5. GOS classification counts
for ds, expected in (("16s", {"NRH": 12190, "generalist": 1095, "intermediate": 3813, "specialist": 5747}),
                     ("its", {"NRH": 2860, "generalist": 234, "intermediate": 522, "specialist": 1302})):
    counts = Counter(r["category"] for r in read_tsv(ROOT / "results" / "tables" / "gos" / f"gos_classified_{ds}.txt"))
    add(f"GOS counts {ds}", "PASS" if dict(counts) == expected else "FAIL", str(dict(counts)))

# 6. Sensitivity sets
rob = read_tsv(ROOT / "server_result" / "function_prediction" / "result_tables" / "sensitivity_robust17.tsv")
dis = read_tsv(ROOT / "server_result" / "function_prediction" / "result_tables" / "sensitivity_discordant9.tsv")
add("sensitivity robust17", "PASS" if len(rob) == 17 else "FAIL", f"n={len(rob)}")
add("sensitivity discordant9", "PASS" if len(dis) == 9 else "FAIL", f"n={len(dis)}")

# 7. Output map
with (ROOT / "manifests" / "figure_map.tsv").open(encoding="utf-8", newline="") as fh:
    fmap = list(csv.DictReader(fh, delimiter="\t"))
add("subplot map", "PASS" if len(fmap) == 83 else "FAIL", f"rows={len(fmap)}")

REPORT.parent.mkdir(parents=True, exist_ok=True)
with REPORT.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["check", "status", "detail"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

for r in rows:
    print(f"{r['status']:4} {r['check']}: {r['detail']}")
print(f"\nWrote {REPORT}")
raise SystemExit(0 if all(r["status"] == "PASS" for r in rows) else 1)
