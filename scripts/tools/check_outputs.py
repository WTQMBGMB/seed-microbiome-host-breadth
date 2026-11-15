# -*- coding: utf-8 -*-
"""Check that every Material_Description2 subplot has a mapped output."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAP = ROOT / "manifests" / "figure_map.tsv"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"


def main() -> int:
    with MAP.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    errors = []
    for r in rows:
        expected = r["Expected_file"]
        if r["Output_kind"] == "figure":
            stem = expected[:-4]
            for ext in (".pdf", ".jpg"):
                p = FIG / (stem + ext)
                if not p.exists() or p.stat().st_size == 0:
                    errors.append(f"missing figure: {p}")
        else:
            p = TAB / expected
            if not p.exists() or p.stat().st_size == 0:
                errors.append(f"missing table: {p}")
    figs = {p.name for p in FIG.glob("*.pdf")}
    mapped = {r["Expected_file"] for r in rows if r["Output_kind"] == "figure"}
    extra = sorted(f for f in figs if f not in mapped)
    if extra:
        errors.append("unmapped figure PDFs: " + ", ".join(extra))
    if errors:
        print("OUTPUT CHECK FAILED")
        print("\n".join(errors))
        return 1
    print(f"OUTPUT CHECK PASS: {len(rows)} subplots, {len(figs)} figure PDFs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
