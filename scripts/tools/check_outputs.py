# -*- coding: utf-8 -*-
"""Check that every expected figure and table was generated."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "manifests" / "output_manifest.tsv"
RESULT = ROOT / "result"


def main() -> int:
    if not MANIFEST.exists():
        print("missing manifest:", MANIFEST)
        return 1
    with MANIFEST.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    errors: list[str] = []
    for row in rows:
        path = RESULT / row["path"]
        if not path.exists():
            errors.append(f"missing: {row['path']}")
        elif path.stat().st_size == 0:
            errors.append(f"empty: {row['path']}")

    expected = {row["path"] for row in rows}
    extras = []
    for ch in ["ch1_dataset", "ch2_host_effect", "ch3_phylogeny",
               "ch4_generality", "ch5_traits"]:
        for kind, pattern in (("figures", "*.pdf"), ("figures", "*.jpg"),
                              ("tables", "*.tsv")):
            for path in (RESULT / ch / kind).glob(pattern):
                rel = path.relative_to(RESULT).as_posix()
                if rel not in expected:
                    extras.append(rel)
    if extras:
        errors.append("unexpected outputs: " + ", ".join(sorted(extras)))

    if errors:
        print("OUTPUT CHECK FAILED")
        print("\n".join(errors))
        return 1

    n_fig = sum(row["kind"] == "figure" for row in rows)
    n_tab = sum(row["kind"] == "table" for row in rows)
    print(f"OUTPUT CHECK PASS: {n_fig} figures, {n_tab} tables")
    return 0


if __name__ == "__main__":
    sys.exit(main())
