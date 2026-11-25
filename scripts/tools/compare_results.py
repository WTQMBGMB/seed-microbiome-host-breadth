# -*- coding: utf-8 -*-
"""Compare figures and tables between a reference result tree and a new run.

Usage / 用法:
  python scripts/tools/compare_results.py --ref <reference/result> \
      --new result [--rtol 1e-6]
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

CHAPTERS = ["ch1_dataset", "ch2_host_effect", "ch3_phylogeny",
            "ch4_generality", "ch5_traits"]


def read_tsv(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def as_float(x: str):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def same_cell(a: str, b: str, rtol: float, atol: float) -> bool:
    if a == b:
        return True
    fa, fb = as_float(a), as_float(b)
    if fa is None or fb is None:
        return False
    return math.isclose(fa, fb, rel_tol=rtol, abs_tol=atol)


def compare_table(ref: Path, new: Path, rtol: float, atol: float) -> list[str]:
    a = read_tsv(ref)
    b = read_tsv(new)
    if not a or not b:
        return ["empty table"]
    if a[0] != b[0]:
        return ["header mismatch"]
    if len(a) != len(b):
        return [f"row count {len(a)} vs {len(b)}"]
    mismatches: list[str] = []
    for i, (ra, rb) in enumerate(zip(a[1:], b[1:]), start=2):
        if len(ra) != len(rb):
            mismatches.append(f"row {i}: column count")
            continue
        for j, (ca, cb) in enumerate(zip(ra, rb)):
            if not same_cell(ca, cb, rtol, atol):
                mismatches.append(f"row {i} col {a[0][j]}: {ca!r} vs {cb!r}")
                break
        if len(mismatches) >= 20:
            mismatches.append("... more than 20 mismatching rows")
            break
    return mismatches


def compare_figures(ref: Path, new: Path) -> list[str]:
    errors = []
    for p in sorted(ref.glob("*")):
        q = new / p.name
        if not q.exists():
            errors.append(f"missing figure: {p.name}")
        elif q.stat().st_size == 0:
            errors.append(f"empty figure: {p.name}")
    return errors


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="reference result directory")
    ap.add_argument("--new", default="result", help="new result directory")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-8)
    args = ap.parse_args()

    ref_root = Path(args.ref)
    new_root = Path(args.new)
    failures = 0

    for ch in CHAPTERS:
        for kind in ("figures", "tables"):
            rdir, ndir = ref_root / ch / kind, new_root / ch / kind
            if not rdir.exists():
                continue
            ref_files = {p.name for p in rdir.glob("*") if p.is_file()}
            new_files = {p.name for p in ndir.glob("*") if p.is_file()} \
                if ndir.exists() else set()
            missing = sorted(ref_files - new_files)
            extra = sorted(new_files - ref_files)
            if missing:
                print(f"[{ch}/{kind}] missing: {', '.join(missing)}")
                failures += len(missing)
            if extra and kind == "figures":
                print(f"[{ch}/{kind}] extra: {', '.join(extra)}")
            if kind == "figures":
                for err in compare_figures(rdir, ndir):
                    print(f"[{ch}/figures] {err}")
                    failures += 1
            else:
                for name in sorted(ref_files & new_files):
                    if not name.endswith(".tsv"):
                        continue
                    diff = compare_table(rdir / name, ndir / name,
                                         args.rtol, args.atol)
                    if diff:
                        print(f"[{ch}/tables] {name}: {'; '.join(diff[:3])}")
                        failures += 1
    if failures:
        print(f"COMPARE FAILED: {failures} issue(s)")
        return 1
    print("COMPARE PASS: file sets and table values match")
    return 0


if __name__ == "__main__":
    sys.exit(main())
