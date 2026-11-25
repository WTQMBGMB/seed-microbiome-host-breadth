# -*- coding: utf-8 -*-
"""Download and verify the archived input data.

The archive URL and SHA256 are listed in manifests/data_download.tsv.
归档地址与校验和见 manifests/data_download.tsv。

Usage / 用法:
  python scripts/tools/download_data.py
  python scripts/tools/download_data.py --manifest manifests/data_download.tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def extract(archive: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tf:
        tf.extractall(dest, filter="data")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=str(ROOT / "manifests" / "data_download.tsv"))
    args = ap.parse_args()
    manifest = Path(args.manifest)
    if not manifest.exists():
        print("missing manifest:", manifest)
        return 1

    with manifest.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    ok = 0
    for row in rows:
        url = (row.get("url") or "").strip()
        expected = (row.get("sha256") or "").strip().lower()
        name = (row.get("archive") or "data_archive").strip()
        if not url or url.startswith("<"):
            print(f"skip {name}: URL not configured in {manifest}")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / name
            print(f"downloading {url}")
            urllib.request.urlretrieve(url, archive)
            got = sha256(archive)
            if got != expected:
                print(f"sha256 mismatch for {name}\n  expected {expected}\n  got      {got}")
                return 1
            print(f"verified {name}: {got}")
            extract(archive, ROOT)
        ok += 1

    if not ok:
        print("No archive downloaded. Edit manifests/data_download.tsv first.")
        return 1
    print(f"data and server results ready under {ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
