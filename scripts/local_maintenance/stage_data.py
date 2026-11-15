# -*- coding: utf-8 -*-
"""Stage the trimmed input data for github2.

The script copies only files that are actually read by the retained
analysis chapters.  It never writes to the source project.
"""
from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path

SRC = Path(r"D:\seed_meta_new\Downstream_analysis")
BRIDGE = Path(r"D:\seed_meta_new\Paper_writing2")
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


def files_for_region(reg: str) -> list[tuple[Path, Path, str]]:
    src = SRC / "filtered_region" / reg
    dst = DATA / "filtered_region" / reg
    out = []
    for suffix, note in [
        ("merged_asvtab_filtered.txt", "ASV abundance table"),
        ("merged_taxonomy_filtered.txt", "ASV taxonomy"),
        ("merged_filtered.fa", "ASV sequences"),
    ]:
        name = f"{reg}_{suffix}"
        out.append((src / name, dst / name, note))
    for name in ["filtering_report.json", "subset_report.json"]:
        if (src / name).exists():
            out.append((src / name, dst / name, "filtering provenance"))
    return out


def build_pairs() -> list[tuple[Path, Path, str]]:
    fm = SRC / "filtered_merge"
    fp = SRC / "Function_prediction"
    pairs: list[tuple[Path, Path, str]] = [
        (SRC / "standard_merge" / "metadata_16S.txt", DATA / "metadata" / "metadata_16S.txt", "authoritative metadata"),
        (SRC / "standard_merge" / "metadata_ITS.txt", DATA / "metadata" / "metadata_ITS.txt", "authoritative metadata"),
        (fm / "16S" / "Merged_16S_table_filtered.txt", DATA / "filtered_merge" / "16S" / "Merged_16S_table_filtered.txt", "filtered merged ASV table"),
        (fm / "16S" / "Merged_16S_taxonomy_filtered.txt", DATA / "filtered_merge" / "16S" / "Merged_16S_taxonomy_filtered.txt", "filtered merged taxonomy"),
        (fm / "16S" / "Merged_16S_Silva_filtered.fa", DATA / "filtered_merge" / "16S" / "Merged_16S_Silva_filtered.fa", "mapped 16S sequences"),
        (fm / "ITS" / "Merged_ITS_table_filtered.txt", DATA / "filtered_merge" / "ITS" / "Merged_ITS_table_filtered.txt", "filtered merged ASV table"),
        (fm / "ITS" / "Merged_ITS_taxonomy_filtered.txt", DATA / "filtered_merge" / "ITS" / "Merged_ITS_taxonomy_filtered.txt", "filtered merged taxonomy"),
        (fm / "ITS" / "Merged_ITS_Unite_filtered.fa", DATA / "filtered_merge" / "ITS" / "Merged_ITS_Unite_filtered.fa", "mapped ITS sequences"),
        (fp / "FAPROTAX2_provenance.json", DATA / "function_prediction" / "FAPROTAX2_provenance.json", "prediction provenance"),
        (fp / "faprotax" / "faprotax_result.txt", DATA / "function_prediction" / "faprotax" / "faprotax_result.txt", "FAPROTAX2 result"),
        (fp / "bugbase" / "bugbase_result.txt", DATA / "function_prediction" / "bugbase" / "bugbase_result.txt", "BugBase result"),
        (fp / "metatraits" / "trait_annotation" / "asv_trait_matrix_clean_boolean.tsv", DATA / "function_prediction" / "metatraits" / "trait_annotation" / "asv_trait_matrix_clean_boolean.tsv", "metaTraits boolean matrix"),
        (fp / "metatraits" / "trait_annotation" / "asv_trait_matrix_clean_continuous.tsv", DATA / "function_prediction" / "metatraits" / "trait_annotation" / "asv_trait_matrix_clean_continuous.tsv", "metaTraits continuous matrix"),
        (fp / "picrust2" / "ko_data_nsti_filtered.txt", DATA / "function_prediction" / "picrust2" / "ko_data_nsti_filtered.txt", "PICRUSt2 KO matrix"),
        (SRC / "plant_tree" / "species_all.nwk", DATA / "plant_tree" / "species_all.nwk", "TimeTree5 export"),
        (SRC / "plant_tree" / "species_all.txt", DATA / "plant_tree" / "species_all.txt", "TimeTree5 species list"),
        (SRC / "db" / "ko_level_table.txt", DATA / "db" / "ko_level_table.txt", "KEGG KO hierarchy"),
        (SRC / "db" / "ko_pgpg_levels.txt", DATA / "db" / "ko_pgpg_levels.txt", "PGPG KO hierarchy"),
        (SRC / "db" / "FAPROTAX2" / "function_config.tsv", DATA / "db" / "FAPROTAX2" / "function_config.tsv", "FAPROTAX2 configuration"),
        (SRC / "db" / "FAPROTAX2" / "taxonomy138_input.tsv", DATA / "db" / "FAPROTAX2" / "taxonomy138_input.tsv", "FAPROTAX2 taxonomy input"),
        (SRC / "db" / "SILVA_138.2_taxonomy.txt", DATA / "db" / "SILVA_138.2_taxonomy.txt", "SILVA reference taxonomy"),
        (BRIDGE / "table" / "Table S1.xlsx", DATA / "manual" / "Table S1.xlsx", "manually curated article list"),
    ]
    for reg in ["V3_V4", "V4", "V5_V6", "V5_V7", "ITS1", "ITS2"]:
        pairs.extend(files_for_region(reg))
    gis = SRC / "scripts" / "python" / "ne_110m_admin_0_countries"
    for name in [
        "ne_110m_admin_0_countries.shp",
        "ne_110m_admin_0_countries.shx",
        "ne_110m_admin_0_countries.dbf",
        "ne_110m_admin_0_countries.prj",
        "ne_110m_admin_0_countries.cpg",
        "ne_110m_admin_0_countries.VERSION.txt",
        "ne_110m_admin_0_countries.README.html",
    ]:
        pairs.append((gis / name, DATA / "gis" / name, "Natural Earth map"))
    return pairs


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    pairs = build_pairs()
    rows = []
    total = 0
    for src, dst, note in pairs:
        if not src.exists():
            raise FileNotFoundError(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists() or dst.stat().st_size != src.stat().st_size:
            shutil.copy2(src, dst)
        total += dst.stat().st_size
        rows.append({
            "role": dst.relative_to(DATA).parts[0],
            "path": dst.relative_to(ROOT).as_posix(),
            "source_path": str(src),
            "bytes": dst.stat().st_size,
            "sha256": sha256(dst),
            "notes": note,
        })
        print(f"staged {dst.relative_to(ROOT)}")
    manifest = DATA / "MANIFEST.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {manifest}; {len(rows)} files; {total / 1024**3:.2f} GiB")


if __name__ == "__main__":
    main()
