# -*- coding: utf-8 -*-
"""Stage precomputed expensive results into server_result/."""
from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path

SRC = Path(r"D:\seed_meta_new\Downstream_analysis")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "server_result"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_file(src: Path, dst: Path, note: str, rows: list[dict]) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists() or dst.stat().st_size != src.stat().st_size:
        shutil.copy2(src, dst)
    rows.append({
        "role": dst.relative_to(OUT).parts[0],
        "path": dst.relative_to(ROOT).as_posix(),
        "source_path": str(src),
        "bytes": dst.stat().st_size,
        "sha256": sha256(dst),
        "notes": note,
    })


def copy_dir(src: Path, dst: Path, note: str, rows: list[dict], patterns: tuple[str, ...] = ("*",)) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    for pattern in patterns:
        for f in sorted(src.rglob(pattern)):
            if not f.is_file():
                continue
            rel = f.relative_to(src)
            copy_file(f, dst / rel, note, rows)


def main() -> None:
    rows: list[dict] = []
    rds = SRC / "result" / "Host_Seed_microbiota" / "rds"
    tab = SRC / "result" / "Host_Seed_microbiota" / "result_table"

    for reg in ["v3_v4", "v4", "v5_v6", "v5_v7", "its1", "its2"]:
        copy_file(rds / f"bray_dist_{reg}.rds", OUT / "beta" / f"bray_dist_{reg}.rds",
                  "Bray-Curtis distance (server/archived)", rows)
    for reg in ["v4", "v5_v7", "its1"]:
        copy_file(rds / "adonis2_source" / f"adonis_Genus_{reg}.rds",
                  OUT / "beta" / "adonis2_source" / f"adonis_Genus_{reg}.rds",
                  "adonis2 999-permutation checkpoint", rows)
    for name in ["boot_lambda_samples.rds", "phylo_perm_nulls.rds", "phylo_perm_null.rds"]:
        copy_file(rds / name, OUT / "phylo" / name, "phylosig bootstrap / permutation checkpoint", rows)
    for name in ["phylo_perm_alpha.tsv", "phylo_perm_beta.tsv", "phylo_perm_summary.tsv",
                 "phylosig_lambda_ci.tsv"]:
        copy_file(tab / name, OUT / "phylo" / name, "phylosig permutation/bootstrap summary", rows)

    bnti_in = tab / "bnti" / "inputs"
    copy_dir(bnti_in, OUT / "bnti" / "inputs", "bNTI regional inputs", rows)
    bnti_rds = rds / "bnti_local"
    copy_dir(bnti_rds, OUT / "bnti", "bNTI local pool/prep checkpoint", rows, ("prep_*.rds", "pool_*.rds"))
    copy_dir(tab / "bnti" / "local", OUT / "bnti" / "local", "bNTI local pair tables and summaries", rows)

    gos_tmp = SRC / "result" / "GOS" / "Temp_data"
    for name in [
        "host_dist.rds", "presence_merged_16s.rds", "presence_merged_its.rds",
        "nhost_null_16s_iter1000.rds", "nhost_null_its_iter1000.rds",
        "mpd_null_16s_rep1000.rds", "mpd_null_its_rep1000.rds",
        "gos_threshold_16s.rds", "gos_threshold_its.rds",
        "presence_species_merge_log_16s.csv", "presence_species_merge_log_its.csv",
    ]:
        copy_file(gos_tmp / name, OUT / "gos" / name, "GOS null-model / presence checkpoint", rows)
    for name in ["presence_merged_16s.txt", "presence_merged_its.txt"]:
        copy_file(SRC / "result" / "GOS" / "result_table" / name, OUT / "gos" / name,
                  "GOS merged presence matrix", rows)
    for name in ["gos_phyla_five_column_sankey_16s_nodes.tsv", "gos_phyla_five_column_sankey_16s_links.tsv",
                 "gos_class_five_column_sankey_its_nodes.tsv", "gos_class_five_column_sankey_its_links.tsv",
                 "phylum_color_map_16s.tsv",
                 "taxon_hbs_results_16s.tsv", "taxon_hbs_results_its.tsv",
                 "taxon_hbs_summary_16s.tsv", "taxon_hbs_summary_its.tsv",
                 "genus_logistic_twopart_16s.tsv", "genus_logistic_twopart_its.tsv",
                 "logistic_3x3_counts_16s.tsv", "logistic_3x3_counts_its.tsv"]:
        copy_file(SRC / "result" / "GOS" / "result_table" / name,
                  OUT / "gos" / "precomputed" / name, "GOS downstream result table", rows)

    tree = SRC / "result" / "GOS" / "tree_16s_allASV"
    for name in [
        "gos_all_hostn_16s_midpoint.nwk", "all_asv_tree_itol.nwk",
        "annotation_16s_allASV.tsv", "identical_sequence_groups.tsv",
        "phylosig_hostn_16s_allASV_all.tsv", "phylosig_hostn_16s_allASV_positive.tsv",
        "phyloD_generalist_all.tsv", "phyloD_generalist_positive.tsv",
    ]:
        copy_file(tree / name, OUT / "gos" / "tree_16s_allASV" / name, "all-ASV tree / phylosig server result", rows)
    copy_dir(tree / "itol_annotation", OUT / "gos" / "tree_16s_allASV" / "itol_annotation",
             "iTOL annotation dataset", rows)

    asv = SRC / "result" / "ASV_host_phylo"
    for sub in ["Temp_data/trees", "Temp_data/tree_input", "Temp_data/candidates", "Temp_data/rds", "result_table", "paper_main"]:
        copy_dir(asv / sub, OUT / "asv_beta_phylo" / sub, "ASV-host phylo server result", rows)

    for name in ["Table_S7A_host_association.tsv", "Table_S7B_alpha_phylogeny.tsv",
                 "Table_S7C_community_phylogeny.tsv", "Table_S7D_dispersion.tsv",
                 "Table_S7_notes.tsv"]:
        copy_file(tab / "sensitivity" / name, OUT / "sensitivity" / name, "sensitivity table source", rows)

    for name in ["beta_source_R2_decomposition.tsv", "beta_source_model.tsv", "beta_source_vpa3.tsv",
                 "beta_source_structure.tsv", "beta_source_robustness.tsv", "beta_source_tests.tsv"]:
        copy_file(tab / name, OUT / "beta" / "source_tables" / name,
                  "beta source decomposition result table", rows)

    trtab = SRC / "result" / "traits_prevalence" / "result_table"
    for name in ["results_numeric_16s.tsv", "results_boolean_16s.tsv", "results_boolean_twopart_16s.tsv",
                 "results_ko_twopart_16s.tsv", "results_ko_twopart_excluded_16s.tsv",
                 "twopart_counts_boolean_16s.tsv", "twopart_counts_ko_16s.tsv",
                 "sensitivity_all_traits.tsv", "sensitivity_robust17.tsv", "sensitivity_discordant9.tsv",
                 "ko_sankey_kegg_lv2_16s_abs_links.tsv", "ko_sankey_kegg_lv2_16s_abs_nodes.tsv",
                 "ko_sankey_kegg_lv2_16s_links.tsv", "ko_sankey_kegg_lv2_16s_nodes.tsv",
                 "ko_sankey_pgpg_lv3_16s_abs_links.tsv", "ko_sankey_pgpg_lv3_16s_abs_nodes.tsv",
                 "ko_sankey_pgpg_lv3_16s_links.tsv", "ko_sankey_pgpg_lv3_16s_nodes.tsv"]:
        copy_file(trtab / name, OUT / "function_prediction" / "result_tables" / name,
                  "trait/KO result table", rows)
    copy_file(SRC / "result" / "taxonomy_16s_traits.txt",
              OUT / "function_prediction" / "taxonomy_16s_traits.txt", "16S trait profile", rows)

    boot = SRC / "result" / "traits_prevalence" / "source_bootstrap" / "independent_20261110"
    copy_dir(boot, OUT / "function_prediction" / "source_bootstrap", "source-level bootstrap", rows)

    manifest = OUT / "MANIFEST.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    total = sum(r["bytes"] for r in rows)
    print(f"Wrote {manifest}; {len(rows)} files; {total / 1024**2:.1f} MiB")


if __name__ == "__main__":
    main()
