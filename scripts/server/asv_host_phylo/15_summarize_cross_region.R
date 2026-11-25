#!/usr/bin/env Rscript
# Summarize cross-region mapped-bin results into a markdown report.

suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("BASE_DIR", "<SERVER_BASE>")
TAB <- file.path(BASE, "result/ASV_host_phylo/result_table/cross_region")
OUT <- file.path(TAB, "RESULTS_cross_region.md")

read_if <- function(f) if (file.exists(f)) fread(f) else NULL
fmt <- function(x, d = 3) {
  if (is.numeric(x)) formatC(x, format = "f", digits = d) else as.character(x)
}
md_table <- function(dt, cols = names(dt), digits = 3) {
  if (is.null(dt) || !nrow(dt)) return("_No rows._\n")
  dt <- dt[, ..cols]
  for (nm in names(dt)) {
    if (!is.numeric(dt[[nm]])) next
    x <- dt[[nm]]
    if (all(is.na(x) | x == round(x))) {
      dt[[nm]] <- as.character(as.integer(x))
    } else {
      dt[[nm]] <- fmt(x, digits)
    }
  }
  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::kable(dt, format = "pipe", row.names = FALSE)
  } else {
    paste(capture.output(print(dt, row.names = FALSE)), collapse = "\n")
  }
}

t1_16 <- read_if(file.path(TAB, "track1_16s_mantel.tsv"))
t1_it <- read_if(file.path(TAB, "track1_its_mantel.tsv"))
t1 <- rbindlist(list(t1_16, t1_it), fill = TRUE)
t1m <- rbindlist(list(
  read_if(file.path(TAB, "track1_16s_mrm.tsv")),
  read_if(file.path(TAB, "track1_its_mrm.tsv"))), fill = TRUE)
t1meta <- rbindlist(list(
  read_if(file.path(TAB, "track1_16s_meta.tsv")),
  read_if(file.path(TAB, "track1_its_meta.tsv"))), fill = TRUE)
t1cf <- rbindlist(list(
  read_if(file.path(TAB, "track1_16s_close_far.tsv")),
  read_if(file.path(TAB, "track1_its_close_far.tsv"))), fill = TRUE)

t2 <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_mantel.tsv")),
  read_if(file.path(TAB, "track2_its_mantel.tsv"))), fill = TRUE)
t2p <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_pairs.tsv")),
  read_if(file.path(TAB, "track2_its_pairs.tsv"))), fill = TRUE)
t2m <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_mrm.tsv")),
  read_if(file.path(TAB, "track2_its_mrm.tsv"))), fill = TRUE)
t2s <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_summary.tsv")),
  read_if(file.path(TAB, "track2_its_summary.tsv"))), fill = TRUE)
t2cf <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_close_far.tsv")),
  read_if(file.path(TAB, "track2_its_close_far.tsv"))), fill = TRUE)
t2sens <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_sensitivity.tsv")),
  read_if(file.path(TAB, "track2_its_sensitivity.tsv"))), fill = TRUE)
if (!is.null(t1m) && !"extra_terms" %in% names(t1m))
  t1m[, extra_terms := NA_character_]
if (!is.null(t2m) && !"extra_terms" %in% names(t2m))
  t2m[, extra_terms := NA_character_]

lines <- c(
  "# Cross-region mapped-reference-bin analysis",
  "",
  "Unit: rows of the cross-region merged tables whose ID is present in",
  "`Merged_16S_Silva_filtered.fa` or `Merged_ITS_Unite_filtered.fa`",
  "(vsearch `--usearch_global --id 0.97`; one row = one 97% reference bin,",
  "not an exact ASV). 16S uses the existing 22,845-tip SILVA reference tree;",
  "ITS uses a 4,918-tip UNITE reference tree built for this analysis.",
  "",
  "## Track 1 — within-region mapped-bin replication",
  "",
  "Filters per region: prevalence >=20, host species >=5 (adaptive to the",
  "region's species pool), sources >=3. Jaccard/Bray distribution distances",
  "vs reference-tree patristic distance. Mantel p is one-sided; p_stratified",
  "uses prevalence-stratified permutations. MRM controls |delta log prevalence|",
  "and |delta log mean relative abundance|; extended MRM adds host breadth and",
  "source-set distance.",
  "",
  "### Mantel",
  md_table(t1[, .(dataset, region, metric, n_bins, r, CI_low, CI_high,
                  p_standard, p_stratified, SES)]),
  "",
  "### MRM, phylo term",
  md_table(t1m[term == "dp", .(dataset, region, metric, model, coef, p,
                                r_squared)]),
  "",
  "### Meta-analysis of Mantel r (Fisher z)",
  md_table(t1meta[, .(dataset, metric, region_set, k, r_RE, CI_low_RE,
                       CI_high_RE, p_RE, I2, Q_p)]),
  "",
  "### Short-range contrast",
  md_table(t1cf[, .(dataset, region, metric, contrast, p, SES, n_close, n_far)]),
  "",
  "## Track 2 — cross-region shared-bin analysis",
  "",
  "For each region pair, only bins detected in both regions; primary filter",
  "total prevalence >=20, total hosts >=5, total sources >=3, and >=5 positive",
  "samples in each region. Within-region Jaccard/Bray distances are averaged",
  "with equal weight; Mantel/MRM use prevalence-stratified permutations.",
  "",
  "### Pair counts",
  md_table(t2p),
  "",
  "### Mantel",
  md_table(t2[, .(dataset, region_a, region_b, pair_type, metric, n_bins, r,
                  CI_low, CI_high, p_standard, p_stratified, SES)]),
  "",
  "### MRM, phylo term",
  md_table(t2m[term == "dp", .(dataset, region_a, region_b, pair_type, metric,
                                model, coef, p, r_squared, extra_terms)]),
  "",
  "### Pair summary (descriptive; pairs are not independent)",
  md_table(t2s),
  "",
  "### Short-range contrast",
  md_table(t2cf[, .(dataset, region_a, region_b, pair_type, metric, contrast,
                     p, SES)]),
  "",
  "## Caveats",
  "",
  "- Mapped bins are 97% reference bins; multiple ASVs from different regions",
  "  and studies can collapse into one row. This is not exact-ASV inference.",
  "- Cross-region shared bins are a selected subset (common/generalist taxa);",
  "  results generalize to cross-region-detectable bins, not all microbes.",
  "- V5_V6 has only 4 host species and 4 sources; the 20/5/3 filter is",
  "  impossible there, so its Track 1 result uses the adaptive 4/3 threshold",
  "  and is exploratory.",
  "- V3_V4/V4 and V5_V6/V5_V7 share overlapping amplicon segments; they are",
  "  reported separately and are not independent evidence.",
  "- Region pairs share bins and regions; no formal meta-analysis across",
  "  Track 2 pairs is performed.",
  "- Cross-region zeros can reflect primer/detection differences, not true",
  "  absence; within-region distances and region-balanced pair distances",
  "  reduce, but do not eliminate, this bias."
)
writeLines(lines, OUT)
cat("Written:", OUT, "\n")
