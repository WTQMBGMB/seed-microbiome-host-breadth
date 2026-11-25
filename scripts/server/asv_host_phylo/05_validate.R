#!/usr/bin/env Rscript
# =============================================================================
# 05_validate.R — consistency checks for the ASV-host-phylogeny outputs
#
# - expected row counts and no missing key statistics
# - recompute the main Mantel r with vegan::mantel and compare to the table
# - check tree/profile dimensions
# - check figure PDF page sizes
# Output: result_table/validation_summary.tsv
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ape)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/asv_host_phylo_utils.R"))

rows <- list()
add <- function(check, value, pass, detail = "") {
  rows[[length(rows) + 1L]] <<- data.table::data.table(
    check = check, value = as.character(value), pass = isTRUE(pass), detail = detail)
}

tab_dir <- file.path(OUT_ROOT, "result_table")
mt  <- fread(file.path(tab_dir, "mantel_results.tsv"))
mrm <- fread(file.path(tab_dir, "mrm_results.tsv"))
cg  <- fread(file.path(tab_dir, "correlogram.tsv"))
sens <- fread(file.path(tab_dir, "sensitivity_results.tsv"))
cp  <- fread(file.path(tab_dir, "close_pair_contrast.tsv"))
ext <- fread(file.path(tab_dir, "mrm_extended_results.tsv"))

add("mantel_rows", nrow(mt), nrow(mt) == 6L)
add("mantel_no_na", sum(is.na(mt$r) | is.na(mt$p_standard)), all(!is.na(mt$r) & !is.na(mt$p_standard)))
add("mrm_rows", nrow(mrm), nrow(mrm) == 24L, "6 region-metric x 4 terms")
add("mrm_phylo_rows", sum(mrm$term == "dp"), sum(mrm$term == "dp") == 6L)
add("correlogram_rows", nrow(cg), nrow(cg) == 60L, "3 regions x 2 metrics x 10 bins")
add("sensitivity_rows", nrow(sens), nrow(sens) >= 12L)
add("close_pair_rows", nrow(cp), nrow(cp) == 6L)
add("extended_mrm_rows", nrow(ext), nrow(ext) >= 18L)

# Recompute main Mantel r with vegan::mantel (999 permutations only for the r;
# the p in the table comes from the 9,999-perm custom test).
for (i in seq_len(nrow(mt))) {
  reg <- mt$region[i]; metric <- mt$metric[i]
  p <- readRDS(profile_file(reg, "main20_5_3"))
  d <- p[[if (metric == "J") "dj" else "db"]]
  vm <- vegan::mantel(d, p$dp, method = "pearson", permutations = 99)
  add(sprintf("vegan_r_%s_%s", reg, metric),
      sprintf("%.8f vs %.8f", vm$statistic, mt$r[i]),
      abs(vm$statistic - mt$r[i]) < 1e-8)
  rm(p); invisible(gc())
}

# Tree/profile dimensions
qc <- fread(file.path(tab_dir, "tree_qc.tsv"))
add("tree_qc_all_match", sum(!qc$all_candidates_in_tree), all(qc$all_candidates_in_tree))
add("tree_qc_no_negative", sum(qc$n_neg_brlen), all(qc$n_neg_brlen == 0))

# Figure page sizes via pdftools (if available)
figs <- c("FigA_candidate_filtering.pdf", "FigB_mantel_scatter.pdf",
          "FigC_correlogram.pdf", "FigD_sensitivity_forest.pdf",
          "FigE_close_pair_contrast.pdf")
if (requireNamespace("pdftools", quietly = TRUE)) {
  for (f in figs) {
    fp <- file.path(OUT_ROOT, "result_figure", f)
    if (file.exists(fp)) {
      info <- pdftools::pdf_info(fp)
      add(sprintf("pdf_size_%s", sub("\\.pdf$", "", f)),
          sprintf("%.0fx%.0f pt", info$width[1], info$height[1]),
          info$width[1] > 100 && info$height[1] > 100)
    }
  }
}

res <- data.table::rbindlist(rows)
data.table::fwrite(res, file.path(tab_dir, "validation_summary.tsv"), sep = "\t")
print(res)
cat("DONE 05_validate\n")
