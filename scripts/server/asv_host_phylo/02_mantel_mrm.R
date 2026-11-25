#!/usr/bin/env Rscript
# =============================================================================
# 02_mantel_mrm.R — main ASV-phylogeny vs host-distribution tests
#
# For each region (v4, v5_v7, its1), main filter (20/5/3):
#   - Mantel: Jaccard / Bray distribution distance ~ ASV patristic distance
#       * standard ASV-label permutation (one-sided greater)
#       * prevalence-stratified permutation (5 prevalence strata)
#   - MRM: distribution distance ~ phylo + |dlog prevalence| + |dlog abundance|
#   - Mantel correlogram: 10 phylo-distance bins, prevalence-stratified null
#
# Outputs:
#   result_table/mantel_results.tsv
#   result_table/mrm_results.tsv
#   result_table/correlogram.tsv
#   Temp_data/rds/mantel_nulls_main.rds
#   Temp_data/rds/correlogram_nulls_main.rds
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ape); library(vegan); library(ecodist)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/asv_host_phylo_utils.R"))

NPERM_MAIN  <- as.integer(Sys.getenv("NPERM_MAIN", "9999"))
NPERM_STRAT <- as.integer(Sys.getenv("NPERM_STRAT", "999"))
NPERM_MRM   <- as.integer(Sys.getenv("NPERM_MRM", "1999"))
NPERM_CG    <- as.integer(Sys.getenv("NPERM_CG", "999"))

regions <- c("v4", "v5_v7", "its1")
filter  <- "main20_5_3"
metric_map <- c(J = "dj", B = "db")

mantel_rows <- list(); mrm_rows <- list(); cg_rows <- list()
mantel_nulls <- list(); cg_nulls <- list()

for (ri in seq_along(regions)) {
  reg <- regions[ri]
  cat(sprintf("\n===== %s / %s =====\n", reg, filter))
  flush.console()
  p <- readRDS(profile_file(reg, filter))
  strata <- prevalence_strata(p$prev, 5)
  cat(sprintf("ASVs=%d; host-level=%d; prevalence strata sizes: %s\n",
              length(p$ids), length(p$host_ids),
              paste(tabulate(strata), collapse = "/")))

  for (mi in seq_along(metric_map)) {
    metric <- names(metric_map)[mi]
    d <- p[[metric_map[mi]]]
    key <- paste(reg, metric, sep = "_")
    cat(sprintf("  [%s] Mantel standard (%d perms, ecodist) ...\n",
                metric, NPERM_MAIN))
    flush.console()
    dp <- p$dp
    e_std <- ecodist::mantel(d ~ dp, nperm = NPERM_MAIN)
    # Custom 999-perm run only for the null mean/SD and SES (same statistic).
    m_ses <- mantel_fast(d, dp, nperm = 999L,
                         seed = 20261111L + 100L * ri + mi,
                         alternative = "greater")
    if (abs(as.numeric(e_std["mantelr"]) - m_ses$r) > 1e-8)
      warning(sprintf("[%s/%s] ecodist and custom Mantel r differ", reg, metric))
    cat(sprintf("  [%s] Mantel prevalence-stratified (%d perms) ...\n",
                metric, NPERM_STRAT))
    m_str <- mantel_fast(d, dp, nperm = NPERM_STRAT, strata = strata,
                         seed = 20261121L + 100L * ri + mi,
                         alternative = "greater")
    mantel_nulls[[key]] <- list(standard = m_ses$null, stratified = m_str$null)

    mantel_rows[[length(mantel_rows) + 1L]] <- data.table::data.table(
      region = reg, filter = filter, metric = metric,
      n_asv = length(p$ids), n_pairs = m_ses$n_pairs,
      r = as.numeric(e_std["mantelr"]),
      p_standard = as.numeric(e_std["pval1"]),
      p_standard_two_sided = as.numeric(e_std["pval3"]),
      CI_low = as.numeric(e_std["llim.2.5%"]),
      CI_high = as.numeric(e_std["ulim.97.5%"]),
      null_standard_mean = m_ses$null_mean,
      null_standard_sd = m_ses$null_sd,
      SES_standard = m_ses$SES,
      p_stratified = m_str$p_greater,
      null_stratified_mean = m_str$null_mean,
      null_stratified_sd = m_str$null_sd,
      SES_stratified = m_str$SES
    )

    cat(sprintf("  [%s] MRM (%d perms) ...\n", metric, NPERM_MRM))
    flush.console()
    d <- p[[metric_map[mi]]]
    dp <- p$dp; dpv <- p$dpv; dab <- p$dab
    mrm <- ecodist::MRM(d ~ dp + dpv + dab, nperm = NPERM_MRM)
    co <- mrm$coef
    mrm_rows[[length(mrm_rows) + 1L]] <- data.table::data.table(
      region = reg, filter = filter, metric = metric,
      n_asv = length(p$ids), n_pairs = length(d),
      term = rownames(co),
      coef = as.numeric(co[, 1L]),
      p = as.numeric(co[, 2L]),
      r_squared = as.numeric(mrm$r.squared[1L]),
      r_squared_p = as.numeric(mrm$r.squared[2L])
    )
    cat(sprintf("  [%s] MRM R2=%.4f; phylo coef=%.4g p=%.4g\n",
                metric, mrm$r.squared[1L],
                co["dp", 1L], co["dp", 2L]))
  }

  cat(sprintf("  correlogram (%d perms, stratified) ...\n", NPERM_CG))
  flush.console()
  cg <- correlogram_fast(list(J = p$dj, B = p$db), p$dp,
                         nbin = 10L, nperm = NPERM_CG,
                         strata = strata, seed = 20261131L + ri)
  for (metric in names(cg)) {
    cg_rows[[length(cg_rows) + 1L]] <-
      cbind(data.table::data.table(region = reg, filter = filter), cg[[metric]])
    cg_nulls[[paste(reg, metric, sep = "_")]] <- attr(cg[[metric]], "null_draws")
  }
  rm(p); invisible(gc())
}

mantel_res <- data.table::rbindlist(mantel_rows)
mantel_res[, p_adj_BH := stats::p.adjust(p_standard, method = "BH")]
mrm_res <- data.table::rbindlist(mrm_rows)
mrm_res[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = .(metric)]
cg_res <- data.table::rbindlist(cg_rows, fill = TRUE)

fwrite_safe(mantel_res, file.path(OUT_ROOT, "result_table", "mantel_results.tsv"))
fwrite_safe(mrm_res, file.path(OUT_ROOT, "result_table", "mrm_results.tsv"))
fwrite_safe(cg_res, file.path(OUT_ROOT, "result_table", "correlogram.tsv"))
saveRDS(mantel_nulls, file.path(OUT_ROOT, "Temp_data/rds", "mantel_nulls_main.rds"))
saveRDS(cg_nulls, file.path(OUT_ROOT, "Temp_data/rds", "correlogram_nulls_main.rds"))

cat("\n=== Mantel results ===\n"); print(mantel_res)
flush.console()
cat("\n=== MRM phylo terms ===\n")
flush.console()
print(mrm_res[term == "dp"])
cat("\nDONE 02_mantel_mrm\n")
