#!/usr/bin/env Rscript
# Source-stratified permutation null for the paper main analysis.
#
# Mimics the logic of Fig_phylo_permutation_beta: keep the observed distribution
# profiles and tree fixed, shuffle bin labels within structured strata, and
# compare the observed Mantel r with 999 null draws.
#
# Strata (preserve study structure + prevalence):
#   Track 1: dominant source x prevalence quintile
#   Track 2: dominant source x total prevalence quintile (V4 + V5-V7)
#
# Outputs in paper_main/:
#   null_model_results.tsv
#   null_draws.tsv
#   null_draws.rds

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ape); library(vegan); library(ecodist)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/mapped_utils.R"))

NPERM <- as.integer(Sys.getenv("NPERM_NULL", "999"))
UNITS <- strsplit(Sys.getenv("NULL_UNITS", "V4,V5_V7,CROSS"), ",")[[1]]
TAB <- file.path(OUT_ROOT, "result_table", "cross_region")
OUT <- file.path(OUT_ROOT, "paper_main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cache <- load_mapped_cache("16S")
tree_file <- file.path(BASE_DIR, MAPPED_CFG[["16S"]]$tree)

dominant_source <- function(mat, meta) {
  so <- meta$source
  so[!nzchar(so)] <- NA_character_
  sm <- model.matrix(~ 0 + factor(so))
  counts <- (mat > 0) %*% sm
  dom <- colnames(counts)[max.col(counts, ties.method = "first")]
  sub("factor(so)", "", dom, fixed = TRUE)
}

paired_dist <- function(dA, dB) {
  ma <- as.matrix(dA); mb <- as.matrix(dB)
  mm <- (ma + mb) / 2
  stats::as.dist(mm)
}

existing_p <- rbindlist(list(
  fread(file.path(TAB, "track1_16s_mantel.tsv"))[region %in% c("V4", "V5_V7"),
    .(unit = fifelse(region == "V4", "V4 (within)", "V5-V7 (within)"),
      metric, p_prevalence = p_stratified)],
  fread(file.path(TAB, "track2_16s_mantel.tsv"))[
    region_a == "V4" & region_b == "V5_V7",
    .(unit = "V4-V5-V7 (shared)", metric, p_prevalence = p_stratified)]
), fill = TRUE)

make_within <- function(reg, unit_label) {
  rs <- region_subset(cache, reg)
  sel <- select_candidates(rs$mat, rs$meta, 20L, 5L, 3L)
  mat <- sel$mat; meta <- sel$meta
  db <- dist_block(mat, meta)
  dp <- phylo_dist(tree_file, sel$ids)
  dom <- dominant_source(mat, meta)
  prev <- db$prev
  strata <- interaction(dom, prevalence_strata(prev, 5), drop = TRUE)
  list(unit = unit_label, ids = sel$ids, d = list(J = db$dj, B = db$db),
       dp = dp, strata = strata, n_bins = length(sel$ids),
       prev = prev)
}

make_cross <- function() {
  A <- region_subset(cache, "V4"); B <- region_subset(cache, "V5_V7")
  shared <- rowSums(A$mat > 0) > 0 & rowSums(B$mat > 0) > 0
  matA0 <- A$mat[shared, , drop = FALSE]
  matB0 <- B$mat[shared, , drop = FALSE]
  matAB0 <- cbind(matA0, matB0)
  metaAB <- rbind(A$meta, B$meta)
  bsAB <- breadth_stats(matAB0, metaAB)
  prevA0 <- rowSums(matA0 > 0)
  prevB0 <- rowSums(matB0 > 0)
  pass <- bsAB$prev >= 20 & bsAB$nhost >= 5 & bsAB$nsrc >= 3 &
          prevA0 >= 5 & prevB0 >= 5
  ids <- rownames(A$mat)[shared][pass]
  matA <- matA0[pass, , drop = FALSE]
  matB <- matB0[pass, , drop = FALSE]
  matAB <- matAB0[pass, , drop = FALSE]
  dbA <- dist_block(matA, A$meta)
  dbB <- dist_block(matB, B$meta)
  dp <- phylo_dist(tree_file, ids)
  dom <- dominant_source(matAB, metaAB)
  prev_total <- prevA0[pass] + prevB0[pass]
  names(prev_total) <- ids
  strata <- interaction(dom, prevalence_strata(prev_total, 5), drop = TRUE)
  list(unit = "V4-V5-V7 (shared)", ids = ids,
       d = list(J = paired_dist(dbA$dj, dbB$dj),
                B = paired_dist(dbA$db, dbB$db)),
       dp = dp, strata = strata, n_bins = length(ids),
       prev = prev_total, n_shared = sum(shared),
       n_pass_global = sum(bsAB$prev >= 20 & bsAB$nhost >= 5 & bsAB$nsrc >= 3))
}

units <- list()
if ("V4" %in% UNITS) units[["V4"]] <- make_within("V4", "V4 (within)")
if ("V5_V7" %in% UNITS) units[["V5_V7"]] <- make_within("V5_V7", "V5-V7 (within)")
if ("CROSS" %in% UNITS) units[["CROSS"]] <- make_cross()

res_rows <- list(); draw_rows <- list()
for (nm in names(units)) {
  u <- units[[nm]]
  cat(sprintf("\n===== null %s: n=%d =====\n", u$unit, u$n_bins))
  tab <- table(u$strata)
  n_perm <- sum(tab[tab >= 2])
  cat(sprintf("strata=%d; strata>=2=%d; bins in strata>=2=%d\n",
              length(tab), sum(tab >= 2), n_perm))
  for (mi in seq_along(c(J = "J", B = "B"))) {
    met <- c(J = "J", B = "B")[mi]
    m <- mantel_fast(u$d[[met]], u$dp, nperm = NPERM,
                     strata = u$strata,
                     seed = 20261231L + 100L * match(nm, names(units)) + mi,
                     alternative = "greater")
    p_prev <- existing_p[unit == u$unit & metric == met, p_prevalence]
    res_rows[[length(res_rows) + 1L]] <- data.table(
      unit = u$unit, metric = met, n_bins = u$n_bins,
      strata_n = length(tab), strata_ge2 = sum(tab >= 2),
      n_permutable = n_perm,
      observed_r = m$r, null_mean = m$null_mean, null_sd = m$null_sd,
      null_q05 = as.numeric(stats::quantile(m$null, 0.05, na.rm = TRUE)),
      null_q95 = as.numeric(stats::quantile(m$null, 0.95, na.rm = TRUE)),
      null_min = min(m$null, na.rm = TRUE),
      null_max = max(m$null, na.rm = TRUE),
      p_source_stratified = m$p_greater,
      SES = m$SES,
      p_prevalence_stratified = p_prev,
      n_perm = NPERM
    )
    draw_rows[[length(draw_rows) + 1L]] <- data.table(
      unit = u$unit, metric = met, draw = seq_along(m$null),
      r = m$null)
    cat(sprintf("  [%s] observed=%.4f; null mean=%.4f; P_source=%.4g; P_prev=%.4g\n",
                met, m$r, m$null_mean, m$p_greater, p_prev))
    flush.console()
  }
}

res <- rbindlist(res_rows)
draws <- rbindlist(draw_rows)
fwrite_safe(res, file.path(OUT, "null_model_results.tsv"))
fwrite_safe(draws, file.path(OUT, "null_draws.tsv"))
saveRDS(draws, file.path(OUT, "null_draws.rds"))
cat("\n=== null model results ===\n"); print(res)
cat("\nDONE 20_null_model_paper\n")
