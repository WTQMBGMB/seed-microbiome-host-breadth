#!/usr/bin/env Rscript
# =============================================================================
# 11_track1_16s.R — Track 1: within-region mapped-bin replication + meta
#
# For each 16S region:
#   - mapped reference bins (97% identity) with within-region prevalence >=20,
#     host species >=5 (adaptive: min(5, region species)), sources >=3
#   - Jaccard / Bray distribution distances vs SILVA reference-tree patristic
#   - Mantel (standard + prevalence-stratified), MRM, correlogram, close-vs-far
# Meta-combine Mantel r across regions (Fisher z, random effects).
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ape); library(vegan); library(ecodist)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/mapped_utils.R"))

NPERM_MAIN  <- as.integer(Sys.getenv("NPERM_MAIN", "9999"))
NPERM_STRAT <- as.integer(Sys.getenv("NPERM_STRAT", "999"))
NPERM_MRM   <- as.integer(Sys.getenv("NPERM_MRM", "1999"))
NPERM_CG    <- as.integer(Sys.getenv("NPERM_CG", "499"))
NPERM_CF    <- as.integer(Sys.getenv("NPERM_CF", "999"))

DS <- Sys.getenv("DS", "16S")
default_regions <- if (DS == "16S") "V3_V4,V4,V5_V6,V5_V7" else "ITS1,ITS2"
REGIONS <- strsplit(Sys.getenv("REGIONS", default_regions), ",")[[1]]
cfg <- MAPPED_CFG[[DS]]
id_map <- if (DS == "ITS")
  file.path(OUT_ROOT, "Temp_data/trees/its_ref_dup_map.tsv") else NULL
cache <- load_mapped_cache(DS)
TAB_DIR <- file.path(OUT_ROOT, "result_table", "cross_region")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

mantel_rows <- list(); mrm_rows <- list(); cg_rows <- list()
cf_rows <- list(); cand_rows <- list()

for (ri in seq_along(REGIONS)) {
  reg <- REGIONS[ri]
  cat(sprintf("\n===== %s / %s =====\n", DS, reg)); flush.console()
  rs <- region_subset(cache, reg)
  sel <- select_candidates(rs$mat, rs$meta, prev_min = 20L,
                           host_min = 5L, src_min = 3L)
  cat(sprintf("samples=%d; candidates=%d (host_thr=%d, src_thr=%d)\n",
              ncol(sel$mat), length(sel$ids), sel$host_thr, sel$src_thr))
  flush.console()
  if (length(sel$ids) < 10L) {
    warning(sprintf("[%s] too few candidates; skipped", reg)); next
  }

  cand_rows[[reg]] <- cbind(
    data.table(dataset = DS, region = reg,
               host_thr = sel$host_thr, src_thr = sel$src_thr),
    sel$stats[keep == TRUE]
  )

  mat <- sel$mat
  db <- dist_block(mat, sel$meta)
  dp <- phylo_dist(file.path(BASE_DIR, cfg$tree), sel$ids, id_map = id_map)
  strata <- prevalence_strata(db$prev, 5L)

  for (mi in seq_along(c(J = "dj", B = "db"))) {
    metric <- names(c(J = "dj", B = "db"))[mi]
    d <- db[[c(J = "dj", B = "db")[mi]]]
    cat(sprintf("  [%s] Mantel ...\n", metric)); flush.console()
    m_std <- mantel_ecodist(d, dp, nperm = NPERM_MAIN)
    m_ses <- mantel_fast(d, dp, nperm = 999L,
                         seed = 20261141L + 100L * ri + mi,
                         alternative = "greater")
    m_str <- mantel_fast(d, dp, nperm = NPERM_STRAT, strata = strata,
                         seed = 20261151L + 100L * ri + mi,
                         alternative = "greater")
    mantel_rows[[length(mantel_rows) + 1L]] <- data.table(
      dataset = DS, region = reg, metric = metric,
      n_bins = length(sel$ids), n_pairs = m_ses$n_pairs,
      r = m_std$r, p_standard = m_std$p, p_two_sided = m_std$p_two,
      CI_low = m_std$CI_low, CI_high = m_std$CI_high,
      null_mean = m_ses$null_mean, null_sd = m_ses$null_sd,
      SES = m_ses$SES,
      p_stratified = m_str$p_greater,
      SES_stratified = m_str$SES
    )

    cat(sprintf("  [%s] MRM ...\n", metric)); flush.console()
    mrm1 <- mrm_fit(d, dp, db$dpv, db$dab, nperm = NPERM_MRM)
    mrm1[, `:=`(dataset = DS, region = reg, metric = metric,
                model = "primary",
                extra_terms = attr(mrm1, "extra_terms"))]
    mrm2 <- mrm_fit(d, dp, db$dpv, db$dab,
                    extra = list(dhost = db$dhost, dsrc = db$dsrc),
                    nperm = NPERM_MRM)
    mrm2[, `:=`(dataset = DS, region = reg, metric = metric,
                model = "extended",
                extra_terms = attr(mrm2, "extra_terms"))]
    mrm_rows[[length(mrm_rows) + 1L]] <- rbind(mrm1, mrm2)

    cf <- close_far_fast(d, dp, strata = strata, nperm = NPERM_CF,
                         seed = 20261161L + 100L * ri + mi)
    if (!is.null(cf)) {
      cf_rows[[length(cf_rows) + 1L]] <- data.table(
        dataset = DS, region = reg, metric = metric,
        n_bins = length(sel$ids),
        contrast = cf$contrast, p = cf$p,
        null_mean = cf$null_mean, null_sd = cf$null_sd, SES = cf$SES,
        n_close = cf$n_close, n_far = cf$n_far
      )
    }
  }

  cat(sprintf("  correlogram (%d perms) ...\n", NPERM_CG)); flush.console()
  cg <- correlogram_fast(list(J = db$dj, B = db$db), dp, nbin = 10L,
                         nperm = NPERM_CG, strata = strata,
                         seed = 20261171L + ri)
  for (metric in names(cg)) {
    cg_rows[[length(cg_rows) + 1L]] <-
      cbind(data.table(dataset = DS, region = reg,
                       n_bins = length(sel$ids)), cg[[metric]])
  }
  rm(rs, sel, mat, db, dp); invisible(gc())
}

mantel_res <- rbindlist(mantel_rows, fill = TRUE)
mrm_res <- rbindlist(mrm_rows, fill = TRUE)
cg_res <- rbindlist(cg_rows, fill = TRUE)
cf_res <- rbindlist(cf_rows, fill = TRUE)
cand_res <- rbindlist(cand_rows, fill = TRUE)

# meta-analysis: core 3 regions vs all 4 (V5_V6 is exploratory)
region_sets <- if (DS == "16S") {
  list(core3 = c("V3_V4", "V4", "V5_V7"),
       all4 = c("V3_V4", "V4", "V5_V6", "V5_V7"))
} else {
  list(all = REGIONS)
}
meta_rows <- list()
for (met in c("J", "B")) {
  for (set_name in names(region_sets)) {
    regs <- region_sets[[set_name]]
    sub <- mantel_res[metric == met & region %in% regs,
                      .(region, metric, r, n_bins)]
    mm <- meta_mantel(sub)
    if (!is.null(mm)) {
      mm[, `:=`(dataset = DS, metric = met, region_set = set_name)]
      meta_rows[[length(meta_rows) + 1L]] <- mm
    }
  }
}
meta_res <- rbindlist(meta_rows, fill = TRUE)

tag <- tolower(DS)
fwrite_safe(cand_res, file.path(TAB_DIR, sprintf("track1_%s_candidates.tsv", tag)))
fwrite_safe(mantel_res, file.path(TAB_DIR, sprintf("track1_%s_mantel.tsv", tag)))
fwrite_safe(mrm_res, file.path(TAB_DIR, sprintf("track1_%s_mrm.tsv", tag)))
fwrite_safe(cg_res, file.path(TAB_DIR, sprintf("track1_%s_correlogram.tsv", tag)))
fwrite_safe(cf_res, file.path(TAB_DIR, sprintf("track1_%s_close_far.tsv", tag)))
fwrite_safe(meta_res, file.path(TAB_DIR, sprintf("track1_%s_meta.tsv", tag)))

cat("\n=== Mantel ===\n"); print(mantel_res)
cat("\n=== MRM phylo term ===\n"); print(mrm_res[term == "dp"])
cat("\n=== close-vs-far ===\n"); print(cf_res)
cat("\n=== meta ===\n"); print(meta_res)
cat("\nDONE 11_track1_16s\n")
