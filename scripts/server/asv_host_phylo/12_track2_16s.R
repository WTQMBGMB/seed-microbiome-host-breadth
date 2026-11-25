#!/usr/bin/env Rscript
# =============================================================================
# 12_track2_16s.R — Track 2: cross-region shared-bin analysis (16S)
#
# For each pair of 16S regions:
#   - shared mapped bins detected in both regions
#   - primary filter: total prevalence >=20, total hosts >=5, total sources >=3,
#     and >=5 positives in EACH region
#   - within-region Jaccard / Bray distances, then equal-weight mean of the two
#     region-specific distances (avoids pooled structural zeros)
#   - Mantel / MRM vs SILVA reference-tree patristic distance
#   - sensitivity: within-region prevalence >=10; min / max instead of mean
#
# Pair-specific results are reported separately; pairs are NOT independent
# (they share bins and regions), so no formal meta-analysis across pairs.
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
NPERM_CF    <- as.integer(Sys.getenv("NPERM_CF", "999"))
NPERM_SENS  <- as.integer(Sys.getenv("NPERM_SENS", "999"))
WITHIN_MIN  <- as.integer(Sys.getenv("WITHIN_MIN", "5"))

DS <- Sys.getenv("DS", "16S")
default_regions <- if (DS == "16S") "V3_V4,V4,V5_V6,V5_V7" else "ITS1,ITS2"
REGIONS <- strsplit(Sys.getenv("REGIONS", default_regions), ",")[[1]]
cfg <- MAPPED_CFG[[DS]]
id_map <- if (DS == "ITS")
  file.path(OUT_ROOT, "Temp_data/trees/its_ref_dup_map.tsv") else NULL
cache <- load_mapped_cache(DS)
TAB_DIR <- file.path(OUT_ROOT, "result_table", "cross_region")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
rs <- setNames(lapply(REGIONS, function(r) region_subset(cache, r)), REGIONS)
prev_all <- lapply(rs, function(x) rowSums(x$mat > 0))
PAIRS <- as.data.table(t(combn(REGIONS, 2)))
setnames(PAIRS, c("region_a", "region_b"))

pair_meta <- list(); mantel_rows <- list(); mrm_rows <- list()
cf_rows <- list(); sens_rows <- list()

for (pi in seq_len(nrow(PAIRS))) {
  ra <- PAIRS$region_a[pi]; rb <- PAIRS$region_b[pi]
  overlap <- (ra == "V3_V4" && rb == "V4") ||
             (ra == "V5_V6" && rb == "V5_V7")
  pair_type <- if (overlap) "overlapping_amplicon" else "non_overlapping"
  cat(sprintf("\n===== %s / %s-%s (%s) =====\n", DS, ra, rb, pair_type))
  flush.console()

  A <- rs[[ra]]; B <- rs[[rb]]
  shared <- prev_all[[ra]] > 0 & prev_all[[rb]] > 0
  ids0 <- rownames(A$mat)[shared]
  matA0 <- A$mat[shared, , drop = FALSE]
  matB0 <- B$mat[shared, , drop = FALSE]
  matAB0 <- cbind(matA0, matB0)
  metaAB <- rbind(A$meta, B$meta)
  bsAB <- breadth_stats(matAB0, metaAB)
  prevA0 <- rowSums(matA0 > 0)
  prevB0 <- rowSums(matB0 > 0)

  pass_global <- bsAB$prev >= 20 & bsAB$nhost >= 5 & bsAB$nsrc >= 3
  pass_pair <- pass_global & prevA0 >= WITHIN_MIN & prevB0 >= WITHIN_MIN
  cat(sprintf("shared=%d; pass global=%d; pass pair=%d\n",
              length(ids0), sum(pass_global), sum(pass_pair)))
  flush.console()

  pair_meta[[pi]] <- data.table(
    dataset = DS, region_a = ra, region_b = rb, pair_type = pair_type,
    n_shared = length(ids0), n_pass_global = sum(pass_global),
    n_pass_pair = sum(pass_pair), within_min = WITHIN_MIN
  )
  if (sum(pass_pair) < 20L) {
    warning(sprintf("[%s-%s] too few shared bins; skipped", ra, rb)); next
  }

  ids <- ids0[pass_pair]
  matA <- matA0[pass_pair, , drop = FALSE]
  matB <- matB0[pass_pair, , drop = FALSE]
  matAB <- matAB0[pass_pair, , drop = FALSE]
  metaA <- A$meta; metaB <- B$meta
  bsA <- breadth_stats(matA, metaA)
  bsB <- breadth_stats(matB, metaB)
  bsAB <- breadth_stats(matAB, metaAB)
  prevA <- prevA0[pass_pair]; prevB <- prevB0[pass_pair]
  names(prevA) <- names(prevB) <- ids

  dA <- dist_block(matA, metaA)
  dB <- dist_block(matB, metaB)
  dAB <- dist_block(matAB, metaAB)

  pair_dist <- function(metric, fun = mean) {
    ma <- as.matrix(dA[[metric]])
    mb <- as.matrix(dB[[metric]])
    mm <- fun(ma, mb)
    rownames(mm) <- colnames(mm) <- ids
    stats::as.dist(mm)
  }
  dj <- pair_dist("dj", function(x, y) (x + y) / 2)
  db <- pair_dist("db", function(x, y) (x + y) / 2)
  if (!identical(attr(dj, "Labels"), ids) || !identical(attr(db, "Labels"), ids))
    stop("pair distance labels misaligned")

  dp <- phylo_dist(file.path(BASE_DIR, cfg$tree), ids, id_map = id_map)
  prev_total <- prevA + prevB
  names(prev_total) <- ids
  mean_rel_bal <- (bsA$mean_rel + bsB$mean_rel) / 2
  names(mean_rel_bal) <- ids
  imbalance <- abs(log10(prevA) - log10(prevB))
  dpv <- stats::as.dist(abs(outer(log10(prev_total), log10(prev_total), "-")))
  dab <- stats::as.dist(abs(outer(log10(mean_rel_bal + 1e-6),
                                  log10(mean_rel_bal + 1e-6), "-")))
  dprev <- stats::as.dist(abs(outer(imbalance, imbalance, "-")))
  dhost <- stats::as.dist(abs(outer(bsAB$nhost, bsAB$nhost, "-")))
  dsrc <- dAB$dsrc
  strata <- prevalence_strata(prev_total, 5L)

  for (mi in seq_along(c(J = "dj", B = "db"))) {
    metric <- names(c(J = "dj", B = "db"))[mi]
    d <- list(J = dj, B = db)[[metric]]
    m_std <- mantel_ecodist(d, dp, nperm = NPERM_MAIN)
    m_ses <- mantel_fast(d, dp, nperm = 999L,
                         seed = 20261201L + 100L * pi + mi,
                         alternative = "greater")
    m_str <- mantel_fast(d, dp, nperm = NPERM_STRAT, strata = strata,
                         seed = 20261211L + 100L * pi + mi,
                         alternative = "greater")
    mantel_rows[[length(mantel_rows) + 1L]] <- data.table(
      dataset = DS, region_a = ra, region_b = rb, pair_type = pair_type,
      metric = metric, n_bins = length(ids), n_pairs = m_ses$n_pairs,
      r = m_std$r, p_standard = m_std$p, p_two_sided = m_std$p_two,
      CI_low = m_std$CI_low, CI_high = m_std$CI_high,
      null_mean = m_ses$null_mean, null_sd = m_ses$null_sd, SES = m_ses$SES,
      p_stratified = m_str$p_greater, SES_stratified = m_str$SES
    )

    mrm1 <- mrm_fit(d, dp, dpv, dab, extra = list(dprev = dprev),
                    nperm = NPERM_MRM)
    mrm1[, `:=`(dataset = DS, region_a = ra, region_b = rb,
                pair_type = pair_type, metric = metric, model = "primary",
                extra_terms = attr(mrm1, "extra_terms"))]
    mrm2 <- mrm_fit(d, dp, dpv, dab,
                    extra = list(dprev = dprev, dhost = dhost, dsrc = dsrc),
                    nperm = NPERM_MRM)
    mrm2[, `:=`(dataset = DS, region_a = ra, region_b = rb,
                pair_type = pair_type, metric = metric, model = "extended",
                extra_terms = attr(mrm2, "extra_terms"))]
    mrm_rows[[length(mrm_rows) + 1L]] <- rbind(mrm1, mrm2)

    cf <- close_far_fast(d, dp, strata = strata, nperm = NPERM_CF,
                         seed = 20261221L + 100L * pi + mi)
    if (!is.null(cf)) {
      cf_rows[[length(cf_rows) + 1L]] <- data.table(
        dataset = DS, region_a = ra, region_b = rb, pair_type = pair_type,
        metric = metric, n_bins = length(ids),
        contrast = cf$contrast, p = cf$p, null_mean = cf$null_mean,
        null_sd = cf$null_sd, SES = cf$SES,
        n_close = cf$n_close, n_far = cf$n_far
      )
    }
  }

  # sensitivity: within-region prevalence >=10
  pass10 <- pass_global & prevA0 >= 10L & prevB0 >= 10L
  if (sum(pass10) >= 20L) {
    ids10 <- ids0[pass10]
    dA10 <- dist_block(matA0[pass10, , drop = FALSE], metaA)
    dB10 <- dist_block(matB0[pass10, , drop = FALSE], metaB)
    dp10 <- phylo_dist(file.path(BASE_DIR, cfg$tree), ids10, id_map = id_map)
    strata10 <- prevalence_strata(prevA0[pass10] + prevB0[pass10], 5L)
    for (metric in c("J", "B")) {
      ma <- as.matrix(dA10[[c(J = "dj", B = "db")[metric]]])
      mb <- as.matrix(dB10[[c(J = "dj", B = "db")[metric]]])
      mm <- (ma + mb) / 2
      rownames(mm) <- colnames(mm) <- ids10
      d10 <- stats::as.dist(mm)
      ms <- mantel_ecodist(d10, dp10, nperm = NPERM_SENS)
      sens_rows[[length(sens_rows) + 1L]] <- data.table(
        dataset = DS, region_a = ra, region_b = rb, pair_type = pair_type,
        metric = metric, analysis = "within_min10", n_bins = length(ids10),
        r = ms$r, p = ms$p
      )
    }
  }

  # sensitivity: min / max instead of mean
  for (metric in c("J", "B")) {
    ma <- as.matrix(dA[[c(J = "dj", B = "db")[metric]]])
    mb <- as.matrix(dB[[c(J = "dj", B = "db")[metric]]])
    for (agg in c("min", "max")) {
      mm <- if (agg == "min") pmin(ma, mb) else pmax(ma, mb)
      rownames(mm) <- colnames(mm) <- ids
      dd <- stats::as.dist(mm)
      ms <- mantel_ecodist(dd, dp, nperm = NPERM_SENS)
      sens_rows[[length(sens_rows) + 1L]] <- data.table(
        dataset = DS, region_a = ra, region_b = rb, pair_type = pair_type,
        metric = metric, analysis = paste0("pair_", agg), n_bins = length(ids),
        r = ms$r, p = ms$p
      )
    }
  }

  rm(A, B, matA0, matB0, matAB0, bsAB, dA, dB, dAB, dj, db, dp)
  invisible(gc())
}

pair_meta_res <- rbindlist(pair_meta, fill = TRUE)
mantel_res <- rbindlist(mantel_rows, fill = TRUE)
mrm_res <- rbindlist(mrm_rows, fill = TRUE)
cf_res <- rbindlist(cf_rows, fill = TRUE)
sens_res <- rbindlist(sens_rows, fill = TRUE)

summary_res <- mantel_res[, .(
  n_pairs = .N,
  n_positive = sum(r > 0),
  median_r = stats::median(r),
  min_r = min(r),
  max_r = max(r),
  n_p_standard_005 = sum(p_standard < 0.05),
  n_p_stratified_005 = sum(p_stratified < 0.05)
), by = .(dataset, metric, pair_type)]

tag <- tolower(DS)
fwrite_safe(pair_meta_res, file.path(TAB_DIR, sprintf("track2_%s_pairs.tsv", tag)))
fwrite_safe(mantel_res, file.path(TAB_DIR, sprintf("track2_%s_mantel.tsv", tag)))
fwrite_safe(mrm_res, file.path(TAB_DIR, sprintf("track2_%s_mrm.tsv", tag)))
fwrite_safe(cf_res, file.path(TAB_DIR, sprintf("track2_%s_close_far.tsv", tag)))
fwrite_safe(sens_res, file.path(TAB_DIR, sprintf("track2_%s_sensitivity.tsv", tag)))
fwrite_safe(summary_res, file.path(TAB_DIR, sprintf("track2_%s_summary.tsv", tag)))

cat("\n=== pair counts ===\n"); print(pair_meta_res)
cat("\n=== Mantel ===\n"); print(mantel_res)
cat("\n=== summary ===\n"); print(summary_res)
cat("\nDONE 12_track2_16s\n")
