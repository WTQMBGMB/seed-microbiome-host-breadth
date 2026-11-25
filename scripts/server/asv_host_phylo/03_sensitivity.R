#!/usr/bin/env Rscript
# =============================================================================
# 03_sensitivity.R — key sensitivity analyses
#
# 1) Alternative ASV filters: loose 10/3/2 and strict 30/5/3.
# 2) Within-largest-source: remove cross-study structure.
# 3) Host-species-level profiles for the main filter.
# 4) Extended MRM on the main filter: add source-set Jaccard and host-breadth
#    differences to the prevalence/abundance controls.
#
# Outputs:
#   result_table/sensitivity_results.tsv
#   result_table/mrm_extended_results.tsv
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ape); library(vegan); library(ecodist)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/asv_host_phylo_utils.R"))

NPERM_SENS     <- as.integer(Sys.getenv("NPERM_SENS", "199"))
NPERM_SENS_MRM <- as.integer(Sys.getenv("NPERM_SENS_MRM", "499"))
NPERM_EXT      <- as.integer(Sys.getenv("NPERM_EXT", "999"))

regions <- strsplit(Sys.getenv("SENS_REGIONS", "v4,v5_v7,its1"), ",")[[1]]
metric_jb <- c(J = "dj", B = "db")
metric_map <- c(J = "dj", B = "db", H = "dhj")  # H = host-level Jaccard

sens_rows <- list(); ext_rows <- list()

add_mantel_row <- function(region, analysis, subset, metric, d, dp,
                           n_asv, n_units, nperm, seed, strata = NULL) {
  m <- mantel_fast(d, dp, nperm = nperm, strata = strata,
                   seed = seed, alternative = "greater")
  data.table::data.table(
    region = region, analysis = analysis, subset = subset, metric = metric,
    n_asv = n_asv, n_units = n_units, n_pairs = m$n_pairs,
    r = m$r, p = m$p_greater, null_mean = m$null_mean,
    null_sd = m$null_sd, SES = m$SES
  )
}

# ---------------------------------------------------------------------------
# 1. Alternative filters (main tree replaced by the filter-specific tree)
# ---------------------------------------------------------------------------
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  cache <- load_region_cache(reg, verbose = FALSE)
  for (fl in c("loose10_3_2", "strict30_5_3")) {
    cand <- data.table::fread(cand_file(reg, fl), sep = "\t")
    ids <- intersect(cand$ASV_ID, cache$ids)
    mat <- cache$mat[ids, , drop = FALSE]
    tax <- cache$tax[match(ids, cache$ids), ]
    cat(sprintf("[alt-filter] %s / %s: %d ASVs\n", reg, fl, length(ids)))
    prof <- build_profile_core(
      mat = mat, meta = cache$meta, tree_file = tree_file(reg, fl),
      region = reg, filter_label = fl, tax = tax, verbose = FALSE
    )
    strata <- prevalence_strata(prof$prev, 5)
    for (mi in seq_along(metric_jb)) {
      metric <- names(metric_jb)[mi]
      sens_rows[[length(sens_rows) + 1L]] <- add_mantel_row(
        reg, "alt_filter", fl, metric, prof[[metric_jb[mi]]], prof$dp,
        length(prof$ids), ncol(mat), NPERM_SENS,
        seed = 20261201L + 100L * ri + mi, strata = NULL
      )
    }
    rm(prof, mat); invisible(gc())
  }
}

# ---------------------------------------------------------------------------
# 2. Within-largest-source
# ---------------------------------------------------------------------------
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  cache <- load_region_cache(reg, verbose = FALSE)
  src_tab <- sort(table(cache$source), decreasing = TRUE)
  src <- names(src_tab)[1L]
  idx <- which(cache$source == src)
  mat_src <- cache$mat[, idx, drop = FALSE]
  meta_src <- cache$meta[idx]
  cand <- data.table::fread(cand_file(reg, "main20_5_3"), sep = "\t")
  ids <- intersect(cand$ASV_ID, cache$ids)
  prev_src <- rowSums(mat_src[ids, , drop = FALSE] > 0)
  thr <- max(5L, ceiling(0.02 * length(idx)))
  ids2 <- ids[prev_src >= thr]
  if (length(ids2) < 30L) {
    thr <- 3L
    ids2 <- ids[prev_src >= thr]
  }
  cat(sprintf("[within-source] %s / %s: %d samples; %d ASVs (prev>=%d)\n",
              reg, src, length(idx), length(ids2), thr))
  tax <- cache$tax[match(ids2, cache$ids), ]
  prof <- build_profile_core(
    mat = mat_src[ids2, , drop = FALSE], meta = meta_src,
    tree_file = tree_file(reg, "main20_5_3"), region = reg,
    filter_label = paste0("within_", src), tax = tax, verbose = FALSE
  )
  strata <- prevalence_strata(prof$prev, 5)
  for (mi in seq_along(metric_jb)) {
    metric <- names(metric_jb)[mi]
    sens_rows[[length(sens_rows) + 1L]] <- add_mantel_row(
      reg, "within_largest_source", src, metric, prof[[metric_jb[mi]]],
      prof$dp, length(prof$ids), length(idx), NPERM_SENS_MRM,
      seed = 20261211L + 100L * ri + mi, strata = strata
    )
    d <- prof[[metric_jb[mi]]]
    dp <- prof$dp; dpv <- prof$dpv; dab <- prof$dab
    mrm <- ecodist::MRM(d ~ dp + dpv + dab, nperm = NPERM_SENS_MRM)
    co <- mrm$coef
    ext_rows[[length(ext_rows) + 1L]] <- data.table::data.table(
      region = reg, analysis = "within_largest_source", subset = src,
      metric = metric, term = rownames(co),
      coef = as.numeric(co[, 1L]), p = as.numeric(co[, 2L]),
      r_squared = as.numeric(mrm$r.squared[1L]),
      r_squared_p = as.numeric(mrm$r.squared[2L])
    )
  }
  rm(prof, mat_src); invisible(gc())
}

# ---------------------------------------------------------------------------
# 3. Host-species-level main profiles
# ---------------------------------------------------------------------------
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  p <- readRDS(profile_file(reg, "main20_5_3"))
  for (metric in c("H", "HB")) {
    d <- if (metric == "H") p$dhj else p$dhb
    sens_rows[[length(sens_rows) + 1L]] <- add_mantel_row(
      reg, "host_level", "main20_5_3", metric, d, p$dp,
      length(p$ids), length(p$host_ids), NPERM_EXT,
      seed = 20261221L + 100L * ri + (metric == "HB"), strata = NULL
    )
  }
  dp <- p$dp; dpv <- p$dpv; dab <- p$dab
  d <- p$dhj
  mrm <- ecodist::MRM(d ~ dp + dpv + dab, nperm = NPERM_EXT)
  co <- mrm$coef
  ext_rows[[length(ext_rows) + 1L]] <- data.table::data.table(
    region = reg, analysis = "host_level", subset = "main20_5_3",
    metric = "H", term = rownames(co),
    coef = as.numeric(co[, 1L]), p = as.numeric(co[, 2L]),
    r_squared = as.numeric(mrm$r.squared[1L]),
    r_squared_p = as.numeric(mrm$r.squared[2L])
  )
  rm(p); invisible(gc())
}

# ---------------------------------------------------------------------------
# 4. Extended MRM on the main sample-level analysis
# ---------------------------------------------------------------------------
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  p <- readRDS(profile_file(reg, "main20_5_3"))
  for (mi in seq_along(metric_jb)) {
    metric <- names(metric_jb)[mi]
    d <- p[[metric_jb[mi]]]
    dp <- p$dp; dpv <- p$dpv; dab <- p$dab; dsrc <- p$dsrc; dhost <- p$dhost
    cat(sprintf("[extended MRM] %s / %s (%d perms)\n", reg, metric, NPERM_EXT))
    mrm <- ecodist::MRM(d ~ dp + dpv + dab + dsrc + dhost, nperm = NPERM_EXT)
    co <- mrm$coef
    ext_rows[[length(ext_rows) + 1L]] <- data.table::data.table(
      region = reg, analysis = "main_extended_mrm", subset = "main20_5_3",
      metric = metric, term = rownames(co),
      coef = as.numeric(co[, 1L]), p = as.numeric(co[, 2L]),
      r_squared = as.numeric(mrm$r.squared[1L]),
      r_squared_p = as.numeric(mrm$r.squared[2L])
    )
  }
  rm(p); invisible(gc())
}


# ---------------------------------------------------------------------------
# 5. Close-vs-far phylogenetic contrast
#    Statistic = mean(dist | far pairs) - mean(dist | closest 10% pairs).
#    Positive values mean close relatives have more similar distributions.
#    Null: prevalence-stratified ASV-label permutation (999).
# ---------------------------------------------------------------------------
close_rows <- list()
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  p <- readRDS(profile_file(reg, "main20_5_3"))
  dpmat <- as.matrix(p$dp); lt <- lower.tri(dpmat); pd <- dpmat[lt]
  idx <- which(lt, arr.ind = TRUE); i <- idx[, 1L]; j <- idx[, 2L]
  q10 <- as.numeric(stats::quantile(pd, 0.10))
  q50 <- as.numeric(stats::quantile(pd, 0.50))
  close <- pd <= q10; far <- pd >= q50
  strata <- prevalence_strata(p$prev, 5)
  n <- length(p$ids)
  for (metric in c("J", "B")) {
    dmat <- as.matrix(p[[if (metric == "J") "dj" else "db"]])
    qd <- dmat[lt]
    obs <- mean(qd[far]) - mean(qd[close])
    set.seed(20261231L + 100L * ri + (metric == "B"))
    null <- numeric(999L)
    for (b in seq_len(999L)) {
      perm <- integer(n)
      for (s in unique(strata)) {
        ii <- which(strata == s)
        perm[ii] <- sample(ii)
      }
      qp <- dmat[cbind(perm[i], perm[j])]
      null[b] <- mean(qp[far]) - mean(qp[close])
    }
    mu <- mean(null); sdv <- stats::sd(null)
    close_rows[[length(close_rows) + 1L]] <- data.table::data.table(
      region = reg, metric = metric,
      n_close = sum(close), n_far = sum(far),
      phylo_q10 = q10, phylo_q50 = q50,
      mean_close = mean(qd[close]), mean_far = mean(qd[far]),
      contrast = obs, null_mean = mu, null_sd = sdv,
      SES = (obs - mu) / sdv,
      p = (1 + sum(null >= obs)) / (length(null) + 1)
    )
    cat(sprintf("[close-vs-far] %s / %s: contrast=%.4f SES=%.2f p=%.4f\n",
                reg, metric, obs, (obs - mu) / sdv,
                (1 + sum(null >= obs)) / (length(null) + 1)))
  }
  rm(p); invisible(gc())
}
close_res <- data.table::rbindlist(close_rows)
fwrite_safe(close_res, file.path(OUT_ROOT, "result_table", "close_pair_contrast.tsv"))

# Save the key sensitivity tables before the optional Spearman/tree-choice
# sections, so a slow optional step cannot hide them.
sens <- data.table::rbindlist(sens_rows, fill = TRUE)
sens[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = .(analysis)]
ext <- data.table::rbindlist(ext_rows, fill = TRUE)
ext[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = .(analysis, metric)]
fwrite_safe(sens, file.path(OUT_ROOT, "result_table", "sensitivity_results.tsv"))
fwrite_safe(ext,  file.path(OUT_ROOT, "result_table", "mrm_extended_results.tsv"))


# ---------------------------------------------------------------------------
# 6. Spearman Mantel sensitivity (outlier-robust)
# ---------------------------------------------------------------------------
spear_rows <- list()
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  p <- readRDS(profile_file(reg, "main20_5_3"))
  for (metric in c("J", "B")) {
    d <- p[[if (metric == "J") "dj" else "db"]]
    dp <- p$dp
    # ecodist mrank=TRUE = Spearman; C implementation is much faster than
    # ranking every permuted vector in R.
    e <- ecodist::mantel(d ~ dp, nperm = 999L, mrank = TRUE)
    spear_rows[[length(spear_rows) + 1L]] <- data.table::data.table(
      region = reg, metric = metric, n_asv = length(p$ids),
      n_pairs = length(d), rho = as.numeric(e["mantelr"]),
      p = as.numeric(e["pval1"]), null_mean = NA_real_,
      null_sd = NA_real_, SES = NA_real_
    )
    cat(sprintf("[spearman] %s / %s: rho=%.4f p=%.4f
",
                reg, metric, as.numeric(e["mantelr"]), as.numeric(e["pval1"])))
  }
  rm(p); invisible(gc())
}
spear <- data.table::rbindlist(spear_rows)
fwrite_safe(spear, file.path(OUT_ROOT, "result_table", "mantel_spearman.tsv"))


# ---------------------------------------------------------------------------
# 7. Tree-choice sensitivity: prune the loose candidate tree to the main ASVs
#    (larger phylogenetic context than the main-only tree).
# ---------------------------------------------------------------------------
tree_rows <- list()
for (ri in seq_along(regions)) {
  reg <- regions[ri]
  p <- readRDS(profile_file(reg, "main20_5_3"))
  tr <- ape::read.tree(tree_file(reg, "loose10_3_2"))
  if (!all(p$ids %in% tr$tip.label))
    stop(sprintf("[%s] loose tree lacks main candidates", reg))
  tr <- ape::keep.tip(tr, p$ids)
  dpm <- stats::cophenetic(tr)[p$ids, p$ids, drop = FALSE]
  dpl <- stats::as.dist(dpm)
  strata <- prevalence_strata(p$prev, 5)
  for (metric in c("J", "B")) {
    d <- p[[if (metric == "J") "dj" else "db"]]
    m <- mantel_fast(d, dpl, nperm = 299L, strata = strata,
                     seed = 20261251L + 100L * ri + (metric == "B"),
                     alternative = "greater")
    tree_rows[[length(tree_rows) + 1L]] <- data.table::data.table(
      region = reg, analysis = "loose_tree_pruned", metric = metric,
      n_asv = length(p$ids), n_pairs = m$n_pairs,
      r = m$r, p = m$p_greater, null_mean = m$null_mean,
      null_sd = m$null_sd, SES = m$SES
    )
    cat(sprintf("[tree-choice] %s / %s: r=%.4f SES=%.2f p=%.4f\n",
                reg, metric, m$r, m$SES, m$p_greater))
  }
  rm(p); invisible(gc())
}
tree_res <- data.table::rbindlist(tree_rows)
fwrite_safe(tree_res, file.path(OUT_ROOT, "result_table", "tree_choice_sensitivity.tsv"))

sens <- data.table::rbindlist(sens_rows, fill = TRUE)
sens[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = .(analysis)]
ext <- data.table::rbindlist(ext_rows, fill = TRUE)
ext[, p_adj_BH := stats::p.adjust(p, method = "BH"), by = .(analysis, metric)]

fwrite_safe(sens, file.path(OUT_ROOT, "result_table", "sensitivity_results.tsv"))
fwrite_safe(ext,  file.path(OUT_ROOT, "result_table", "mrm_extended_results.tsv"))

cat("\n=== sensitivity Mantel ===\n"); print(sens)
cat("\n=== sensitivity MRM phylo terms ===\n")
print(ext[term == "dp"])
cat("\nDONE 03_sensitivity\n")
