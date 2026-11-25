#!/usr/bin/env Rscript
# =============================================================================
# asv_host_phylo_utils.R — shared helpers for the ASV-host-phylogeny analysis
#
# Design:
#   - Cache per region the union of candidate ASVs (loose filter) and aligned
#     metadata, so repeated filter/sensitivity runs do not reread 150 MB tables.
#   - build_profile_core() turns an ASV x sample matrix into:
#       * ASV x ASV Jaccard (presence) and Bray (sample-wise relative abundance)
#       * ASV patristic distances from the candidate tree
#       * covariate distances: prevalence, mean abundance, host breadth,
#         source-set Jaccard
#       * host-species-level profiles (Jaccard / Bray)
#   - mantel_fast() computes one-sided permutation Mantel tests, optionally
#     permuting ASV labels within prevalence strata.
#   - correlogram_fast() computes phylo-distance-bin means and
#     prevalence-stratified null envelopes.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(vegan)
})

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "<SERVER_BASE>")
OUT_ROOT <- file.path(BASE_DIR, "result/ASV_host_phylo")

REGION_CFG <- list(
  v4 = list(
    css  = "standard_region/v4/css_v4_asvtab.txt",
    meta = "standard_merge/metadata_16S.txt", rf = "V4",
    tax  = "filtered_region/V4/V4_merged_taxonomy_filtered.txt"),
  v5_v7 = list(
    css  = "standard_region/v5_v7/css_v5_v7_asvtab.txt",
    meta = "standard_merge/metadata_16S.txt", rf = "V5_V7",
    tax  = "filtered_region/V5_V7/V5_V7_merged_taxonomy_filtered.txt"),
  its1 = list(
    css  = "standard_region/its1/css_its1_asvtab.txt",
    meta = "standard_merge/metadata_ITS.txt", rf = "ITS1",
    tax  = "filtered_region/ITS1/ITS1_merged_taxonomy_filtered.txt")
)

tree_file <- function(region, filter) {
  file.path(OUT_ROOT, "Temp_data/trees", sprintf("%s_%s.nwk", region, filter))
}
cand_file <- function(region, filter) {
  file.path(OUT_ROOT, "Temp_data/candidates",
            sprintf("%s_%s_candidates.tsv", region, filter))
}
profile_file <- function(region, filter) {
  file.path(OUT_ROOT, "Temp_data/rds", sprintf("profiles_%s_%s.rds", region, filter))
}

# ---------------------------------------------------------------------------
# Region cache: loose-filter ASV rows + aligned metadata + taxonomy
# ---------------------------------------------------------------------------
load_region_cache <- function(region, force = FALSE, verbose = TRUE) {
  stopifnot(region %in% names(REGION_CFG))
  cache_file <- file.path(OUT_ROOT, "Temp_data/rds", sprintf("cache_%s.rds", region))
  if (file.exists(cache_file) && !force) {
    if (verbose) cat(sprintf("[cache] load %s\n", basename(cache_file)))
    return(readRDS(cache_file))
  }
  g <- REGION_CFG[[region]]
  if (verbose) cat(sprintf("[cache] build %s from CSS table ...\n", region))
  raw <- data.table::fread(file.path(BASE_DIR, g$css), sep = "\t", header = TRUE)
  ids <- as.character(raw[[1L]])
  mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(mat) <- ids
  rm(raw); invisible(gc())

  meta <- data.table::fread(file.path(BASE_DIR, g$meta), sep = "\t")
  meta <- meta[Region == g$rf]
  common <- intersect(colnames(mat), as.character(meta$SampleID))
  if (!length(common)) {
    colnames(mat) <- sub("^X", "", colnames(mat))
    common <- intersect(colnames(mat), as.character(meta$SampleID))
  }
  if (length(common) < 10L) stop(sprintf("[%s] too few aligned samples", region))
  mat <- mat[, common, drop = FALSE]
  meta <- meta[match(common, meta$SampleID)]
  rownames(meta) <- meta$SampleID

  # Keep the union of all three candidate filters (loose is the largest).
  loose <- data.table::fread(cand_file(region, "loose10_3_2"), sep = "\t")$ASV_ID
  keep <- intersect(loose, rownames(mat))
  mat <- mat[keep, , drop = FALSE]

  tax <- data.table::fread(file.path(BASE_DIR, g$tax), sep = "\t")
  tax <- tax[match(rownames(mat), tax$taxid)]
  tax_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                        names(tax))
  tax <- tax[, tax_cols, with = FALSE]

  cache <- list(
    region = region,
    mat = mat,
    meta = meta,
    tax = tax,
    ids = rownames(mat),
    sample_ids = colnames(mat),
    species = as.character(meta$species),
    source = as.character(meta$source)
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, cache_file)
  if (verbose) cat(sprintf("[cache] saved %s: %d ASVs x %d samples\n",
                           basename(cache_file), nrow(mat), ncol(mat)))
  cache
}

# ---------------------------------------------------------------------------
# Build all ASV-level and host-level distance objects from one matrix
# ---------------------------------------------------------------------------
build_profile_core <- function(mat, meta, tree_file, region = NA_character_,
                               filter_label = NA_character_,
                               tax = NULL, min_host_samples = 3L,
                               verbose = TRUE) {
  stopifnot(is.matrix(mat), nrow(mat) >= 10L, ncol(mat) >= 10L)
  ids <- rownames(mat)
  if (anyDuplicated(ids)) stop("duplicated ASV IDs")
  if (verbose) cat(sprintf("[profile] %s/%s: %d ASVs x %d samples\n",
                           region, filter_label, nrow(mat), ncol(mat)))

  prev <- rowSums(mat > 0)
  col_tot <- colSums(mat)
  col_tot[!is.finite(col_tot) | col_tot <= 0] <- 1
  rel <- sweep(mat, 2L, col_tot, "/")
  mean_rel <- rowMeans(rel)
  sp <- as.character(meta$species)
  so <- as.character(meta$source)
  nhost <- vapply(seq_len(nrow(mat)), function(i)
    data.table::uniqueN(sp[mat[i, ] > 0]), integer(1))
  nsrc <- vapply(seq_len(nrow(mat)), function(i)
    data.table::uniqueN(so[mat[i, ] > 0]), integer(1))
  names(prev) <- names(mean_rel) <- names(nhost) <- names(nsrc) <- ids

  # ---- ASV x ASV distribution distances ----
  dj <- vegan::vegdist(mat > 0, method = "jaccard", binary = TRUE)
  if (requireNamespace("parallelDist", quietly = TRUE)) {
    db <- parallelDist::parDist(rel, method = "bray",
                                threads = as.integer(Sys.getenv("DIST_THREADS", "4")))
    if (is.null(attr(db, "Labels"))) attr(db, "Labels") <- ids
  } else {
    db <- vegan::vegdist(rel, method = "bray")
  }
  if (!identical(attr(dj, "Labels"), ids) || !identical(attr(db, "Labels"), ids))
    stop("distance labels are not aligned to ASV IDs")

  # ---- microbial tree ----
  tr <- ape::read.tree(tree_file)
  if (!all(ids %in% tr$tip.label))
    stop(sprintf("tree %s lacks %d candidate tips", basename(tree_file),
                 sum(!ids %in% tr$tip.label)))
  if (length(tr$tip.label) != length(ids)) tr <- ape::keep.tip(tr, ids)
  dpm <- stats::cophenetic(tr)[ids, ids, drop = FALSE]
  dp <- stats::as.dist(dpm)
  if (!identical(attr(dp, "Labels"), ids)) stop("phylo distance labels misaligned")

  # ---- covariate distances ----
  dpv    <- stats::as.dist(abs(outer(log10(prev), log10(prev), "-")))
  dab    <- stats::as.dist(abs(outer(log10(mean_rel + 1e-6),
                                     log10(mean_rel + 1e-6), "-")))
  dhost  <- stats::as.dist(abs(outer(nhost, nhost, "-")))
  if (length(unique(so)) >= 2L) {
    src_fac <- factor(so, levels = sort(unique(so)))
    src_model <- stats::model.matrix(~ 0 + src_fac)
    src_pa <- (mat > 0) %*% src_model
    dsrc   <- vegan::vegdist(src_pa > 0, method = "jaccard", binary = TRUE)
  } else {
    # Within a single source, the source-set distance is identically zero.
    dsrc <- stats::as.dist(matrix(0, nrow(mat), nrow(mat),
                                  dimnames = list(ids, ids)))
  }
  for (dd in list(dpv, dab, dhost, dsrc))
    if (!identical(attr(dd, "Labels"), ids)) stop("covariate distance labels misaligned")

  # ---- host-species-level profiles ----
  host_tab <- table(sp)
  host_tab <- host_tab[!is.na(names(host_tab)) & nzchar(names(host_tab))]
  min_host <- min_host_samples
  repeat {
    host_keep <- names(host_tab)[host_tab >= min_host]
    if (length(host_keep) >= 10L || min_host <= 1L) break
    min_host <- min_host - 1L
  }
  host_idx <- split(seq_len(ncol(mat)), factor(sp, levels = host_keep))
  host_idx <- host_idx[lengths(host_idx) > 0L]
  host_ids <- names(host_idx)
  host_n <- lengths(host_idx)
  host_rel <- do.call(cbind, lapply(host_idx, function(ii)
    rowMeans(rel[, ii, drop = FALSE])))
  host_npos <- do.call(cbind, lapply(host_idx, function(ii)
    rowSums(mat[, ii, drop = FALSE] > 0)))
  host_rel <- matrix(host_rel, nrow = nrow(mat), dimnames = list(ids, host_ids))
  host_npos <- matrix(host_npos, nrow = nrow(mat), dimnames = list(ids, host_ids))
  det_thr <- ifelse(host_n >= 5L, 2L, 1L)
  host_pa <- sweep(host_npos, 2L, det_thr, ">=")
  dhj <- vegan::vegdist(host_pa, method = "jaccard", binary = TRUE)
  dhb <- vegan::vegdist(host_rel, method = "bray")
  if (verbose) cat(sprintf("[profile] host level: %d host species (min samples %d)\n",
                           length(host_ids), min_host))

  list(
    region = region, filter = filter_label,
    ids = ids, sample_ids = colnames(mat),
    host_ids = host_ids, host_n_samples = host_n,
    prev = prev, mean_rel = mean_rel, nhost = nhost, nsrc = nsrc,
    dj = dj, db = db, dp = dp,
    dpv = dpv, dab = dab, dhost = dhost, dsrc = dsrc,
    dhj = dhj, dhb = dhb,
    tax = tax
  )
}

# ---------------------------------------------------------------------------
# Fast one-sided Mantel test with optional within-stratum label permutation
# ---------------------------------------------------------------------------
mantel_fast <- function(dx, dy, nperm = 999L, strata = NULL, seed = 20261111L,
                        alternative = c("greater", "less", "two.sided"),
                        method = c("pearson", "spearman")) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  if (!identical(attr(dx, "Labels"), attr(dy, "Labels")))
    stop("dx and dy labels differ")
  n <- attr(dx, "Size")
  dm <- as.matrix(dx); ym <- as.matrix(dy)
  lt <- lower.tri(dm)
  idx <- which(lt, arr.ind = TRUE)
  i <- idx[, 1L]; j <- idx[, 2L]
  x <- dm[lt]
  y <- ym[lt]
  obs <- stats::cor(x, y, method = method)

  if (!is.null(strata)) {
    strata <- as.integer(factor(strata))
    if (length(strata) != n) stop("strata length must equal matrix size")
  }

  set.seed(seed)
  null <- numeric(nperm)
  for (b in seq_len(nperm)) {
    if (is.null(strata)) {
      perm <- sample.int(n)
    } else {
      perm <- integer(n)
      for (s in unique(strata)) {
        ii <- which(strata == s)
        perm[ii] <- sample(ii)
      }
    }
    yp <- ym[cbind(perm[i], perm[j])]
    null[b] <- stats::cor(x, yp, method = method)
  }

  p_greater <- (1 + sum(null >= obs)) / (nperm + 1)
  p_less    <- (1 + sum(null <= obs)) / (nperm + 1)
  p_two     <- min(1, 2 * min(p_greater, p_less))
  p <- switch(alternative, greater = p_greater, less = p_less, two.sided = p_two)
  mu <- mean(null); sdv <- stats::sd(null)
  list(
    n = n, n_pairs = length(x), r = obs, p = p,
    p_greater = p_greater, p_less = p_less, p_two_sided = p_two,
    null_mean = mu, null_sd = sdv,
    SES = if (is.finite(sdv) && sdv > 0) (obs - mu) / sdv else NA_real_,
    null = null
  )
}

# ---------------------------------------------------------------------------
# Mantel correlogram: phylo-distance bins, prevalence-stratified null envelope
# ---------------------------------------------------------------------------
correlogram_fast <- function(dist_list, dp, nbin = 10L, nperm = 999L,
                             strata = NULL, seed = 20261112L) {
  n <- attr(dp, "Size")
  if (!is.null(strata)) {
    strata <- as.integer(factor(strata))
    if (length(strata) != n) stop("strata length mismatch")
  }
  pm <- as.matrix(dp)
  lt <- lower.tri(pm)
  idx <- which(lt, arr.ind = TRUE)
  i <- idx[, 1L]; j <- idx[, 2L]
  pd <- pm[lt]

  br <- unique(stats::quantile(pd, probs = seq(0, 1, length.out = nbin + 1L),
                               na.rm = TRUE))
  if (length(br) < 3L) stop("too few unique phylo distances for bins")
  bin_id <- cut(pd, breaks = br, include.lowest = TRUE, labels = FALSE)
  nb <- max(bin_id, na.rm = TRUE)
  bin_mid <- tapply(pd, bin_id, mean)
  bin_n <- tabulate(bin_id, nbins = nb)

  out <- list()
  for (metric in names(dist_list)) {
    dm <- as.matrix(dist_list[[metric]])
    if (!identical(rownames(dm), rownames(pm))) stop("metric/phylo labels differ")
    y <- dm[lt]
    obs <- tapply(y, bin_id, mean)
    nmat <- matrix(NA_real_, nrow = nperm, ncol = nb)
    set.seed(seed)
    for (b in seq_len(nperm)) {
      if (is.null(strata)) {
        perm <- sample.int(n)
      } else {
        perm <- integer(n)
        for (s in unique(strata)) {
          ii <- which(strata == s)
          perm[ii] <- sample(ii)
        }
      }
      yp <- dm[cbind(perm[i], perm[j])]
      nmat[b, ] <- tapply(yp, bin_id, mean)
    }
    nm <- colMeans(nmat, na.rm = TRUE)
    ns <- apply(nmat, 2L, stats::sd, na.rm = TRUE)
    ses <- (obs - nm) / ns
    p_lo <- vapply(seq_len(nb), function(b)
      (1 + sum(nmat[, b] <= obs[b], na.rm = TRUE)) / (sum(!is.na(nmat[, b])) + 1),
      numeric(1))
    p_hi <- vapply(seq_len(nb), function(b)
      (1 + sum(nmat[, b] >= obs[b], na.rm = TRUE)) / (sum(!is.na(nmat[, b])) + 1),
      numeric(1))
    out[[metric]] <- data.table::data.table(
      metric = metric, bin = seq_len(nb),
      phylo_mid = as.numeric(bin_mid[seq_len(nb)]),
      n_pairs = as.integer(bin_n),
      obs_mean = as.numeric(obs[seq_len(nb)]),
      null_mean = as.numeric(nm),
      null_sd = as.numeric(ns),
      SES = as.numeric(ses),
      p_lower = p_lo, p_upper = p_hi,
      null_q025 = apply(nmat, 2L, stats::quantile, probs = 0.025, na.rm = TRUE),
      null_q975 = apply(nmat, 2L, stats::quantile, probs = 0.975, na.rm = TRUE)
    )
    attr(out[[metric]], "null_draws") <- nmat
  }
  out
}

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
prevalence_strata <- function(prev, k = 5L) {
  r <- rank(prev, ties.method = "first")
  as.integer(cut(r, breaks = k, labels = FALSE))
}

as_matrix_dist <- function(d) {
  m <- as.matrix(d)
  m
}

fwrite_safe <- function(x, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file, sep = "\t", ...)
}
