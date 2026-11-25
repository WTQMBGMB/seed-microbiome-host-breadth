#!/usr/bin/env Rscript
# =============================================================================
# mapped_utils.R — helpers for cross-region mapped-reference-bin analyses
#
# Unit: mapped reference bin = row ID present in Merged_*_{Silva,Unite}.fa
# (vsearch --usearch_global at 97% identity; one tip per reference sequence).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(vegan)
  library(ecodist)
})

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "<SERVER_BASE>")
OUT_ROOT <- file.path(BASE_DIR, "result/ASV_host_phylo")
source(file.path(OUT_ROOT, "scripts/asv_host_phylo_utils.R"))

MAPPED_CFG <- list(
  `16S` = list(
    table = "filtered_merge/16S/Merged_16S_table_filtered.txt",
    ref   = "filtered_merge/16S/Merged_16S_Silva_filtered.fa",
    tax   = "filtered_merge/16S/Merged_16S_taxonomy_filtered.txt",
    meta  = "standard_merge/metadata_16S.txt",
    tree  = "result/GOS/tree_16s_allASV/all_asv_tree_itol.nwk"
  ),
  ITS = list(
    table = "filtered_merge/ITS/Merged_ITS_table_filtered.txt",
    ref   = "filtered_merge/ITS/Merged_ITS_Unite_filtered.fa",
    tax   = "filtered_merge/ITS/Merged_ITS_taxonomy_filtered.txt",
    meta  = "standard_merge/metadata_ITS.txt",
    tree  = "result/ASV_host_phylo/Temp_data/trees/its_ref_unite.nwk"
  )
)

norm_id <- function(x) {
  x <- gsub("-", ".", x, fixed = TRUE)
  ifelse(grepl("^[0-9]", x), paste0("X", x), x)
}

read_fa_ids <- function(f) {
  x <- readLines(f, warn = FALSE)
  sub("^>([^ ]+).*$", "\\1", x[startsWith(x, ">")])
}

# ---------------------------------------------------------------------------
# Load merged mapped table once and cache matrix + aligned metadata + taxonomy
# ---------------------------------------------------------------------------
load_mapped_cache <- function(ds, force = FALSE, verbose = TRUE) {
  stopifnot(ds %in% names(MAPPED_CFG))
  cfg <- MAPPED_CFG[[ds]]
  cache_file <- file.path(OUT_ROOT, "Temp_data/rds",
                          sprintf("mapped_cache_%s.rds", ds))
  if (file.exists(cache_file) && !force) {
    if (verbose) cat(sprintf("[cache] load %s\n", basename(cache_file)))
    return(readRDS(cache_file))
  }

  if (verbose) cat(sprintf("[cache] reading %s ...\n", cfg$table))
  ref_ids <- read_fa_ids(file.path(BASE_DIR, cfg$ref))
  tab <- fread(file.path(BASE_DIR, cfg$table), sep = "\t", header = TRUE,
               check.names = FALSE)
  setnames(tab, 1L, "ASV_ID")
  keep <- as.character(tab$ASV_ID) %in% ref_ids
  tab <- tab[keep]
  ids <- as.character(tab$ASV_ID)
  mat <- as.matrix(tab[, -1, with = FALSE])
  rownames(mat) <- ids
  rm(tab); invisible(gc())

  sample_raw <- colnames(mat)
  sample_norm <- norm_id(sample_raw)
  if (anyDuplicated(sample_norm)) {
    stop(sprintf("%s: duplicated sample names after normalization", ds))
  }

  meta <- fread(file.path(BASE_DIR, cfg$meta), sep = "\t")
  meta[, SampleID := norm_id(SampleID)]
  idx <- match(sample_norm, meta$SampleID)
  if (anyNA(idx)) {
    stop(sprintf("%s: %d mapped-table samples not matched to metadata",
                 ds, sum(is.na(idx))))
  }
  meta <- meta[idx]
  meta[, species := as.character(species)]
  meta[, source := as.character(source)]
  meta[, Region := as.character(Region)]

  tax <- unique(fread(file.path(BASE_DIR, cfg$tax), sep = "\t"), by = "taxid")
  tax <- tax[match(ids, tax$taxid)]
  tax_cols <- intersect(c("taxid", "Kingdom", "Phylum", "Class", "Order",
                          "Family", "Genus", "Species"), names(tax))
  tax <- tax[, ..tax_cols]

  cache <- list(ds = ds, ids = ids, mat = mat, meta = meta, tax = tax,
                sample_raw = sample_raw, sample_norm = sample_norm,
                ref_ids = ref_ids)
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, cache_file)
  if (verbose) cat(sprintf("[cache] saved %s: %d mapped bins x %d samples\n",
                           basename(cache_file), nrow(mat), ncol(mat)))
  cache
}

# ---------------------------------------------------------------------------
# Region subset and candidate filter
# ---------------------------------------------------------------------------
region_subset <- function(cache, region) {
  cols <- which(cache$meta$Region == region)
  if (!length(cols)) stop(sprintf("region %s not found", region))
  list(mat = cache$mat[, cols, drop = FALSE],
       meta = cache$meta[cols])
}

breadth_stats <- function(mat, meta) {
  sp <- meta$species
  so <- meta$source
  sp[!nzchar(sp)] <- NA_character_
  so[!nzchar(so)] <- NA_character_
  prev <- rowSums(mat > 0)
  col_tot <- colSums(mat)
  col_tot[!is.finite(col_tot) | col_tot <= 0] <- 1
  rel <- sweep(mat, 2L, col_tot, "/")
  mean_rel <- rowMeans(rel)
  nhost <- vapply(seq_len(nrow(mat)), function(i)
    uniqueN(sp[mat[i, ] > 0], na.rm = TRUE), integer(1))
  nsrc <- vapply(seq_len(nrow(mat)), function(i)
    uniqueN(so[mat[i, ] > 0], na.rm = TRUE), integer(1))
  ids <- rownames(mat)
  names(prev) <- names(mean_rel) <- names(nhost) <- names(nsrc) <- ids
  list(prev = prev, mean_rel = mean_rel, nhost = nhost, nsrc = nsrc,
       rel = rel)
}

select_candidates <- function(mat, meta, prev_min = 20L,
                              host_min = 5L, src_min = 3L) {
  bs <- breadth_stats(mat, meta)
  sp <- meta$species; sp[!nzchar(sp)] <- NA_character_
  so <- meta$source; so[!nzchar(so)] <- NA_character_
  host_thr <- min(host_min, uniqueN(sp, na.rm = TRUE))
  src_thr <- min(src_min, uniqueN(so, na.rm = TRUE))
  keep <- bs$prev >= prev_min & bs$nhost >= host_thr & bs$nsrc >= src_thr
  list(
    mat = mat[keep, , drop = FALSE],
    meta = meta,
    ids = rownames(mat)[keep],
    stats = data.table(
      ASV_ID = rownames(mat),
      prevalence = bs$prev,
      n_hosts = bs$nhost,
      n_sources = bs$nsrc,
      keep = keep
    ),
    host_thr = host_thr, src_thr = src_thr, prev_min = prev_min
  )
}

# ---------------------------------------------------------------------------
# Distribution distances and covariates
# ---------------------------------------------------------------------------
dist_block <- function(mat, meta) {
  ids <- rownames(mat)
  bs <- breadth_stats(mat, meta)
  rel <- bs$rel
  dj <- vegan::vegdist(mat > 0, method = "jaccard", binary = TRUE)
  if (requireNamespace("parallelDist", quietly = TRUE)) {
    db <- parallelDist::parDist(rel, method = "bray",
                                threads = as.integer(Sys.getenv("DIST_THREADS", "4")))
    if (is.null(attr(db, "Labels"))) attr(db, "Labels") <- ids
  } else {
    db <- vegan::vegdist(rel, method = "bray")
  }
  if (!identical(attr(dj, "Labels"), ids)) stop("Jaccard labels misaligned")
  if (!identical(attr(db, "Labels"), ids)) stop("Bray labels misaligned")

  dpv <- stats::as.dist(abs(outer(log10(bs$prev), log10(bs$prev), "-")))
  dab <- stats::as.dist(abs(outer(log10(bs$mean_rel + 1e-6),
                                  log10(bs$mean_rel + 1e-6), "-")))
  dhost <- stats::as.dist(abs(outer(bs$nhost, bs$nhost, "-")))
  so <- meta$source; so[!nzchar(so)] <- NA_character_
  if (uniqueN(so, na.rm = TRUE) >= 2L) {
    so_fac <- factor(so)
    src_model <- stats::model.matrix(~ 0 + so_fac)
    src_pa <- (mat > 0) %*% src_model
    dsrc <- vegan::vegdist(src_pa > 0, method = "jaccard", binary = TRUE)
  } else {
    dsrc <- stats::as.dist(matrix(0, nrow(mat), nrow(mat),
                                  dimnames = list(ids, ids)))
  }
  for (dd in list(dpv, dab, dhost, dsrc))
    if (!identical(attr(dd, "Labels"), ids)) stop("covariate labels misaligned")
  list(dj = dj, db = db, prev = bs$prev, mean_rel = bs$mean_rel,
       nhost = bs$nhost, nsrc = bs$nsrc,
       dpv = dpv, dab = dab, dhost = dhost, dsrc = dsrc)
}

# ---------------------------------------------------------------------------
# Reference-tree patristic distance
# ---------------------------------------------------------------------------
phylo_dist <- function(tree_file, ids, id_map = NULL) {
  tr <- ape::read.tree(tree_file)
  if (!is.null(id_map)) {
    map <- data.table::fread(id_map)
    if (!all(c("id", "rep") %in% names(map)))
      stop("id_map must have columns id and rep")
    reps <- map$rep[match(ids, map$id)]
    if (anyNA(reps)) {
      stop(sprintf("id_map lacks %d requested ids; first: %s",
                   sum(is.na(reps)), paste(head(ids[is.na(reps)], 3),
                                           collapse = ", ")))
    }
    need <- unique(reps)
    miss <- setdiff(need, tr$tip.label)
    if (length(miss)) {
      stop(sprintf("tree lacks %d mapped representatives; first: %s",
                   length(miss), paste(head(miss, 3), collapse = ", ")))
    }
    if (length(tr$tip.label) != length(need)) tr <- ape::keep.tip(tr, need)
    dpm <- stats::cophenetic(tr)[reps, reps, drop = FALSE]
    dimnames(dpm) <- list(ids, ids)
    dp <- stats::as.dist(dpm)
  } else {
    miss <- setdiff(ids, tr$tip.label)
    if (length(miss)) {
      stop(sprintf("tree lacks %d of %d requested tips; first: %s",
                   length(miss), length(ids), paste(head(miss, 3), collapse = ", ")))
    }
    if (length(tr$tip.label) != length(ids)) tr <- ape::keep.tip(tr, ids)
    dpm <- stats::cophenetic(tr)[ids, ids, drop = FALSE]
    dp <- stats::as.dist(dpm)
  }
  if (!identical(attr(dp, "Labels"), ids)) stop("phylo labels misaligned")
  dp
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
mantel_ecodist <- function(d, dp, nperm = 9999L) {
  e <- ecodist::mantel(d ~ dp, nperm = nperm)
  list(r = as.numeric(e["mantelr"]), p = as.numeric(e["pval1"]),
       p_two = as.numeric(e["pval3"]),
       CI_low = as.numeric(e["llim.2.5%"]),
       CI_high = as.numeric(e["ulim.97.5%"]))
}

mrm_fit <- function(d, dp, dpv, dab, extra = NULL, nperm = 1999L) {
  lt <- function(dd) { m <- as.matrix(dd); m[lower.tri(m)] }
  base_vecs <- list(dp = lt(dp), dpv = lt(dpv), dab = lt(dab))
  used_extra <- character(0)
  if (!is.null(extra) && length(extra)) {
    for (nm in names(extra)) {
      if (!inherits(extra[[nm]], "dist")) next
      v <- lt(extra[[nm]])
      if (stats::sd(v) == 0) next
      cors <- vapply(base_vecs, function(b)
        abs(stats::cor(v, b, use = "complete.obs")), numeric(1))
      if (any(!is.finite(cors)) || max(cors) > 0.995) next
      if (length(used_extra)) {
        cors2 <- vapply(used_extra, function(nm2)
          abs(stats::cor(v, lt(extra[[nm2]]), use = "complete.obs")), numeric(1))
        if (any(!is.finite(cors2)) || max(cors2) > 0.995) next
      }
      used_extra <- c(used_extra, nm)
      base_vecs[[nm]] <- v
      assign(nm, extra[[nm]])
    }
  }
  if (!length(used_extra)) {
    form <- stats::as.formula("d ~ dp + dpv + dab")
  } else {
    form <- stats::as.formula(
      paste("d ~ dp + dpv + dab +", paste(used_extra, collapse = " + ")))
  }
  m <- ecodist::MRM(form, nperm = nperm)
  co <- m$coef
  out <- data.table(term = rownames(co), coef = as.numeric(co[, 1L]),
                    p = as.numeric(co[, 2L]),
                    r_squared = as.numeric(m$r.squared[1L]),
                    r_squared_p = as.numeric(m$r.squared[2L]))
  attr(out, "extra_terms") <- paste(used_extra, collapse = "+")
  out
}

close_far_fast <- function(d, dp, strata = NULL, nperm = 999L, seed = 20261112L) {
  if (!identical(attr(d, "Labels"), attr(dp, "Labels")))
    stop("close_far labels differ")
  dm <- as.matrix(d); pm <- as.matrix(dp)
  lt <- lower.tri(dm)
  y <- dm[lt]; x <- pm[lt]
  q_close <- stats::quantile(x, 0.10, na.rm = TRUE)
  q_far   <- stats::quantile(x, 0.50, na.rm = TRUE)
  close <- x <= q_close
  far   <- x >= q_far
  if (!any(close) || !any(far)) return(NULL)
  obs <- mean(y[far]) - mean(y[close])
  n <- nrow(dm)
  idx <- which(lt, arr.ind = TRUE)
  ii <- idx[, 1L]; jj <- idx[, 2L]
  if (!is.null(strata)) strata <- as.integer(factor(strata))
  set.seed(seed)
  null <- numeric(nperm)
  for (b in seq_len(nperm)) {
    if (is.null(strata)) perm <- sample.int(n) else {
      perm <- integer(n)
      for (s in unique(strata)) {
        ii_s <- which(strata == s)
        perm[ii_s] <- sample(ii_s)
      }
    }
    yp <- dm[cbind(perm[ii], perm[jj])]
    null[b] <- mean(yp[far]) - mean(yp[close])
  }
  p <- (1 + sum(null >= obs)) / (nperm + 1)
  mu <- mean(null); sdv <- stats::sd(null)
  list(contrast = obs, p = p, null_mean = mu, null_sd = sdv,
       SES = if (is.finite(sdv) && sdv > 0) (obs - mu) / sdv else NA_real_,
       n_close = sum(close), n_far = sum(far))
}

# ---------------------------------------------------------------------------
# Random-effects meta-analysis of Mantel r (Fisher z)
# ---------------------------------------------------------------------------
meta_mantel <- function(dt) {
  dt <- dt[is.finite(r) & n_bins > 4L]
  if (nrow(dt) < 2L) return(NULL)
  yi <- atanh(dt$r)
  vi <- 1 / (dt$n_bins - 3)
  fit_re <- metafor::rma.uni(yi = yi, vi = vi, method = "REML")
  fit_fe <- metafor::rma.uni(yi = yi, vi = vi, method = "FE")
  ci <- function(fit) {
    b <- as.numeric(fit$beta[1L]); se <- fit$se[1L]
    c(tanh(b), tanh(b - 1.96 * se), tanh(b + 1.96 * se))
  }
  re <- ci(fit_re); fe <- ci(fit_fe)
  data.table(
    k = nrow(dt),
    r_RE = re[1], CI_low_RE = re[2], CI_high_RE = re[3],
    p_RE = as.numeric(fit_re$pval),
    r_FE = fe[1], CI_low_FE = fe[2], CI_high_FE = fe[3],
    p_FE = as.numeric(fit_fe$pval),
    tau2 = as.numeric(fit_re$tau2),
    I2 = as.numeric(fit_re$I2),
    Q = as.numeric(fit_re$QE),
    Q_p = as.numeric(fit_re$QEp)
  )
}
