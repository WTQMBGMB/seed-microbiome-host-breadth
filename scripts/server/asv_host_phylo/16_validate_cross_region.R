#!/usr/bin/env Rscript
# Validation checks for the cross-region mapped-bin analysis.

suppressPackageStartupMessages({ library(data.table); library(ape) })
BASE <- Sys.getenv("BASE_DIR", "<SERVER_BASE>")
TAB <- file.path(BASE, "result/ASV_host_phylo/result_table/cross_region")
OUT <- file.path(TAB, "validation_cross_region.tsv")

checks <- list()
add <- function(check, value, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(
    check = check, value = as.character(value),
    pass = as.logical(pass), detail = detail)
}

for (ds in c("16S", "ITS")) {
  tag <- tolower(ds)
  f1 <- file.path(TAB, sprintf("track1_%s_mantel.tsv", tag))
  f2 <- file.path(TAB, sprintf("track2_%s_mantel.tsv", tag))
  f1m <- file.path(TAB, sprintf("track1_%s_mrm.tsv", tag))
  f2m <- file.path(TAB, sprintf("track2_%s_mrm.tsv", tag))
  f1c <- file.path(TAB, sprintf("track1_%s_candidates.tsv", tag))
  f1g <- file.path(TAB, sprintf("track1_%s_correlogram.tsv", tag))
  f1f <- file.path(TAB, sprintf("track1_%s_close_far.tsv", tag))
  f2p <- file.path(TAB, sprintf("track2_%s_pairs.tsv", tag))

  if (file.exists(f1)) {
    m1 <- fread(f1)
    nreg <- uniqueN(m1$region); nmet <- uniqueN(m1$metric)
    add(paste0(ds, "_track1_rows"), nrow(m1), nrow(m1) == nreg * nmet)
    add(paste0(ds, "_track1_no_na"),
        sum(is.na(m1$r) | is.na(m1$p_standard) | is.na(m1$CI_low) |
            is.na(m1$CI_high)),
        !anyNA(m1$r) && !anyNA(m1$p_standard) && !anyNA(m1$CI_low) &&
          !anyNA(m1$CI_high))
    add(paste0(ds, "_track1_r_range"),
        paste(range(m1$r), collapse = ";"),
        all(m1$r >= -1 & m1$r <= 1))
    add(paste0(ds, "_track1_p_range"),
        paste(range(m1$p_standard), collapse = ";"),
        all(m1$p_standard >= 0 & m1$p_standard <= 1))
    if (file.exists(f1c)) {
      c1 <- fread(f1c)
      add(paste0(ds, "_track1_candidate_counts"),
          paste(c1[, .N, by = region]$N, collapse = "/"),
          all(c1[, .N, by = region]$N ==
                m1[order(region), .(n = unique(n_bins)), by = region]$n))
    }
    if (file.exists(f1g)) {
      cg <- fread(f1g)
      add(paste0(ds, "_track1_correlogram_rows"), nrow(cg),
          nrow(cg) == nreg * nmet * 10)
    }
    if (file.exists(f1f)) {
      cf <- fread(f1f)
      add(paste0(ds, "_track1_close_far_rows"), nrow(cf),
          nrow(cf) == nreg * nmet)
    }
    if (file.exists(f1m)) {
      mrm <- fread(f1m)
      nfit <- uniqueN(paste(mrm$region, mrm$metric, mrm$model))
      add(paste0(ds, "_track1_mrm_fits"), nfit, nfit == nreg * nmet * 2)
      add(paste0(ds, "_track1_mrm_dp_terms"), sum(mrm$term == "dp"),
          sum(mrm$term == "dp") == nreg * nmet * 2)
    }
  }

  if (file.exists(f2)) {
    m2 <- fread(f2)
    npair <- uniqueN(paste(m2$region_a, m2$region_b))
    add(paste0(ds, "_track2_rows"), nrow(m2), nrow(m2) == npair * 2)
    add(paste0(ds, "_track2_no_na"),
        sum(is.na(m2$r) | is.na(m2$p_standard) | is.na(m2$CI_low) |
            is.na(m2$CI_high)),
        !anyNA(m2$r) && !anyNA(m2$p_standard) && !anyNA(m2$CI_low) &&
          !anyNA(m2$CI_high))
    add(paste0(ds, "_track2_r_range"),
        paste(range(m2$r), collapse = ";"),
        all(m2$r >= -1 & m2$r <= 1))
    if (file.exists(f2p)) {
      p2 <- fread(f2p)
      add(paste0(ds, "_track2_pair_rows"), nrow(p2), nrow(p2) == npair)
      add(paste0(ds, "_track2_pair_counts"),
          paste(p2$n_pass_pair, collapse = "/"),
          all(p2$n_pass_pair >= 20))
    }
    if (file.exists(f2m)) {
      mrm <- fread(f2m)
      nfit <- uniqueN(paste(mrm$region_a, mrm$region_b, mrm$metric,
                            mrm$model))
      add(paste0(ds, "_track2_mrm_fits"), nfit, nfit == npair * 2 * 2)
      add(paste0(ds, "_track2_mrm_dp_terms"), sum(mrm$term == "dp"),
          sum(mrm$term == "dp") == npair * 2 * 2)
    }
  }
}

# tree QC
tr16 <- file.path(BASE, "result/GOS/tree_16s_allASV/all_asv_tree_itol.nwk")
if (file.exists(tr16)) {
  tr <- read.tree(tr16)
  add("tree_16s_tips", length(tr$tip.label), length(tr$tip.label) == 22845)
  add("tree_16s_negative_branches", sum(tr$edge.length < 0),
      sum(tr$edge.length < 0) == 0)
}
trits <- file.path(BASE, "result/ASV_host_phylo/Temp_data/trees/its_ref_unite.nwk")
if (file.exists(trits)) {
  tr <- read.tree(trits)
  add("tree_its_tips", length(tr$tip.label), length(tr$tip.label) == 4918)
  add("tree_its_negative_branches", sum(tr$edge.length < 0),
      sum(tr$edge.length < 0) == 0)
}

res <- rbindlist(checks)
fwrite(res, OUT, sep = "\t")
cat("Validation:", sum(res$pass), "/", nrow(res), "checks passed\n")
print(res[pass == FALSE])
cat("Written:", OUT, "\n")
