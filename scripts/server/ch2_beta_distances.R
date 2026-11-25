#!/usr/bin/env Rscript
# Per-region Bray-Curtis distances / 区段 Bray-Curtis 距离
# Usage: Rscript ch2_beta_distances.R <work_dir>
# Reads <work_dir>/standard_region/<reg>/css_<reg>_asvtab.txt and writes
# <work_dir>/bray_dist_<reg>.rds for v4, v5_v7 and its1.
suppressPackageStartupMessages({
  library(data.table)
  library(Rfast)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: ch2_beta_distances.R <work_dir>")
work <- args[1]

for (reg in c("v4", "v5_v7", "its1")) {
  f <- file.path(work, "standard_region", reg, paste0("css_", reg, "_asvtab.txt"))
  raw <- fread(f, sep = "\t", header = TRUE)
  mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(mat) <- as.character(raw[[1]])
  otu <- t(mat)
  rm(raw, mat); invisible(gc())
  rs <- rowSums(otu)
  d <- Rfast::Dist(otu, method = "manhattan", parallel = TRUE,
                   result = "matrix") / outer(rs, rs, "+")
  rownames(d) <- colnames(d) <- rownames(otu)
  diag(d) <- 0
  stopifnot(all(is.finite(d)))
  saveRDS(as.dist(d), file.path(work, paste0("bray_dist_", reg, ".rds")))
  cat(reg, ":", nrow(d), "samples\n")
}
