#!/usr/bin/env Rscript
# =============================================================================
# gos_hostn_pos_phylosig_fast.R — phylogenetic signal on the tree built from
# host breadth >= 1 ASVs (de novo, 16S).
#
# Inputs (GOS_TREE_DIR):
#   gos_hostn_pos_16s_midpoint.nwk
#   annotation_16s_hostn_pos.tsv
#   scripts/R/phylosig_fast.R
#
# Env:
#   GOS_TREE_DIR, GOS_NSIM (default 1000), FAST_LAMBDA (default 1)
#
# Output:
#   phylosig_hostn_16s_hostn_pos.tsv
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
})

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "<SERVER_BASE>2")
TREE_DIR <- Sys.getenv("GOS_TREE_DIR", unset = BASE_DIR)
FAST_LAMBDA <- Sys.getenv("FAST_LAMBDA", unset = "1") == "1"
NSIM <- as.integer(Sys.getenv("GOS_NSIM", unset = "1000"))

source(file.path(BASE_DIR, "scripts/R/phylosig_fast.R"))

tree <- ape::read.tree(file.path(TREE_DIR, "gos_hostn_pos_16s_midpoint.nwk"))
ann <- data.table::fread(file.path(TREE_DIR, "annotation_16s_hostn_pos.tsv"))
stopifnot(ape::is.rooted(tree))
stopifnot(!anyDuplicated(tree$tip.label), !anyDuplicated(ann$taxid))
stopifnot(setequal(tree$tip.label, ann$taxid))

ann <- ann[match(tree$tip.label, taxid)]
stopifnot(all(ann$taxid == tree$tip.label), !anyNA(ann$host_breadth))
stopifnot(all(ann$host_breadth >= 1))
x <- stats::setNames(as.numeric(ann$host_breadth), ann$taxid)

p <- prepare_vcv(tree, x)
if (FAST_LAMBDA) {
  fit_lambda <- pagel_lambda_pruning(tree, x)
} else {
  fit_lambda <- pagel_lambda_fast(p$C, p$x)
}
fit_K <- blomberg_k_fast(p$C, p$x, nsim = NSIM, seed = 20261107L)

sig <- data.table::data.table(
  trait = "host breadth (raw)",
  subset = "host breadth >= 1",
  method = c("Pagel_lambda", "Blomberg_K"),
  estimate = c(unname(fit_lambda$lambda), unname(fit_K$K)),
  P = c(unname(fit_lambda$P), unname(fit_K$P)),
  n_tips = length(tree$tip.label),
  rooting = "midpoint",
  note = c(if (FAST_LAMBDA)
             "likelihood-ratio test; O(n) pruning (phylolm BM on lambda-transformed tree)"
           else
             "likelihood-ratio test; fast eigen implementation",
           sprintf("%d permutations; vectorized BLAS%s", NSIM,
                   if (fit_K$ridge > 0)
                     sprintf("; VCV ridge=%.3g (zero-length tips)",
                             fit_K$ridge) else "")))
out_file <- file.path(TREE_DIR, "phylosig_hostn_16s_hostn_pos.tsv")
data.table::fwrite(sig, out_file, sep = "\t")

cat(sprintf("[phylosig-hostn-pos] tips=%d | lambda=%.6f (P=%.4g) | K=%.6f (P=%.4g)\n",
            length(tree$tip.label), fit_lambda$lambda, fit_lambda$P,
            fit_K$K, fit_K$P))
print(sig)
cat("DONE.\n")
