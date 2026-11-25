#!/usr/bin/env Rscript
# =============================================================================
# gos_allASV_phylosig_fast.R — all-ASV tree 的 host_n 系统发育信号
#
# 输入：
#   result/GOS/tree_16s_allASV/gos_all_hostn_16s_midpoint.nwk
#   result/GOS/tree_16s_allASV/annotation_16s_allASV.tsv
#   scripts/R/phylosig_fast.R
#
# 环境变量：
#   GOS_TREE_DIR  树目录
#   SUBSET        all（默认，全部 22,845 ASV）或 positive（仅 n_hosts>=1，剪枝）
#   GOS_NSIM      K 置换次数，默认 1000
#
# 输出：
#   phylosig_hostn_16s_allASV_{SUBSET}.tsv
#
# 用法：
#   OPENBLAS_NUM_THREADS=32 OMP_NUM_THREADS=32 \
#     SUBSET=all Rscript gos_allASV_phylosig_fast.R
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
})

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "<SERVER_BASE>")
TREE_DIR <- Sys.getenv("GOS_TREE_DIR",
                       unset = file.path(BASE_DIR, "result/GOS/tree_16s_allASV"))
SUBSET <- Sys.getenv("SUBSET", unset = "all")
FAST_LAMBDA <- Sys.getenv("FAST_LAMBDA", unset = "0") == "1"
NSIM <- as.integer(Sys.getenv("GOS_NSIM", unset = "1000"))
stopifnot(SUBSET %in% c("all", "positive"))

source(file.path(BASE_DIR, "scripts/R/phylosig_fast.R"))

tree <- ape::read.tree(file.path(TREE_DIR, "gos_all_hostn_16s_midpoint.nwk"))
ann <- data.table::fread(file.path(TREE_DIR, "annotation_16s_allASV.tsv"))
stopifnot(ape::is.rooted(tree))
stopifnot(!anyDuplicated(tree$tip.label), !anyDuplicated(ann$taxid))
stopifnot(setequal(tree$tip.label, ann$taxid))

if (SUBSET == "positive") {
  keep <- ann[taxid %in% tree$tip.label & n_hosts >= 1, taxid]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  ann <- ann[taxid %in% tree$tip.label]
}

ann <- ann[match(tree$tip.label, taxid)]
stopifnot(all(ann$taxid == tree$tip.label), !anyNA(ann$n_hosts))
x <- stats::setNames(as.numeric(ann$n_hosts), ann$taxid)

p <- prepare_vcv(tree, x)
if (FAST_LAMBDA) {
  fit_lambda <- pagel_lambda_pruning(tree, x)
} else {
  fit_lambda <- pagel_lambda_fast(p$C, p$x)
}
fit_K <- blomberg_k_fast(p$C, p$x, nsim = NSIM, seed = 20261107L)

sig <- data.table::data.table(
  trait = "n_hosts (raw)",
  subset = SUBSET,
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
out_file <- file.path(TREE_DIR,
                      sprintf("phylosig_hostn_16s_allASV_%s.tsv", SUBSET))
data.table::fwrite(sig, out_file, sep = "\t")

cat(sprintf("[phylosig-allASV] subset=%s tips=%d | lambda=%.6f (P=%.4g) | K=%.6f (P=%.4g)\n",
            SUBSET, length(tree$tip.label), fit_lambda$lambda, fit_lambda$P,
            fit_K$K, fit_K$P))
print(sig)
cat("DONE.\n")
