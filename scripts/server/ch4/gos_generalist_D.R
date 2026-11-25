#!/usr/bin/env Rscript
# =============================================================================
# gos_generalist_D.R — Fritz & Purvis D for binary generalist trait
#
# 与 λ/K 相同的三套 tip 集合：
#   old tree (host_n>=1), all-ASV tree (all tips), all-ASV tree (host_n>=1)
#
# 树处理：collapse.singles() 去掉单子节点；若非二叉，multi2di(random=FALSE)
# 以零长度枝解析多分枝；随后 phylo_d_fast() 把零枝设为 eps。
# 多分枝/零枝在 all-ASV 树上较多，D 的解释需谨慎。
#
# 环境变量：
#   GOS_TREE_DIR, TREE_FILE, ANN_FILE, CAT_COL, GEN_VALUES, SUBSET,
#   LABEL, OUT_FILE, D_PERMUT（默认 1000）, D_EPS（默认 1e-8）
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
})

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "<SERVER_BASE>")
TREE_DIR <- Sys.getenv("GOS_TREE_DIR", unset = BASE_DIR)
TREE_FILE <- Sys.getenv("TREE_FILE", unset = "gos_all_hostn_16s_midpoint.nwk")
ANN_FILE <- Sys.getenv("ANN_FILE", unset = "annotation_16s_allASV.tsv")
CAT_COL <- Sys.getenv("CAT_COL", unset = "category")
GEN_VALUES <- strsplit(Sys.getenv("GEN_VALUES", unset = "generalist,Generalist"),
                       ",")[[1]]
SUBSET <- Sys.getenv("SUBSET", unset = "all")
LABEL <- Sys.getenv("LABEL", unset = "all-ASV tree (all tips)")
OUT_FILE <- Sys.getenv("OUT_FILE", unset = "phyloD_generalist.tsv")
PERMUT <- as.integer(Sys.getenv("D_PERMUT", unset = "1000"))
EPS <- as.numeric(Sys.getenv("D_EPS", unset = "1e-8"))

source(file.path(BASE_DIR, "scripts/R/phylo_d_fast.R"))

tree <- ape::read.tree(file.path(TREE_DIR, TREE_FILE))
tree$node.label <- NULL
ann <- data.table::fread(file.path(TREE_DIR, ANN_FILE))
stopifnot(!anyDuplicated(ann$taxid), CAT_COL %in% names(ann))

if (SUBSET == "positive") {
  keep <- ann[taxid %in% tree$tip.label & n_hosts >= 1, taxid]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
}
ann <- ann[taxid %in% tree$tip.label]
ann <- ann[match(tree$tip.label, taxid)]
stopifnot(all(ann$taxid == tree$tip.label))

ds <- as.integer(ann[[CAT_COL]] %in% GEN_VALUES)
names(ds) <- ann$taxid
stopifnot(sum(ds) > 0, sum(ds) < length(ds))

n_unary <- sum(lengths(split(tree$edge[, 2], tree$edge[, 1])) < 2)
n_poly <- sum(lengths(split(tree$edge[, 2], tree$edge[, 1])) > 2)
if (n_unary > 0) tree <- ape::collapse.singles(tree)
if (!ape::is.binary(tree)) tree <- ape::multi2di(tree, random = FALSE)
n_zero <- sum(tree$edge.length <= 0)
stopifnot(ape::is.binary(tree), ape::is.rooted(tree))

cat(sprintf("[D] %s | tips=%d generalist=%d (%.1f%%) | unary=%d polytomy=%d zero_edges=%d\n",
            LABEL, length(ds), sum(ds), 100 * mean(ds),
            n_unary, n_poly, n_zero))

tm <- system.time(fit <- phylo_d_fast(tree, ds, permut = PERMUT,
                                      seed = 20261108L, eps = EPS))

out <- data.table::data.table(
  trait = "generalist (binary)",
  analysis = LABEL,
  subset = SUBSET,
  method = "Fritz_Purvis_D",
  D = fit$D,
  P_random = fit$P_random,
  P_BM = fit$P_BM,
  observed = fit$observed,
  mean_random = fit$mean_random,
  mean_brownian = fit$mean_brownian,
  n_tips = length(ds),
  n_generalist = sum(ds),
  pct_generalist = round(100 * mean(ds), 2),
  polytomies_resolved = n_poly,
  zero_edges = n_zero,
  eps = EPS,
  permut = PERMUT,
  elapsed_s = round(tm[["elapsed"]], 1),
  note = "Fritz & Purvis D; P_random = P(D=1), P_BM = P(D=0)"
)
data.table::fwrite(out, file.path(TREE_DIR, OUT_FILE), sep = "\t")
print(out)
cat("DONE.\n")
