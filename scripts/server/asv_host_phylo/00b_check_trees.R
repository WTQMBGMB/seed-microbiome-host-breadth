#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(ape); library(data.table)})
BASE <- "<SERVER_BASE>"
OUT <- file.path(BASE, "result/ASV_host_phylo")
regs <- c("v4","v5_v7","its1"); fls <- c("main20_5_3","loose10_3_2","strict30_5_3")
rows <- list()
for (reg in regs) for (fl in fls) {
  tr <- read.tree(file.path(OUT, "Temp_data/trees", sprintf("%s_%s.nwk", reg, fl)))
  dt <- fread(file.path(OUT, "Temp_data/candidates", sprintf("%s_%s_candidates.tsv", reg, fl)))
  pd <- cophenetic(tr)
  rows[[length(rows)+1L]] <- data.table(
    region=reg, filter=fl, n_candidates=nrow(dt), n_tips=length(tr$tip.label),
    all_candidates_in_tree=all(dt$ASV_ID %in% tr$tip.label),
    n_missing=sum(!dt$ASV_ID %in% tr$tip.label),
    n_extra=sum(!tr$tip.label %in% dt$ASV_ID),
    rooted=is.rooted(tr), has_brlen=!is.null(tr$edge.length),
    n_neg_brlen=if (is.null(tr$edge.length)) NA_integer_ else sum(tr$edge.length < 0),
    n_zero_brlen=if (is.null(tr$edge.length)) NA_integer_ else sum(tr$edge.length == 0),
    cophenetic_median=median(pd[upper.tri(pd)]), cophenetic_q90=as.numeric(quantile(pd[upper.tri(pd)], .9))
  )
}
res <- rbindlist(rows)
fwrite(res, file.path(OUT, "result_table", "tree_qc.tsv"), sep="\t")
print(res)
