#!/usr/bin/env Rscript
# Recompute Track 1 MRM rows after fixing select_candidates() to return meta.
# The previous run passed sel$meta = NULL to dist_block(), so dhost/dsrc were
# all zero and the "extended" models were effectively the primary model.
#
# This script rebuilds candidates, computes correct dhost/dsrc, and rewrites
# the MRM tables for all Track 1 regions (16S + ITS).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ape); library(vegan); library(ecodist)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/mapped_utils.R"))

NPERM_MRM <- as.integer(Sys.getenv("NPERM_MRM", "1999"))
TAB_DIR <- file.path(OUT_ROOT, "result_table", "cross_region")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

run_ds <- function(ds) {
  cfg <- MAPPED_CFG[[ds]]
  id_map <- if (ds == "ITS")
    file.path(OUT_ROOT, "Temp_data/trees/its_ref_dup_map.tsv") else NULL
  cache <- load_mapped_cache(ds)
  regions <- if (ds == "16S") c("V3_V4", "V4", "V5_V6", "V5_V7") else c("ITS1", "ITS2")
  rows <- list()
  for (reg in regions) {
    cat(sprintf("\n===== MRM %s / %s =====\n", ds, reg)); flush.console()
    rs <- region_subset(cache, reg)
    sel <- select_candidates(rs$mat, rs$meta, prev_min = 20L,
                             host_min = 5L, src_min = 3L)
    mat <- sel$mat
    db <- dist_block(mat, sel$meta)
    dp <- phylo_dist(file.path(BASE_DIR, cfg$tree), sel$ids, id_map = id_map)
    for (metric in c("J", "B")) {
      d <- db[[c(J = "dj", B = "db")[metric]]]
      m1 <- mrm_fit(d, dp, db$dpv, db$dab, nperm = NPERM_MRM)
      m1[, `:=`(dataset = ds, region = reg, metric = metric,
                model = "primary",
                extra_terms = attr(m1, "extra_terms"))]
      m2 <- mrm_fit(d, dp, db$dpv, db$dab,
                    extra = list(dhost = db$dhost, dsrc = db$dsrc),
                    nperm = NPERM_MRM)
      m2[, `:=`(dataset = ds, region = reg, metric = metric,
                model = "extended",
                extra_terms = attr(m2, "extra_terms"))]
      rows[[length(rows) + 1L]] <- rbind(m1, m2)
      cat(sprintf("  [%s] primary dp p=%.4g; extended dp p=%.4g; extra=%s\n",
                  metric, m1[term == "dp", p], m2[term == "dp", p],
                  attr(m2, "extra_terms")))
      flush.console()
    }
    rm(rs, sel, mat, db, dp); invisible(gc())
  }
  out <- rbindlist(rows, fill = TRUE)
  fwrite_safe(out, file.path(TAB_DIR,
             sprintf("track1_%s_mrm.tsv", tolower(ds))))
  out
}

res16 <- run_ds("16S")
resit <- run_ds("ITS")
cat("\nDONE 19_fix_track1_mrm\n")
print(rbind(res16[term == "dp"], resit[term == "dp"]))
