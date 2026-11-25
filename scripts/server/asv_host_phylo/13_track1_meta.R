#!/usr/bin/env Rscript
# Recompute Track 1 meta-analysis from saved Mantel tables.
# Kept separate so the meta step can be rerun without recomputing permutations.

source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/mapped_utils.R"))

TAB_DIR <- file.path(OUT_ROOT, "result_table", "cross_region")
for (ds in c("16S", "ITS")) {
  tag <- tolower(ds)
  f <- file.path(TAB_DIR, sprintf("track1_%s_mantel.tsv", tag))
  if (!file.exists(f)) next
  m <- fread(f)
  region_sets <- if (ds == "16S") {
    list(core3 = c("V3_V4", "V4", "V5_V7"),
         all4 = c("V3_V4", "V4", "V5_V6", "V5_V7"))
  } else {
    list(all = unique(m$region))
  }
  rows <- list()
  for (met in c("J", "B")) {
    for (set_name in names(region_sets)) {
      regs <- region_sets[[set_name]]
      sub <- m[metric == met & region %in% regs,
               .(region, metric, r, n_bins)]
      mm <- meta_mantel(sub)
      if (!is.null(mm)) {
        mm[, `:=`(dataset = ds, metric = met, region_set = set_name)]
        rows[[length(rows) + 1L]] <- mm
      }
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  fwrite_safe(out, file.path(TAB_DIR, sprintf("track1_%s_meta.tsv", tag)))
  cat("\n=== meta", ds, "===\n"); print(out)
}
