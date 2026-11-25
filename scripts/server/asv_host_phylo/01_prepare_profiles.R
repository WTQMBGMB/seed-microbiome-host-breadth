#!/usr/bin/env Rscript
# =============================================================================
# 01_prepare_profiles.R — build main-filter ASV profiles and distance objects
#
# Outputs:
#   Temp_data/rds/cache_{region}.rds              (loose ASV matrix + metadata)
#   Temp_data/rds/profiles_{region}_main20_5_3.rds
#   result_table/main_profile_summary.tsv
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ape); library(vegan)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/asv_host_phylo_utils.R"))

regions <- c("v4", "v5_v7", "its1")
filter  <- "main20_5_3"
rows <- list()

for (reg in regions) {
  cache <- load_region_cache(reg, verbose = TRUE)
  cand <- data.table::fread(cand_file(reg, filter), sep = "\t")
  ids <- intersect(cand$ASV_ID, cache$ids)
  if (length(ids) != nrow(cand))
    cat(sprintf("[%s] WARNING: %d candidate ASVs absent from cache; using %d\n",
                reg, nrow(cand) - length(ids), length(ids)))
  mat <- cache$mat[ids, , drop = FALSE]
  tax <- cache$tax[match(ids, cache$ids), ]
  prof <- build_profile_core(
    mat = mat, meta = cache$meta, tree_file = tree_file(reg, filter),
    region = reg, filter_label = filter, tax = tax
  )
  saveRDS(prof, profile_file(reg, filter))
  cat(sprintf("[%s] saved %s\n", reg, basename(profile_file(reg, filter))))

  rows[[length(rows) + 1L]] <- data.table::data.table(
    region = reg, filter = filter, n_asv = length(ids),
    n_samples = ncol(mat), n_sources = data.table::uniqueN(cache$source),
    n_host_species = data.table::uniqueN(cache$species),
    n_host_level = length(prof$host_ids),
    host_level_min_samples = min(prof$host_n_samples),
    prev_median = median(prof$prev), prev_max = max(prof$prev),
    mean_rel_median = median(prof$mean_rel)
  )
  rm(mat, prof); invisible(gc())
}

summ <- data.table::rbindlist(rows)
fwrite_safe(summ, file.path(OUT_ROOT, "result_table", "main_profile_summary.tsv"))
print(summ)
cat("DONE 01_prepare_profiles\n")
