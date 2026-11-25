#!/usr/bin/env Rscript
# =============================================================================
# 00_prepare_candidates.R — ASV-host-phylogeny exploratory analysis
#
# Goal: build the candidate ASV sets used for the microbe-phylogeny vs
#       host-distribution analysis, and export per-region/per-filter FASTAs
#       for server-side candidate trees (MAFFT -> trimAl -> FastTree).
#
# Main filter (user-approved): prevalence >= 20 samples, >= 5 host species,
#                              >= 3 sources.
# Sensitivity filters: loose 10/3/2 and strict 30/5/3.
#
# Inputs (relative to Downstream_analysis):
#   standard_region/{v4,v5_v7,its1}/css_*_asvtab.txt
#   standard_merge/metadata_16S.txt | metadata_ITS.txt
#   filtered_region/{V4,V5_V7,ITS1}/*_merged_filtered.fa
#   filtered_region/{V4,V5_V7,ITS1}/*_merged_taxonomy_filtered.txt
#
# Outputs (result/ASV_host_phylo):
#   Temp_data/candidates/{region}_{filter}_candidates.tsv
#   Temp_data/tree_input/{region}_{filter}.fa
#   result_table/candidate_filter_summary.tsv
#   result_table/tree_build_manifest.tsv
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "<SERVER_BASE>")
OUT_ROOT <- file.path(BASE_DIR, "result/ASV_host_phylo")
for (d in c("Temp_data/candidates", "Temp_data/tree_input", "Temp_data/trees",
            "Temp_data/rds", "result_table", "result_figure", "logs")) {
  dir.create(file.path(OUT_ROOT, d), recursive = TRUE, showWarnings = FALSE)
}

REGIONS <- list(
  v4 = list(
    css  = "standard_region/v4/css_v4_asvtab.txt",
    meta = "standard_merge/metadata_16S.txt", rf = "V4",
    fa   = "filtered_region/V4/V4_merged_filtered.fa",
    tax  = "filtered_region/V4/V4_merged_taxonomy_filtered.txt"),
  v5_v7 = list(
    css  = "standard_region/v5_v7/css_v5_v7_asvtab.txt",
    meta = "standard_merge/metadata_16S.txt", rf = "V5_V7",
    fa   = "filtered_region/V5_V7/V5_V7_merged_filtered.fa",
    tax  = "filtered_region/V5_V7/V5_V7_merged_taxonomy_filtered.txt"),
  its1 = list(
    css  = "standard_region/its1/css_its1_asvtab.txt",
    meta = "standard_merge/metadata_ITS.txt", rf = "ITS1",
    fa   = "filtered_region/ITS1/ITS1_merged_filtered.fa",
    tax  = "filtered_region/ITS1/ITS1_merged_taxonomy_filtered.txt")
)

FILTERS <- list(
  main20_5_3   = list(prev = 20L, host = 5L, src = 3L),
  loose10_3_2  = list(prev = 10L, host = 3L, src = 2L),
  strict30_5_3 = list(prev = 30L, host = 5L, src = 3L)
)

# ---- FASTA parser (wrapped sequences; headers start with '>') --------------
read_fasta <- function(file) {
  ln <- readLines(file, warn = FALSE)
  ln <- ln[nzchar(ln)]
  hdr <- which(startsWith(ln, ">"))
  if (!length(hdr)) stop("empty FASTA: ", file)
  ids <- sub("^>", "", ln[hdr])
  ids <- sub("\\s.*$", "", ids)               # drop description if any
  ends <- c(hdr[-1] - 1L, length(ln))
  seqs <- vapply(seq_along(hdr), function(i) {
    if (hdr[i] >= ends[i]) "" else paste(ln[(hdr[i] + 1L):ends[i]], collapse = "")
  }, character(1))
  if (anyDuplicated(ids)) stop("duplicated FASTA headers in ", file)
  stats::setNames(seqs, ids)
}

summary_rows <- list()
manifest_rows <- list()

for (reg in names(REGIONS)) {
  g <- REGIONS[[reg]]
  cat(sprintf("\n=== [%s] loading CSS table ===\n", reg))
  raw <- data.table::fread(file.path(BASE_DIR, g$css), sep = "\t", header = TRUE)
  ids <- as.character(raw[[1L]])
  mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(mat) <- ids
  rm(raw); invisible(gc())

  meta <- data.table::fread(file.path(BASE_DIR, g$meta), sep = "\t")
  meta <- meta[Region == g$rf]
  common <- intersect(colnames(mat), as.character(meta$SampleID))
  if (!length(common)) {
    colnames(mat) <- sub("^X", "", colnames(mat))
    common <- intersect(colnames(mat), as.character(meta$SampleID))
    ids <- rownames(mat)
  }
  if (length(common) < 10L) stop(sprintf("[%s] too few aligned samples", reg))
  mat <- mat[, common, drop = FALSE]
  meta <- meta[match(common, meta$SampleID)]

  col_tot <- colSums(mat)
  col_tot[!is.finite(col_tot) | col_tot <= 0] <- 1
  rel <- sweep(mat, 2L, col_tot, "/")

  prev   <- rowSums(mat > 0)
  total  <- rowSums(mat)
  mean_r <- rowMeans(rel)
  sp     <- as.character(meta$species)
  so     <- as.character(meta$source)
  nhost  <- vapply(seq_len(nrow(mat)), function(i)
    data.table::uniqueN(sp[mat[i, ] > 0]), integer(1))
  nsrc   <- vapply(seq_len(nrow(mat)), function(i)
    data.table::uniqueN(so[mat[i, ] > 0]), integer(1))
  names(prev) <- names(total) <- names(mean_r) <- names(nhost) <- names(nsrc) <- ids

  cat(sprintf("[%s] table %d ASVs x %d samples; %d host species; %d sources\n",
              reg, nrow(mat), ncol(mat), data.table::uniqueN(sp), data.table::uniqueN(so)))

  tax <- data.table::fread(file.path(BASE_DIR, g$tax), sep = "\t")
  stopifnot("taxid" %in% names(tax))
  tax_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                        names(tax))

  fa <- read_fasta(file.path(BASE_DIR, g$fa))
  cat(sprintf("[%s] filtered FASTA: %d sequences\n", reg, length(fa)))

  for (fl in names(FILTERS)) {
    f <- FILTERS[[fl]]
    keep_all <- ids[prev[ids] >= f$prev & nhost[ids] >= f$host & nsrc[ids] >= f$src]
    in_fa <- keep_all %in% names(fa)
    keep <- keep_all[in_fa]
    n_missing <- sum(!in_fa)
    if (n_missing > 0L) {
      cat(sprintf("  [%s/%s] WARNING: %d candidate ASVs missing from FASTA; dropped\n",
                  reg, fl, n_missing))
    }
    if (length(keep) < 10L) stop(sprintf("[%s/%s] too few candidates", reg, fl))

    tax_idx <- match(keep, tax$taxid)
    dt <- data.table::data.table(
      ASV_ID               = keep,
      region               = reg,
      filter               = fl,
      prevalence_samples   = as.integer(prev[keep]),
      prevalence_frac      = as.numeric(prev[keep] / ncol(mat)),
      n_host_species       = as.integer(nhost[keep]),
      n_sources            = as.integer(nsrc[keep]),
      total_css            = as.numeric(total[keep]),
      mean_rel_abund       = as.numeric(mean_r[keep]),
      sequence_length      = nchar(fa[keep])
    )
    for (cc in tax_cols) {
      v <- rep(NA_character_, length(keep))
      ok <- !is.na(tax_idx)
      v[ok] <- as.character(tax[[cc]][tax_idx[ok]])
      dt[[cc]] <- v
    }
    data.table::setorder(dt, -prevalence_samples, ASV_ID)

    cand_file <- file.path(OUT_ROOT, "Temp_data/candidates",
                           sprintf("%s_%s_candidates.tsv", reg, fl))
    fa_file <- file.path(OUT_ROOT, "Temp_data/tree_input",
                         sprintf("%s_%s.fa", reg, fl))
    data.table::fwrite(dt, cand_file, sep = "\t")
    writeLines(as.vector(rbind(paste0(">", dt$ASV_ID), fa[dt$ASV_ID])), fa_file)

    cat(sprintf("  [%s/%s] %d candidates (prev>=%d, hosts>=%d, sources>=%d); FASTA -> %s\n",
                reg, fl, nrow(dt), f$prev, f$host, f$src, basename(fa_file)))

    summary_rows[[length(summary_rows) + 1L]] <- data.table::data.table(
      region = reg, filter = fl,
      n_candidates = nrow(dt), n_missing_fasta = n_missing,
      n_samples = ncol(mat), n_host_species_total = data.table::uniqueN(sp),
      n_sources_total = data.table::uniqueN(so),
      prevalence_min = min(dt$prevalence_samples),
      prevalence_median = as.numeric(stats::median(dt$prevalence_samples)),
      prevalence_max = max(dt$prevalence_samples),
      n_host_min = min(dt$n_host_species), n_sources_min = min(dt$n_sources)
    )
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table::data.table(
      region = reg, filter = fl,
      fasta = fa_file,
      expected_tree = file.path(OUT_ROOT, "Temp_data/trees",
                                sprintf("%s_%s.nwk", reg, fl)),
      n_sequences = nrow(dt)
    )
  }
  rm(mat, rel, fa); invisible(gc())
}

summ <- data.table::rbindlist(summary_rows)
man  <- data.table::rbindlist(manifest_rows)
data.table::fwrite(summ, file.path(OUT_ROOT, "result_table", "candidate_filter_summary.tsv"), sep = "\t")
data.table::fwrite(man,  file.path(OUT_ROOT, "result_table", "tree_build_manifest.tsv"), sep = "\t")

cat("\n=== candidate summary ===\n")
print(summ)
cat("\nDONE 00_prepare_candidates\n")
