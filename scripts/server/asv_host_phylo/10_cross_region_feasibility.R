#!/usr/bin/env Rscript
# Feasibility check: can mapped reference bins in the cross-region merged
# 16S / ITS abundance tables support a cross-region ASV-distribution analysis?
#
# "mapped" = row ID present in Merged_{16S_Silva,ITS_Unite}_filtered.fa
# (vsearch --usearch_global at 97% identity; multiple ASVs from any region
#  collapse onto the same SILVA / UNITE reference ID).
#
# Outputs: result_table/cross_region_feasibility/
#   cross_region_sample_summary.tsv
#   cross_region_region_span.tsv
#   cross_region_mapped_overlap.tsv
#   cross_region_mapped_shared_rows.tsv

suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

BASE <- "<SERVER_BASE>"
setwd(BASE)
OUT <- file.path("result", "ASV_host_phylo", "result_table", "cross_region_feasibility")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

norm_id <- function(x) {
  x <- gsub("-", ".", x, fixed = TRUE)
  ifelse(grepl("^[0-9]", x), paste0("X", x), x)
}

read_fa_ids <- function(f) {
  x <- readLines(f, warn = FALSE)
  sub("^>([^ ]+).*$", "\\1", x[startsWith(x, ">")])
}

summarise_dataset <- function(ds) {
  cat("\n=== ", ds, " ===\n", sep = "")
  tab_file <- file.path("filtered_merge", ds,
                        sprintf("Merged_%s_table_filtered.txt", ds))
  ref_file <- file.path("filtered_merge", ds,
                        sprintf("Merged_%s_%s_filtered.fa", ds,
                                if (ds == "16S") "Silva" else "Unite"))
  tax_file <- file.path("filtered_merge", ds,
                        sprintf("Merged_%s_taxonomy_filtered.txt", ds))
  meta_file <- file.path("standard_merge", sprintf("metadata_%s.txt", ds))

  ref_ids <- read_fa_ids(ref_file)
  tax <- unique(fread(tax_file), by = "taxid")
  meta <- fread(meta_file)
  meta[, SampleID := norm_id(SampleID)]

  cat("reading table ...\n")
  tab <- fread(tab_file, sep = "\t", header = TRUE, check.names = FALSE)
  setnames(tab, 1L, "ASV_ID")
  sample_cols <- names(tab)[-1L]
  norm_samples <- norm_id(sample_cols)
  idx <- match(norm_samples, meta$SampleID)
  if (anyNA(idx)) {
    stop(sprintf("%s: %d table samples not matched to metadata",
                 ds, sum(is.na(idx))))
  }
  region <- meta$Region[idx]
  species <- meta$species[idx]
  source <- meta$source[idx]
  regs <- sort(unique(region))

  sample_summary <- data.table(
    dataset = ds,
    region = regs,
    n_samples = as.integer(table(region)[regs])
  )
  n_sp <- vapply(regs, function(r) uniqueN(species[region == r]), integer(1))
  n_src <- vapply(regs, function(r) uniqueN(source[region == r]), integer(1))
  sample_summary[, n_species := n_sp]
  sample_summary[, n_sources := n_src]

  # presence / abundance summaries
  cat("computing region spans ...\n")
  region_pos <- sapply(regs, function(r) {
    cols <- sample_cols[region == r]
    sub <- as.matrix(tab[, ..cols])
    rowSums(sub > 0) > 0
  })
  if (!is.matrix(region_pos)) region_pos <- matrix(region_pos, ncol = 1)
  colnames(region_pos) <- regs
  n_regions <- rowSums(region_pos)

  total_prev <- rowSums(tab[, ..sample_cols] > 0)
  total_abund <- rowSums(tab[, ..sample_cols])
  is_mapped <- tab$ASV_ID %in% ref_ids

  span <- data.table(
    dataset = ds,
    ASV_ID = tab$ASV_ID,
    mapped = is_mapped,
    n_regions = n_regions,
    total_prevalence = total_prev,
    total_abundance = total_abund
  )
  for (r in regs) {
    cols <- sample_cols[region == r]
    sub <- as.matrix(tab[, ..cols])
    span[, (paste0("prev_", r)) := rowSums(sub > 0)]
    span[, (paste0("abund_", r)) := rowSums(sub)]
  }

  # host-species / source breadth only for rows detected in >=2 regions
  shared <- which(n_regions >= 2)
  cat("rows detected in >=2 regions:", length(shared), "\n")
  if (length(shared)) {
    sub_tab <- tab[shared, ..sample_cols]
    pres <- sub_tab > 0
    sp_levels <- sort(unique(species))
    sp_counts <- matrix(0L, nrow = nrow(pres), ncol = length(sp_levels))
    for (j in seq_along(sp_levels)) {
      cols <- which(species == sp_levels[j])
      sp_counts[, j] <- rowSums(pres[, cols, drop = FALSE])
    }
    n_hosts <- rowSums(sp_counts > 0)
    src_levels <- sort(unique(source))
    n_sources <- integer(nrow(pres))
    for (s in src_levels) {
      cols <- which(source == s)
      n_sources <- n_sources + (rowSums(pres[, cols, drop = FALSE]) > 0)
    }
    span[shared, n_hosts := n_hosts]
    span[shared, n_sources := n_sources]
  }

  # attach taxonomy (only for rows that exist in taxonomy)
  span <- merge(span, tax, by.x = "ASV_ID", by.y = "taxid",
                all.x = TRUE, sort = FALSE)

  # --- output 1: region span distribution ---
  span_dist <- span[, .(n_rows = .N), by = .(dataset, mapped, n_regions)]
  span_dist <- span_dist[order(mapped, n_regions)]
  fwrite(span_dist, file.path(OUT, sprintf("cross_region_region_span_%s.tsv", ds)),
         sep = "\t")

  # --- output 2: mapped overlap among regions ---
  mapped_span <- span[mapped == TRUE]
  overlap <- data.table(
    dataset = ds,
    n_mapped = nrow(mapped_span),
    n_mapped_shared2 = sum(mapped_span$n_regions >= 2),
    pct_mapped_shared2 = 100 * mean(mapped_span$n_regions >= 2),
    n_mapped_all_regions = sum(mapped_span$n_regions == length(regs)),
    abundance_share_shared2 = 100 * sum(mapped_span$total_abundance[
      mapped_span$n_regions >= 2]) / sum(mapped_span$total_abundance),
    n_unmapped = sum(!span$mapped),
    n_unmapped_shared2 = sum(!span$mapped & span$n_regions >= 2)
  )
  fwrite(overlap, file.path(OUT, sprintf("cross_region_mapped_overlap_%s.tsv", ds)),
         sep = "\t")

  # --- output 3: shared-row detail ---
  keep_cols <- c("dataset", "ASV_ID", "mapped", "n_regions",
                 "total_prevalence", "total_abundance", "n_hosts", "n_sources",
                 paste0("prev_", regs), paste0("abund_", regs),
                 "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                 "Species")
  keep_cols <- intersect(keep_cols, names(span))
  shared_detail <- span[shared, ..keep_cols]
  shared_detail[, abund_share_pooled := 100 * total_abundance /
                  sum(span$total_abundance)]
  shared_detail[, filter_pass := total_prevalence >= 20 & n_hosts >= 5 &
                  n_sources >= 3]
  fwrite(shared_detail,
         file.path(OUT, sprintf("cross_region_mapped_shared_rows_%s.tsv", ds)),
         sep = "\t")

  # region-pair overlap among mapped bins
  pair_rows <- list()
  if (length(regs) >= 2) {
    for (i in seq_len(length(regs) - 1)) {
      for (j in (i + 1):length(regs)) {
        a <- regs[i]; b <- regs[j]
        both <- mapped_span[[paste0("prev_", a)]] > 0 &
                mapped_span[[paste0("prev_", b)]] > 0
        pair_rows[[length(pair_rows) + 1]] <- data.table(
          dataset = ds, region_a = a, region_b = b,
          n_shared = sum(both),
          n_shared_pass = sum(both & mapped_span$total_prevalence >= 20 &
                              mapped_span$n_hosts >= 5 &
                              mapped_span$n_sources >= 3),
          n_abund_shared = sum(mapped_span$total_abundance[both])
        )
      }
    }
  }
  pair_tab <- rbindlist(pair_rows)

  # host-species overlap across regions
  host_rows <- list()
  if (length(regs) >= 2) {
    sp_sets <- lapply(regs, function(r) unique(species[region == r]))
    names(sp_sets) <- regs
    for (i in seq_len(length(regs) - 1)) {
      for (j in (i + 1):length(regs)) {
        a <- regs[i]; b <- regs[j]
        host_rows[[length(host_rows) + 1]] <- data.table(
          dataset = ds, region_a = a, region_b = b,
          n_species_a = length(sp_sets[[a]]),
          n_species_b = length(sp_sets[[b]]),
          n_shared = length(intersect(sp_sets[[a]], sp_sets[[b]])),
          jaccard = length(intersect(sp_sets[[a]], sp_sets[[b]])) /
                    length(union(sp_sets[[a]], sp_sets[[b]]))
        )
      }
    }
  }
  host_tab <- rbindlist(host_rows)

  # within-region candidate counts for a region-stratified meta-analysis
  reg_cand <- list()
  for (r in regs) {
    cols_r <- sample_cols[region == r]
    sp_r <- species[region == r]
    src_r <- source[region == r]
    prev_r <- span[[paste0("prev_", r)]]
    sel <- which(is_mapped & prev_r >= 20)
    sp_lv <- unique(sp_r)
    src_lv <- unique(src_r)
    if (length(sel)) {
      sub <- as.matrix(tab[sel, ..cols_r])
      pres <- sub > 0
      nh <- integer(length(sel))
      for (s in sp_lv) {
        cc <- which(sp_r == s)
        nh <- nh + (rowSums(pres[, cc, drop = FALSE]) > 0)
      }
      ns <- integer(length(sel))
      for (s in src_lv) {
        cc <- which(src_r == s)
        ns <- ns + (rowSums(pres[, cc, drop = FALSE]) > 0)
      }
      host_thr <- min(5L, length(sp_lv))
      src_thr <- min(3L, length(src_lv))
      reg_cand[[r]] <- data.table(
        dataset = ds, region = r, n_samples = length(cols_r),
        n_species = length(sp_lv), n_sources = length(src_lv),
        n_mapped_prev20 = length(sel),
        host_thr = host_thr, source_thr = src_thr,
        n_pass_adaptive = sum(nh >= host_thr & ns >= src_thr),
        n_pass_5_3 = sum(nh >= 5L & ns >= 3L)
      )
    } else {
      reg_cand[[r]] <- data.table(
        dataset = ds, region = r, n_samples = length(cols_r),
        n_species = length(sp_lv), n_sources = length(src_lv),
        n_mapped_prev20 = 0L, host_thr = min(5L, length(sp_lv)),
        source_thr = min(3L, length(src_lv)),
        n_pass_adaptive = 0L, n_pass_5_3 = 0L
      )
    }
  }
  reg_cand <- rbindlist(reg_cand)

  list(sample_summary = sample_summary, overlap = overlap,
       span_dist = span_dist, span = span,
       pair_tab = pair_tab, host_tab = host_tab,
       reg_cand = reg_cand)
}

res16 <- summarise_dataset("16S")
resit <- summarise_dataset("ITS")

sample_summary <- rbind(res16$sample_summary, resit$sample_summary)
overlap <- rbind(res16$overlap, resit$overlap)
pair_tab <- rbind(res16$pair_tab, resit$pair_tab)
host_tab <- rbind(res16$host_tab, resit$host_tab)
reg_cand <- rbind(res16$reg_cand, resit$reg_cand)
fwrite(sample_summary, file.path(OUT, "cross_region_sample_summary.tsv"),
       sep = "\t")
fwrite(overlap, file.path(OUT, "cross_region_mapped_overlap.tsv"), sep = "\t")
fwrite(pair_tab, file.path(OUT, "cross_region_pair_overlap.tsv"), sep = "\t")
fwrite(host_tab, file.path(OUT, "cross_region_host_overlap.tsv"), sep = "\t")
fwrite(reg_cand, file.path(OUT, "cross_region_region_candidates.tsv"), sep = "\t")

cat("\n=== sample summary ===\n")
print(sample_summary)
cat("\n=== mapped overlap ===\n")
print(overlap)
cat("\n=== within-region candidate counts (mapped, prev>=20) ===\n")
print(reg_cand)
cat("\nDone. Outputs in ", normalizePath(OUT), "\n", sep = "")
