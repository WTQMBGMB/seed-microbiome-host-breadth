#!/usr/bin/env Rscript
# Prepare the ITS reference tree and a duplicate-ID -> tree-tip map.
#
# FastTree collapses identical trimmed sequences. We map every mapped UNITE
# reference ID to a tree tip carrying the identical trimmed sequence; members
# of the same group then have patristic distance 0, preserving the distance
# geometry for all 4,918 mapped bins.

suppressPackageStartupMessages({ library(data.table); library(ape) })
BASE <- Sys.getenv("BASE_DIR", "<SERVER_BASE>")
TREE_DIR <- file.path(BASE, "result/ASV_host_phylo/Temp_data/trees")
dir.create(TREE_DIR, recursive = TRUE, showWarnings = FALSE)

TREE_IN  <- file.path(TREE_DIR, "its_ref_unite_raw.nwk")
FA_IN    <- file.path(TREE_DIR, "its_ref_gt80.fa")
TREE_OUT <- file.path(TREE_DIR, "its_ref_unite.nwk")
MAP_OUT  <- file.path(TREE_DIR, "its_ref_dup_map.tsv")

stopifnot(file.exists(TREE_IN), file.exists(FA_IN))
tr <- read.tree(TREE_IN)
seqs <- list()
id <- NULL
for (line in readLines(FA_IN, warn = FALSE)) {
  if (startsWith(line, ">")) {
    id <- sub("^>([^ ]+).*$", "\\1", line)
    seqs[[id]] <- character(0)
  } else {
    seqs[[id]] <- c(seqs[[id]], line)
  }
}
seq_str <- vapply(seqs, paste0, character(1), collapse = "")
grp <- split(names(seq_str), seq_str)

rows <- list()
for (g in grp) {
  tip <- intersect(g, tr$tip.label)
  rep <- if (length(tip)) tip[1] else NA_character_
  rows[[length(rows) + 1L]] <- data.table(
    id = g, rep = rep, group_size = length(g),
    rep_in_tree = !is.na(rep)
  )
}
map <- rbindlist(rows)

all_ids <- sub("^>([^ ]+).*$", "\\1",
               grep("^>", readLines(file.path(BASE,
                 "filtered_merge/ITS/Merged_ITS_Unite_filtered.fa"),
                 warn = FALSE), value = TRUE))
missing <- setdiff(all_ids, map$id)
if (length(missing)) {
  cat("IDs absent from trimmed FASTA:", length(missing), "\n")
}
unmapped <- map[is.na(rep)]
cat(sprintf("tree tips=%d; mapped IDs=%d; unique groups=%d; unmapped groups=%d (%d IDs)\n",
            length(tr$tip.label), nrow(map), uniqueN(map$rep[!is.na(map$rep)]),
            nrow(unmapped), sum(unmapped$group_size)))

# Keep the tree as-is; phylo_dist() expands duplicate IDs to representatives.
file.copy(TREE_IN, TREE_OUT, overwrite = TRUE)
fwrite(map[, .(id, rep, group_size, rep_in_tree)], MAP_OUT, sep = "\t")
cat("Wrote", TREE_OUT, "and", MAP_OUT, "\n")
