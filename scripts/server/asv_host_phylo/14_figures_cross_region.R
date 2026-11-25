#!/usr/bin/env Rscript
# Figures for the cross-region mapped-bin analysis (from result TSVs only).

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
BASE <- Sys.getenv("BASE_DIR", "<SERVER_BASE>")
TAB <- file.path(BASE, "result/ASV_host_phylo/result_table/cross_region")
FIG <- file.path(BASE, "result/ASV_host_phylo/result_figure/cross_region")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

read_if <- function(f) if (file.exists(f)) fread(f) else NULL

d1 <- rbindlist(list(
  read_if(file.path(TAB, "track1_16s_mantel.tsv")),
  read_if(file.path(TAB, "track1_its_mantel.tsv"))
), fill = TRUE)
if (!is.null(d1) && nrow(d1)) {
  d1[, metric := factor(metric, levels = c("J", "B"),
                        labels = c("Jaccard", "Bray-Curtis"))]
  p1 <- ggplot(d1, aes(x = r, y = reorder(region, r), colour = metric)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                    position = position_dodge(width = 0.45),
                    linewidth = 0.4, size = 0.25) +
    facet_wrap(~dataset, scales = "free_y") +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = "Mantel r (95% bootstrap CI)", y = NULL,
         colour = "Distribution distance",
         title = "Track 1: within-region mapped-bin phylogenetic signal") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(FIG, "FigCrossA_track1_forest.pdf"), p1,
         width = 8, height = 5, device = "pdf")
  ggsave(file.path(FIG, "FigCrossA_track1_forest.jpg"), p1,
         width = 8, height = 5, dpi = 300)
}

d2 <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_mantel.tsv")),
  read_if(file.path(TAB, "track2_its_mantel.tsv"))
), fill = TRUE)
if (!is.null(d2) && nrow(d2)) {
  d2[, pair := paste(region_a, region_b, sep = " - ")]
  d2[, metric := factor(metric, levels = c("J", "B"),
                        labels = c("Jaccard", "Bray-Curtis"))]
  p2 <- ggplot(d2, aes(x = r, y = reorder(pair, r), colour = pair_type)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                    linewidth = 0.4, size = 0.25) +
    facet_grid(dataset ~ metric, scales = "free_y") +
    scale_colour_brewer(palette = "Set2") +
    labs(x = "Mantel r (95% CI)", y = NULL, colour = "Region pair",
         title = "Track 2: cross-region shared-bin phylogenetic signal") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(FIG, "FigCrossB_track2_forest.pdf"), p2,
         width = 10, height = 5.5, device = "pdf")
  ggsave(file.path(FIG, "FigCrossB_track2_forest.jpg"), p2,
         width = 10, height = 5.5, dpi = 300)
}

d3a <- rbindlist(list(
  read_if(file.path(TAB, "track1_16s_close_far.tsv")),
  read_if(file.path(TAB, "track1_its_close_far.tsv"))
), fill = TRUE)
d3b <- rbindlist(list(
  read_if(file.path(TAB, "track2_16s_close_far.tsv")),
  read_if(file.path(TAB, "track2_its_close_far.tsv"))
), fill = TRUE)
if (!is.null(d3a) && nrow(d3a)) {
  d3a[, unit := region]
  d3a[, analysis := "Track 1: within-region"]
}
if (!is.null(d3b) && nrow(d3b)) {
  d3b[, unit := paste(region_a, region_b, sep = " - ")]
  d3b[, analysis := "Track 2: cross-region pair"]
}
d3 <- rbindlist(list(d3a, d3b), fill = TRUE)
if (!is.null(d3) && nrow(d3)) {
  d3[, metric := factor(metric, levels = c("J", "B"),
                        labels = c("Jaccard", "Bray-Curtis"))]
  p3 <- ggplot(d3, aes(x = contrast, y = reorder(unit, contrast),
                       colour = metric)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_point(size = 1.8) +
    facet_wrap(~dataset + analysis, scales = "free_y") +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = "Far 50% mean distance - closest 10% mean distance",
         y = NULL, colour = "Distribution distance",
         title = "Short-range phylogenetic signal") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(FIG, "FigCrossC_close_far.pdf"), p3,
         width = 9, height = 5, device = "pdf")
  ggsave(file.path(FIG, "FigCrossC_close_far.jpg"), p3,
         width = 9, height = 5, dpi = 300)
}
cat("Figures written to", FIG, "\n")
