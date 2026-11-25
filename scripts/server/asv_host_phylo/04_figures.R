#!/usr/bin/env Rscript
# =============================================================================
# 04_figures.R — figures for the ASV-host-phylogeny analysis
#
# All figures are built from result_table/*.tsv; no conclusions are read from
# the rendered images (per the project's paper-writing rule).
#
# Outputs (result_figure):
#   FigA_candidate_filtering.{pdf,jpg}
#   FigB_mantel_scatter.{pdf,jpg}
#   FigC_correlogram.{pdf,jpg}
#   FigD_sensitivity_forest.{pdf,jpg}
#   FigE_close_pair_contrast.{pdf,jpg}
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/asv_host_phylo_utils.R"))

FIG <- file.path(OUT_ROOT, "result_figure")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(FIG, paste0(stem, ".jpg")), p, width = w, height = h,
         dpi = 300, bg = "white", limitsize = FALSE)
  message(sprintf("[figure] %s (%.1f x %.1f in)", stem, w, h))
}

region_lab <- c(v4 = "V4", v5_v7 = "V5-V7", its1 = "ITS1")
metric_lab <- c(J = "Jaccard", B = "Bray-Curtis", H = "Jaccard", HB = "Bray-Curtis")
metric_levels <- c("Jaccard", "Bray-Curtis")
filter_lab <- c(main20_5_3 = "Main 20/5/3",
                loose10_3_2 = "Loose 10/3/2",
                strict30_5_3 = "Strict 30/5/3")

# ---------------------------------------------------------------------------
# Fig A: candidate filtering
# ---------------------------------------------------------------------------
summ <- fread(file.path(OUT_ROOT, "result_table", "candidate_filter_summary.tsv"))
summ[, `:=`(region_l = factor(region_lab[region], levels = region_lab),
            filter_l = factor(filter_lab[filter], levels = filter_lab))]

pA1 <- ggplot(summ, aes(x = region_l, y = n_candidates, fill = filter_l)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n_candidates), position = position_dodge(width = 0.75),
            vjust = -0.3, size = 2.7) +
  scale_fill_manual(values = c("Main 20/5/3" = "#3D5A80",
                               "Loose 10/3/2" = "#98C1D9",
                               "Strict 30/5/3" = "#EE6C4D"), name = NULL) +
  labs(title = "Candidate ASV sets", x = NULL, y = "Number of ASVs") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

cand <- data.table::rbindlist(lapply(c("v4", "v5_v7", "its1"), function(r) {
  d <- fread(cand_file(r, "loose10_3_2")); d[, region := r]; d
}))
cand[, `:=`(region_l = factor(region_lab[region], levels = region_lab),
            log_prev = log10(prevalence_samples))]
pA2 <- ggplot(cand, aes(x = log_prev, y = n_host_species, colour = n_sources)) +
  geom_point(size = 0.55, alpha = 0.45) +
  geom_vline(xintercept = log10(20), linetype = 2, colour = "grey35") +
  geom_hline(yintercept = 5, linetype = 2, colour = "grey35") +
  scale_colour_viridis_c(trans = "log1p", name = "Sources") +
  facet_wrap(~ region_l, nrow = 1) +
  labs(title = "Loose-filter candidate space (filter thresholds dashed)",
       x = "log10(prevalence, samples)", y = "Host species breadth") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

pA <- pA1 / pA2 + plot_layout(heights = c(1, 1.6))
save_fig(pA, "FigA_candidate_filtering", 10, 7)

# ---------------------------------------------------------------------------
# Fig B: Mantel scatter (subsampled pairs) + MRM phylo p
# ---------------------------------------------------------------------------
mt <- fread(file.path(OUT_ROOT, "result_table", "mantel_results.tsv"))
mrm <- fread(file.path(OUT_ROOT, "result_table", "mrm_results.tsv"))
mrm_phylo <- mrm[term == "dp", .(region, metric, mrm_coef = coef, mrm_p = p)]
mt <- merge(mt, mrm_phylo, by = c("region", "metric"), all.x = TRUE)

scatter_rows <- list()
for (reg in c("v4", "v5_v7", "its1")) {
  p <- readRDS(profile_file(reg, "main20_5_3"))
  dp <- as.matrix(p$dp); lt <- lower.tri(dp)
  idx <- which(lt, arr.ind = TRUE)
  set.seed(20261301L + match(reg, c("v4", "v5_v7", "its1")))
  sel <- sample.int(nrow(idx), min(30000L, nrow(idx)))
  for (metric in c("J", "B")) {
    dm <- as.matrix(p[[if (metric == "J") "dj" else "db"]])
    scatter_rows[[length(scatter_rows) + 1L]] <- data.table(
      region = reg, metric = metric,
      phylo = dp[idx[sel, , drop = FALSE]],
      dist = dm[idx[sel, , drop = FALSE]]
    )
  }
  rm(p); invisible(gc())
}
sc <- rbindlist(scatter_rows)
sc <- sc[phylo > 0]
sc[, `:=`(region_l = factor(region_lab[region], levels = region_lab),
          metric_l = factor(metric_lab[metric], levels = metric_levels))]
ann_pos <- sc[, .(x = min(phylo, na.rm = TRUE),
                   y = max(dist, na.rm = TRUE)),
               by = .(region_l, metric_l)]
ann <- merge(
  ann_pos,
  mt[, .(region_l = factor(region_lab[region], levels = region_lab),
         metric_l = factor(metric_lab[metric], levels = metric_levels),
         label = sprintf("r = %.3f\nP = %.2g\nSES = %.1f\nMRM P = %.2g",
                         r, p_standard, SES_standard, mrm_p))],
  by = c("region_l", "metric_l"))

pB <- ggplot(sc, aes(x = phylo, y = dist)) +
  geom_hex(bins = 55) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#C0392B", linewidth = 0.7) +
  geom_text(data = ann, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 2.6) +
  scale_fill_viridis_c(trans = "log10", name = "Pairs") +
  scale_x_log10() +
  facet_grid(metric_l ~ region_l, scales = "free") +
  labs(x = "ASV patristic distance (log scale)",
       y = "ASV distribution dissimilarity") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")
save_fig(pB, "FigB_mantel_scatter", 12, 7)

# ---------------------------------------------------------------------------
# Fig C: Mantel correlogram
# ---------------------------------------------------------------------------
cg <- fread(file.path(OUT_ROOT, "result_table", "correlogram.tsv"))
cg[, `:=`(region_l = factor(region_lab[region], levels = region_lab),
          metric_l = factor(metric_lab[metric], levels = metric_levels))]
cg_plot <- cg[phylo_mid > 0]
cg1 <- merge(
  cg_plot[, .(x = min(phylo_mid, na.rm = TRUE),
              y = max(null_q975, obs_mean, na.rm = TRUE)),
          by = .(region_l, metric_l)],
  cg[bin == 1L, .(region_l, metric_l, SES, p_lower,
                  label = sprintf("bin1 SES = %.1f\nP(lower) = %.3f", SES, p_lower))],
  by = c("region_l", "metric_l"))
pC <- ggplot(cg_plot, aes(x = phylo_mid)) +
  geom_ribbon(data = cg_plot, aes(ymin = null_q025, ymax = null_q975),
              alpha = 0.22, fill = "grey55") +
  geom_line(data = cg_plot, aes(y = null_mean), linetype = 2, colour = "grey30") +
  geom_line(data = cg_plot, aes(y = obs_mean), colour = "#D55E00", linewidth = 0.7) +
  geom_point(data = cg_plot, aes(y = obs_mean), colour = "#D55E00", size = 1.3) +
  geom_text(data = cg1, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 2.4) +
  scale_x_log10() +
  facet_grid(metric_l ~ region_l, scales = "free") +
  labs(x = "Phylogenetic distance bin midpoint (log scale)",
       y = "Mean distribution dissimilarity") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
save_fig(pC, "FigC_correlogram", 12, 7)

# ---------------------------------------------------------------------------
# Fig D: sensitivity forest
# ---------------------------------------------------------------------------
sens <- fread(file.path(OUT_ROOT, "result_table", "sensitivity_results.tsv"))
main <- copy(mt)
main[, `:=`(analysis = "main", subset = filter, n_units = NA_integer_,
            p = p_standard)]
main <- main[, .(region, analysis, subset, metric, n_asv, n_units, r, p)]
forest <- rbindlist(list(main, sens[, .(region, analysis, subset, metric,
                                        n_asv, n_units, r, p)]), fill = TRUE)
forest[, analysis_l := ifelse(
  analysis == "main", "Main 20/5/3",
  ifelse(analysis == "alt_filter" & subset == "loose10_3_2", "Loose 10/3/2",
    ifelse(analysis == "alt_filter" & subset == "strict30_5_3", "Strict 30/5/3",
      ifelse(analysis == "within_largest_source", "Within largest source",
        ifelse(analysis == "host_level", "Host-species level", analysis)))))]
forest[, `:=`(region_l = factor(region_lab[region], levels = region_lab),
              metric_l = factor(metric_lab[metric], levels = metric_levels),
              analysis_l = factor(analysis_l,
                levels = c("Main 20/5/3", "Loose 10/3/2", "Strict 30/5/3",
                           "Within largest source", "Host-species level")))]

pD <- ggplot(forest, aes(x = analysis_l, y = r, colour = region_l)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  geom_text(aes(label = ifelse(p < 0.001, "***",
                               ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))),
            position = position_dodge(width = 0.55), vjust = -0.8, size = 2.6,
            show.legend = FALSE) +
  facet_wrap(~ metric_l, scales = "free_y") +
  scale_colour_manual(values = c("V4" = "#4E79A7", "V5-V7" = "#F28E2B",
                                 "ITS1" = "#59A14F"), name = NULL) +
  labs(x = NULL, y = "Mantel r") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank(), legend.position = "top")
save_fig(pD, "FigD_sensitivity_forest", 10, 5)

# ---------------------------------------------------------------------------
# Fig E: close-vs-far contrast
# ---------------------------------------------------------------------------
cp_file <- file.path(OUT_ROOT, "result_table", "close_pair_contrast.tsv")
if (file.exists(cp_file)) {
  cp <- fread(cp_file)
  cp[, `:=`(region_l = factor(region_lab[region], levels = region_lab),
            metric_l = factor(metric_lab[metric], levels = metric_levels))]
  pE <- ggplot(cp, aes(x = region_l, y = contrast, colour = metric_l)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(ymin = contrast - 1.96 * null_sd,
                        ymax = contrast + 1.96 * null_sd),
                    position = position_dodge(width = 0.4), size = 0.45) +
    facet_wrap(~ metric_l, scales = "free_y") +
    scale_colour_manual(values = c("Jaccard" = "#0072B2",
                                   "Bray-Curtis" = "#D55E00"), name = NULL) +
    labs(x = NULL,
         y = "Mean dissimilarity (far 50% - close 10%)\npositive = close relatives more similar") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  save_fig(pE, "FigE_close_pair_contrast", 8, 4)
}

cat("DONE 04_figures\n")
