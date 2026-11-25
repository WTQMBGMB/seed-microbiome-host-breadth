#!/usr/bin/env Rscript
# Single main figure + single main table for the paper:
# only 16S V4 and V5-V7 (within-region) and V4-V5-V7 shared bins.
#
# Outputs:
#   paper_main/Fig1_cross_region.{pdf,jpg}
#   paper_main/Table1_cross_region.tsv
#   paper_main/Table1_cross_region.md
#   paper_main/source_data_Fig1.tsv

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
BASE <- Sys.getenv("BASE_DIR", "<SERVER_BASE>")
TAB <- file.path(BASE, "result/ASV_host_phylo/result_table/cross_region")
OUT <- file.path(BASE, "result/ASV_host_phylo/paper_main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

t1m <- fread(file.path(TAB, "track1_16s_mantel.tsv"))
t1r <- fread(file.path(TAB, "track1_16s_mrm.tsv"))
t1c <- fread(file.path(TAB, "track1_16s_close_far.tsv"))
t2m <- fread(file.path(TAB, "track2_16s_mantel.tsv"))
t2r <- fread(file.path(TAB, "track2_16s_mrm.tsv"))
t2c <- fread(file.path(TAB, "track2_16s_close_far.tsv"))
t2s <- fread(file.path(TAB, "track2_16s_sensitivity.tsv"))

t1m <- t1m[region %in% c("V4", "V5_V7")]
t1r <- t1r[region %in% c("V4", "V5_V7") & term == "dp" & model == "primary"]
t1c <- t1c[region %in% c("V4", "V5_V7")]
t2m <- t2m[region_a == "V4" & region_b == "V5_V7"]
t2r <- t2r[region_a == "V4" & region_b == "V5_V7" &
           term == "dp" & model == "primary"]
t2c <- t2c[region_a == "V4" & region_b == "V5_V7"]
t2s <- t2s[region_a == "V4" & region_b == "V5_V7"]

null_file <- file.path(OUT, "null_model_results.tsv")
null_draw_file <- file.path(OUT, "null_draws.rds")
null_res <- if (file.exists(null_file)) fread(null_file) else NULL
null_draws <- if (file.exists(null_draw_file)) readRDS(null_draw_file) else NULL

unit_levels <- c("V4 (within)", "V5-V7 (within)", "V4-V5-V7 (shared)")
unit_lab <- c(V4 = "V4 (within)", V5_V7 = "V5-V7 (within)")

# ---- long effect table -----------------------------------------------------
eff1 <- t1m[, .(unit = unit_lab[region], n_bins = n_bins,
                metric = metric, r = r, CI_low = CI_low,
                CI_high = CI_high, p = p_standard,
                p_strat = p_stratified)]
eff2 <- t2m[, .(unit = "V4-V5-V7 (shared)", n_bins = n_bins,
                metric = metric, r = r, CI_low = CI_low,
                CI_high = CI_high, p = p_standard,
                p_strat = p_stratified)]
eff <- rbind(eff1, eff2)
eff[, unit := factor(unit, levels = unit_levels)]
eff[, metric := factor(metric, levels = c("J", "B"))]

short1 <- t1c[, .(unit = unit_lab[region], metric = metric,
                  contrast = contrast, p = p, SES = SES)]
short2 <- t2c[, .(unit = "V4-V5-V7 (shared)", metric = metric,
                  contrast = contrast, p = p, SES = SES)]
short <- rbind(short1, short2)
short[, unit := factor(unit, levels = unit_levels)]
short[, metric := factor(metric, levels = c("J", "B"))]

# cross-region robustness: primary + sensitivity scenarios
prim <- t2m[, .(scenario = "primary", metric = metric, r = r,
                p = p_standard, n_bins = n_bins)]
sens <- t2s[, .(scenario = analysis, metric = metric, r = r,
                p = p, n_bins = n_bins)]
sens <- rbind(prim, sens)
sens[, scenario := factor(
  scenario,
  levels = c("primary", "within_min10", "pair_min", "pair_max"),
  labels = c("Primary", "Within-region n >= 10", "Pair min", "Pair max"))]
sens[, metric := factor(metric, levels = c("J", "B"))]

if (!is.null(null_res)) {
  null_res[, unit := factor(unit, levels = unit_levels)]
  null_res[, metric := factor(metric, levels = c("J", "B"))]
}
if (!is.null(null_draws)) {
  null_draws[, unit := factor(unit, levels = unit_levels)]
  null_draws[, metric := factor(metric, levels = c("J", "B"))]
}

# ---- panel A ---------------------------------------------------------------
star <- function(p) ifelse(
  is.na(p), "",
  ifelse(p <= 0.001, "***",
         ifelse(p <= 0.01, "**",
                ifelse(p <= 0.05, "*", "ns"))))

if (!is.null(null_res) && !is.null(null_draws)) {
  null_res[, sig := star(p_source_stratified)]
  pA <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    geom_jitter(data = null_draws,
                aes(x = unit, y = r, colour = metric),
                position = position_jitterdodge(jitter.width = 0.16,
                                                dodge.width = 0.60),
                size = 0.5, alpha = 0.25, stroke = 0) +
    geom_pointrange(data = eff,
                    aes(x = unit, y = r, ymin = CI_low, ymax = CI_high,
                        colour = metric),
                    position = position_dodge(width = 0.60),
                    linewidth = 0.45, size = 0.25) +
    geom_point(data = eff,
               aes(x = unit, y = r, colour = metric),
               position = position_dodge(width = 0.60),
               shape = 18, size = 3.8) +
    geom_text(data = null_res,
              aes(x = unit, y = observed_r + 0.014, label = sig,
                  colour = metric),
              position = position_dodge(width = 0.60),
              size = 3.0, vjust = 0, show.legend = FALSE) +
    scale_colour_manual(values = c(J = "#2C6E9B", B = "#D17B2F"),
                        labels = c(J = "Jaccard", B = "Bray-Curtis")) +
    scale_x_discrete(labels = c("V4 (within)" = "V4\n(within)",
                                "V5-V7 (within)" = "V5-V7\n(within)",
                                "V4-V5-V7 (shared)" = "V4-V5-V7\n(shared)")) +
    scale_y_continuous(limits = c(-0.02, 0.235),
                       breaks = seq(0, 0.2, by = 0.05)) +
    labs(x = NULL, y = "Mantel r", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
} else {
  pA <- ggplot(eff, aes(x = r, y = unit, colour = metric)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                    position = position_dodge(width = 0.35),
                    linewidth = 0.45, size = 0.35) +
    scale_colour_manual(values = c(J = "#2C6E9B", B = "#D17B2F"),
                        labels = c(J = "Jaccard", B = "Bray-Curtis")) +
    scale_x_continuous(limits = c(0, 0.22),
                       breaks = seq(0, 0.2, by = 0.05)) +
    labs(x = "Mantel r (95% CI)", y = NULL, colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
}

# ---- panel B ---------------------------------------------------------------
pB <- ggplot(short, aes(x = contrast, y = unit, colour = metric)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c(J = "#2C6E9B", B = "#D17B2F"),
                      labels = c(J = "Jaccard", B = "Bray-Curtis")) +
  scale_x_continuous(limits = c(0, 0.06),
                     breaks = seq(0, 0.05, by = 0.01)) +
  labs(x = "Far 50% - closest 10% mean distance", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

# ---- panel C ---------------------------------------------------------------
pC <- ggplot(sens, aes(x = r, y = scenario, colour = metric)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c(J = "#2C6E9B", B = "#D17B2F"),
                      labels = c(J = "Jaccard", B = "Bray-Curtis")) +
  scale_x_continuous(limits = c(0, 0.17),
                     breaks = seq(0, 0.15, by = 0.05)) +
  labs(x = "Mantel r", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

fig <- (pA | pB) / pC +
  plot_layout(heights = c(1, 0.9))

ggsave(file.path(OUT, "Fig1_cross_region.pdf"), fig,
       width = 10, height = 7.5, device = "pdf")
ggsave(file.path(OUT, "Fig1_cross_region.jpg"), fig,
       width = 10, height = 7.5, dpi = 300)

# ---- source data -----------------------------------------------------------
src <- rbind(
  eff[, .(panel = "A", unit = as.character(unit),
          metric = as.character(metric), value = r,
          low = CI_low, high = CI_high, p = p, n = n_bins,
          scenario = NA_character_)],
  short[, .(panel = "B", unit = as.character(unit),
            metric = as.character(metric), value = contrast,
            low = NA_real_, high = NA_real_, p = p, n = NA_integer_,
            scenario = NA_character_)],
  sens[, .(panel = "C", unit = "V4-V5-V7 (shared)",
           metric = as.character(metric), value = r,
           low = NA_real_, high = NA_real_, p = p, n = n_bins,
           scenario = as.character(scenario))]
)
fwrite(src, file.path(OUT, "source_data_Fig1.tsv"), sep = "\t")

# ---- main table ------------------------------------------------------------
fmt_p <- function(x) ifelse(is.na(x), NA_character_,
                            ifelse(x < 0.001, "<0.001",
                                   formatC(x, format = "f", digits = 3)))
fmt_r <- function(r, lo, hi) sprintf("%.3f (%.3f-%.3f)", r, lo, hi)
fmt_far <- function(x, p) sprintf("%.3f (%s)", x, fmt_p(p))

tab_long <- merge(
  eff[unit != "V4-V5-V7 (shared)",
      .(unit, n_bins, metric, r, CI_low, CI_high,
        p_mantel = p, p_strat)],
  t1r[, .(unit = unit_lab[region], metric = metric,
          p_mrm = p, r2_mrm = r_squared)],
  by = c("unit", "metric"), all.x = TRUE)
# cross-region MRM rows
t2r_rows <- t2r[, .(unit = "V4-V5-V7 (shared)", metric = metric,
                    p_mrm = p, r2_mrm = r_squared)]
tab_long <- rbind(
  tab_long,
  merge(eff[unit == "V4-V5-V7 (shared)"],
        t2r_rows, by = c("unit", "metric"), all.x = TRUE),
  fill = TRUE)
tab_long <- merge(
  tab_long,
  short[, .(unit, metric, contrast, p_short = p)],
  by = c("unit", "metric"), all.x = TRUE)
if (!is.null(null_res)) {
  tab_long <- merge(
    tab_long,
    null_res[, .(unit, metric, p_null = p_source_stratified,
                 null_mean, null_q05, null_q95)],
    by = c("unit", "metric"), all.x = TRUE)
} else {
  tab_long[, `:=`(p_null = NA_real_, null_mean = NA_real_,
                  null_q05 = NA_real_, null_q95 = NA_real_)]
}
tab_long[, unit := factor(unit, levels = unit_levels)]
tab_long[, metric := factor(metric, levels = c("J", "B"))]
setorder(tab_long, unit, metric)

make_row <- function(u) {
  j <- tab_long[unit == u & metric == "J"]
  b <- tab_long[unit == u & metric == "B"]
  data.table(
    Analysis = as.character(u),
    `n bins` = j$n_bins,
    `Jaccard r (95% CI)` = fmt_r(j$r, j$CI_low, j$CI_high),
    `Jaccard MRM P (R2)` = sprintf("%s (%.3f)", fmt_p(j$p_mrm), j$r2_mrm),
    `Jaccard null P` = fmt_p(j$p_null),
    `Jaccard far-close (P)` = fmt_far(j$contrast, j$p_short),
    `Bray-Curtis r (95% CI)` = fmt_r(b$r, b$CI_low, b$CI_high),
    `Bray-Curtis MRM P (R2)` = sprintf("%s (%.3f)", fmt_p(b$p_mrm), b$r2_mrm),
    `Bray-Curtis null P` = fmt_p(b$p_null),
    `Bray-Curtis far-close (P)` = fmt_far(b$contrast, b$p_short)
  )
}
tab <- rbindlist(lapply(unit_levels, make_row))
fwrite(tab, file.path(OUT, "Table1_cross_region.tsv"), sep = "\t")

# markdown with footnotes
md <- c(
  "**Table. Phylogenetic signal in host-distribution features of mapped 16S reference bins (V4 and V5-V7).**",
  "",
  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::kable(tab, format = "pipe", row.names = FALSE)
  } else {
    paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n")
  },
  "",
  "Notes:",
  "- Mapped bin = 97% identity SILVA reference bin (vsearch --usearch_global --id 0.97), not an exact ASV; multiple ASVs can collapse into one bin.",
  "- Within-region filters: prevalence >=20 samples, >=5 host species, >=3 sources; n = 955 (V4) and 1,487 (V5-V7).",
  "- Cross-region row: bins detected in both V4 and V5-V7; total filters as above plus >=5 positive samples in each region; n = 355 (983 shared before the per-region >=5 filter; 466 passed the total 20/5/3 filter).",
  "- Mantel tests: one-sided, 9,999 permutations; prevalence-stratified P = 0.001 for all rows; standard Mantel P < 1e-4 for all rows.",
  "- MRM controls |delta log prevalence| and |delta log mean relative abundance|; primary MRM P = 0.001 for all rows. Extended MRM adds host breadth and source-set distance; P = 0.001 for all rows.",
  "- Source-stratified null: observed Mantel r compared with 999 label permutations within strata of dominant source x prevalence quintile (within-region rows) or dominant source x total-prevalence quintile (shared-bin row). P values are shown in the null P columns; the null distribution is plotted in Fig_ASV_beta_phylo (top left).",
  "- Far-close = mean distance for the far 50% of phylogenetic pairs minus the closest 10%; P = 0.001 for all rows (prevalence-stratified label permutation).",
  "- Cross-region distribution distance is the equal-weight mean of the within-region Jaccard/Bray distances. Sensitivity (within-region n >= 10; pair min; pair max): Jaccard r = 0.097-0.136, Bray-Curtis r = 0.075-0.105; all P = 0.001.",
  "- P values are one-sided; 95% CI is the Mantel bootstrap interval."
)
writeLines(md, file.path(OUT, "Table1_cross_region.md"))

cat("Written figure and table to", OUT, "\n")
print(tab)
