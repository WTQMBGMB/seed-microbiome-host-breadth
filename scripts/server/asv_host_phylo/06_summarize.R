#!/usr/bin/env Rscript
# =============================================================================
# 06_summarize.R — compact, human-readable summary of the result tables
# Output: result_table/key_results_summary.tsv
#         result_table/RESULTS.md
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("BASE_DIR", "<SERVER_BASE>"),
                 "result/ASV_host_phylo/scripts/asv_host_phylo_utils.R"))

TAB <- file.path(OUT_ROOT, "result_table")
mt  <- fread(file.path(TAB, "mantel_results.tsv"))
mrm <- fread(file.path(TAB, "mrm_results.tsv"))
cg  <- fread(file.path(TAB, "correlogram.tsv"))
cp  <- fread(file.path(TAB, "close_pair_contrast.tsv"))
sens <- fread(file.path(TAB, "sensitivity_results.tsv"))
spear <- fread(file.path(TAB, "mantel_spearman.tsv"))
ext <- fread(file.path(TAB, "mrm_extended_results.tsv"))

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt <- function(x, d = 3) sprintf(paste0("%.", d, "f"), x)

main <- merge(
  mt[, .(region, metric, n_asv, n_pairs, r, p_standard, p_stratified,
         SES_standard, SES_stratified, p_adj_BH)],
  mrm[term == "dp", .(region, metric, mrm_coef = coef, mrm_p = p,
                      r_squared, r_squared_p)],
  by = c("region", "metric"), all.x = TRUE)
main[, `:=`(p_std_f = fmt_p(p_standard), p_strat_f = fmt_p(p_stratified),
            mrm_p_f = fmt_p(mrm_p), r_f = fmt(r), SES_f = fmt(SES_standard, 2),
            R2_f = fmt(r_squared))]
setorder(main, region, metric)

bin1 <- cg[bin == 1, .(region, metric, phylo_mid, n_pairs, obs_mean, null_mean,
                       diff = obs_mean - null_mean, SES, p_lower)]
bin1[, `:=`(diff_f = fmt(diff), SES_f = fmt(SES, 1), p_f = fmt_p(p_lower),
            obs_f = fmt(obs_mean), null_f = fmt(null_mean))]
setorder(bin1, region, metric)

cp2 <- copy(cp)
cp2[, `:=`(contrast_f = fmt(contrast), SES_f = fmt(SES, 1), p_f = fmt_p(p))]
setorder(cp2, region, metric)

ext_phylo <- ext[term == "dp", .(region, analysis, subset, metric,
                                 coef, p, r_squared)]
ext_phylo[, `:=`(p_f = fmt_p(p), coef_f = fmt(coef, 4))]
setorder(ext_phylo, region, analysis, metric)

# Long-format key summary
long <- rbindlist(list(
  main[, .(section = "main_mantel_mrm", region, metric,
           term = "Mantel_r", value = r, p = p_standard,
           note = sprintf("p_strat=%s; SES=%.2f; MRM_phylo_p=%s; R2=%.3f",
                          fmt_p(p_stratified), SES_standard,
                          fmt_p(mrm_p), r_squared))],
  bin1[, .(section = "correlogram_bin1", region, metric,
           term = "obs_minus_null", value = diff, p = p_lower,
           note = sprintf("obs=%.3f; null=%.3f; SES=%.1f; n_pairs=%d",
                          obs_mean, null_mean, SES, n_pairs))],
  cp2[, .(section = "close_vs_far", region, metric,
          term = "contrast", value = contrast, p = p,
          note = sprintf("SES=%.1f; mean_close=%.3f; mean_far=%.3f",
                         SES, mean_close, mean_far))],
  sens[, .(section = "sensitivity", region, metric,
           term = paste(analysis, subset, sep = ":"), value = r, p = p,
           note = sprintf("n_asv=%d; n_units=%d; SES=%.1f",
                          n_asv, n_units, SES))],
  spear[, .(section = "spearman", region, metric,
            term = "Spearman_rho", value = rho, p = p,
            note = sprintf("SES=%.1f", SES))],
  ext_phylo[, .(section = "extended_mrm_phylo", region, metric,
                term = paste(analysis, subset, sep = ":"), value = coef, p = p,
                note = sprintf("R2=%.3f", r_squared))]
), fill = TRUE)
fwrite(long, file.path(TAB, "key_results_summary.tsv"), sep = "\t")

# Markdown summary
md <- c(
  "# ASV distribution × microbial phylogeny — key results",
  "",
  "## Main Mantel / MRM (main filter 20/5/3)",
  "",
  "| Region | Metric | n ASV | Mantel r | P | P (prevalence-stratified) | SES | MRM phylo P | MRM R² |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|",
  sprintf("| %s | %s | %d | %.3f | %s | %s | %.2f | %s | %.3f |",
          main$region, main$metric, main$n_asv, main$r,
          main$p_std_f, main$p_strat_f, main$SES_standard,
          main$mrm_p_f, main$r_squared),
  "",
  "## Correlogram, closest phylo bin (bin 1)",
  "",
  "| Region | Metric | Phylo mid | Observed | Null | Obs − null | SES | P(lower) |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  sprintf("| %s | %s | %.4f | %.3f | %.3f | %.3f | %.1f | %s |",
          bin1$region, bin1$metric, bin1$phylo_mid, bin1$obs_mean,
          bin1$null_mean, bin1$diff, bin1$SES, bin1$p_f),
  "",
  "## Close-vs-far contrast (far 50% − closest 10%)",
  "",
  "| Region | Metric | Contrast | SES | P |",
  "|---|---|---:|---:|---:|",
  sprintf("| %s | %s | %.4f | %.1f | %s |",
          cp2$region, cp2$metric, cp2$contrast, cp2$SES, cp2$p_f),
  "",
  "## Sensitivity (Mantel r)",
  "",
  "| Region | Metric | Analysis | Subset | n ASV | r | P |",
  "|---|---|---|---|---:|---:|---:|",
  sprintf("| %s | %s | %s | %s | %d | %.3f | %s |",
          sens$region, sens$metric, sens$analysis, sens$subset,
          sens$n_asv, sens$r, fmt_p(sens$p)),
  "",
  "## Spearman Mantel",
  "",
  "| Region | Metric | rho | P |",
  "|---|---|---:|---:|",
  sprintf("| %s | %s | %.3f | %s |",
          spear$region, spear$metric, spear$rho, fmt_p(spear$p)),
  "",
  "## Extended MRM, phylo term",
  "",
  "| Region | Metric | Analysis | Coef | P | R² |",
  "|---|---|---|---:|---:|---:|",
  sprintf("| %s | %s | %s | %.4f | %s | %.3f |",
          ext_phylo$region, ext_phylo$metric, ext_phylo$analysis,
          ext_phylo$coef, ext_phylo$p_f, ext_phylo$r_squared)
)
writeLines(md, file.path(TAB, "RESULTS.md"))
cat(md, sep = "\n")
cat("\nDONE 06_summarize\n")
