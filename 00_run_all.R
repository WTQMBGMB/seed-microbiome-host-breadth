#!/usr/bin/env Rscript
# =============================================================================
# Run the complete pipeline / 运行完整流程
#
# Usage / 用法:
#   Rscript 00_run_all.R
#   SEED_RECOMPUTE=css,alphabeta Rscript 00_run_all.R
#   SEED_ALLOW_SERVER_RECOMPUTE=true SEED_RECOMPUTE=all Rscript 00_run_all.R
# =============================================================================

root <- normalizePath(Sys.getenv("SEED_PROJECT_ROOT", unset = getwd()),
                      winslash = "/", mustWork = TRUE)
setwd(root)
Sys.setenv(SEED_PROJECT_ROOT = root)
source(file.path(root, "R", "config.R"))

rmd_order <- c(
  "00_preprocessing.Rmd",
  "Ch1_dataset.Rmd",
  "Ch2_host_effect.Rmd",
  "Ch3_phylogeny.Rmd",
  "Ch4_generality.Rmd",
  "Ch5_traits.Rmd"
)

# ---- 1. Preflight / 预检 ----------------------------------------------------
packages <- c("data.table", "ggplot2", "dplyr", "tidyr", "patchwork", "ape",
              "vegan", "lme4", "metafor", "phytools", "ecodist", "ggtree",
              "ggtreeExtra", "itol.toolkit", "Matrix", "sf", "ragg", "ggrepel",
              "ggnewscale",
              "ggpubr", "cowplot", "gridExtra", "scales", "rmarkdown",
              "Rfast", "metagenomeSeq", "nortest")
missing_pkg <- packages[!vapply(packages, requireNamespace, logical(1),
                                quietly = TRUE)]
if (length(missing_pkg)) {
  stop("Missing R packages / 缺少 R 包: ", paste(missing_pkg, collapse = ", "),
       "\nSee environment/README.md.", call. = FALSE)
}
if (!rmarkdown::pandoc_available()) {
  stop("pandoc is required / 需要 pandoc")
}

required <- c(
  file.path(root, "data", "metadata", "metadata_16S.txt"),
  file.path(root, "data", "metadata", "metadata_ITS.txt"),
  file.path(root, "data", "filtered_merge", "16S", "Merged_16S_table_filtered.txt"),
  file.path(root, "data", "filtered_merge", "ITS", "Merged_ITS_table_filtered.txt"),
  file.path(root, "server_result", "MANIFEST.tsv"),
  file.path(root, "manifests", "output_manifest.tsv")
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing inputs / 缺少输入:\n  ", paste(missing, collapse = "\n  "),
       "\nSee data/README.md.", call. = FALSE)
}

ensure_v2_dirs()
dir.create(file.path(root, "result", "html"), recursive = TRUE,
           showWarnings = FALSE)

# ---- 2. Render in order / 按顺序渲染 ---------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
log_file <- file.path(root, "result", "logs", "run_all.log")
cat("Pipeline start:", format(Sys.time()), "\n", file = log_file)
for (rmd in rmd_order) {
  message("\n=== ", rmd, " ===")
  t0 <- Sys.time()
  status <- system2(
    rscript,
    c("--vanilla", shQuote(file.path(root, "scripts", "render_one.R")), rmd,
      shQuote(file.path(root, "result", "html"))),
    stdout = TRUE, stderr = TRUE
  )
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  status_code <- attr(status, "status")
  if (is.null(status_code)) status_code <- 0L
  cat(sprintf("%s\t%s\t%.1f min\tstatus=%d\n", rmd, format(Sys.time()),
              elapsed, status_code), file = log_file, append = TRUE)
  if (status_code != 0L) {
    cat(status, sep = "\n")
    stop("Rendering failed / 渲染失败: ", rmd, call. = FALSE)
  }
  message("done in ", elapsed, " min")
}

# ---- 3. Output check / 输出检查 --------------------------------------------
checker <- file.path(root, "scripts", "tools", "check_outputs.py")
if (file.exists(checker)) {
  message("\n=== output check ===")
  status <- system2("python", shQuote(checker), stdout = TRUE, stderr = TRUE)
  cat(status, sep = "\n")
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop("Output check failed / 输出检查失败", call. = FALSE)
  }
}
cat("Pipeline finished:", format(Sys.time()), "\n", file = log_file, append = TRUE)
message("\nAll chapters completed / 全部章节完成")
