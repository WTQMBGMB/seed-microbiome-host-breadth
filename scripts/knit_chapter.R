#!/usr/bin/env Rscript
# Knit one chapter fragment for chapter-level validation.
# 单独 Knit 一个章节片段，便于逐章验证。
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/knit_chapter.R analysis/01_dataset_overview.Rmd")
chapter <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
wrapper <- file.path(root, "results", "logs", paste0("_", basename(chapter), ".Rmd"))
dir.create(dirname(wrapper), recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "---",
  paste0("title: \"", basename(chapter), "\""),
  "output:",
  "  html_document:",
  "    self_contained: true",
  "    code_folding: hide",
  "---",
  "",
  paste0("```{r child='", chapter, "'}"),
  "```"
), wrapper, useBytes = TRUE)
Sys.setenv(SEED_PROJECT_ROOT = root)
rmarkdown::render(wrapper, quiet = FALSE, envir = new.env(parent = globalenv()))
