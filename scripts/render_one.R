#!/usr/bin/env Rscript
# Render one Rmd in an isolated R process / 在独立 R 进程中渲染单个 Rmd
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: render_one.R <rmd> [output_dir]")
rmd <- args[1]
out <- if (length(args) >= 2L) args[2] else file.path("result", "html")
root <- Sys.getenv("SEED_PROJECT_ROOT", unset = getwd())
setwd(root)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Figure links are written relative to the rendered document.
# 图件链接按渲染文档所在目录写成相对路径。
options(seed.html_dir = normalizePath(out, winslash = "/"))
rmarkdown::render(rmd, output_dir = out, quiet = TRUE)
