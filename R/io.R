# Shared I/O, checkpoint and figure export helpers / 共享 I/O、检查点与图件导出

server_path <- function(...) file.path(SERVER_DIR, ...)
derived_path <- function(...) file.path(DERIVED_DIR, ...)
ckpt_path <- function(chapter, ...) file.path(CKPT_DIR, chapter, ...)

assert_files <- function(paths, what = "input") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing ", what, " / 缺少", what, ":\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

# Verify a read-only server asset; never write into server_result.
use_server <- function(rel, required = TRUE) {
  p <- server_path(rel)
  if (!file.exists(p)) {
    if (required) {
      stop("Missing server result / 缺少 server_result 文件: ", p,
           "\nSet ALLOW_SERVER_RECOMPUTE=true and run the documented server script if a local recomputation is required.",
           call. = FALSE)
    }
    return(NULL)
  }
  p
}

# Run a block only when its output files are absent or the step is forced.
# The expression is evaluated in the caller's environment.
checkpoint <- function(paths, step, code, meta = list()) {
  if (all(file.exists(paths)) && !recompute_step(step)) {
    message("[checkpoint] reuse / 复用: ", paste(basename(paths), collapse = ", "))
    return(invisible(FALSE))
  }
  for (p in paths) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  force(code)
  assert_files(paths, "checkpoint output")
  meta_file <- paste0(paths[[1]], ".meta.json")
  info <- c(list(step = step, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 r_version = R.version.string, n_cores = N_CORES), meta)
  jsonlite::write_json(info, meta_file, auto_unbox = TRUE, pretty = TRUE)
  invisible(TRUE)
}

# One figure export implementation for every manuscript panel.
save_panel <- function(plot, stem, width, height, dpi = 300) {
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  jpg_file <- file.path(FIG_DIR, paste0(stem, ".jpg"))
  ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                  device = "pdf", bg = "white")
  ggplot2::ggsave(jpg_file, plot, width = width, height = height,
                  device = "jpeg", dpi = dpi, quality = 95, bg = "white")
  message("[save_panel] ", basename(pdf_file), " / ", basename(jpg_file))
  invisible(c(pdf = pdf_file, jpg = jpg_file))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}
