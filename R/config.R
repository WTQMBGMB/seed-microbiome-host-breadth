# =============================================================================
# Project paths and run switches / 项目路径与运行开关
#
# All paths are resolved from the repository root.  No machine-specific or
# server paths are allowed here. / 所有路径相对于仓库根目录，不包含本机或服务器绝对路径。
# =============================================================================

PROJECT_ROOT <- normalizePath(
  Sys.getenv("SEED_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

DATA_DIR      <- file.path(PROJECT_ROOT, "data")
SERVER_DIR    <- file.path(PROJECT_ROOT, "server_result")
RESULT_DIR    <- file.path(PROJECT_ROOT, "result")
PREPARED_DIR  <- file.path(RESULT_DIR, "prepared")
MANIFEST_DIR  <- file.path(PROJECT_ROOT, "manifests")

CHAPTERS <- c("ch1_dataset", "ch2_host_effect", "ch3_phylogeny",
              "ch4_generality", "ch5_traits")

# ---- Run switches / 运行开关 -------------------------------------------------
# SEED_RECOMPUTE: comma-separated steps (css, alphabeta, traits, all).
# 重算指定步骤；默认读取仓库内已准备好的文件。
RECOMPUTE_STEPS <- strsplit(tolower(Sys.getenv("SEED_RECOMPUTE", "")), ",",
                            fixed = TRUE)[[1]]
RECOMPUTE_STEPS <- trimws(RECOMPUTE_STEPS[nzchar(RECOMPUTE_STEPS)])

want_recompute <- function(step) {
  "all" %in% RECOMPUTE_STEPS || step %in% RECOMPUTE_STEPS
}

# Server-equivalent steps are read-only unless explicitly allowed.
# 服务器等价步骤默认只读，只有显式允许时才本地重算。
ALLOW_SERVER_RECOMPUTE <- tolower(Sys.getenv("SEED_ALLOW_SERVER_RECOMPUTE", "false")) %in%
  c("true", "1", "yes")

SEED_NCORE <- suppressWarnings(as.integer(Sys.getenv("SEED_NCORE", "16")))
if (is.na(SEED_NCORE) || SEED_NCORE < 1L) SEED_NCORE <- 1L
SEED_NCORE <- min(SEED_NCORE, 16L)

# ---- Directory helpers / 目录助手 -------------------------------------------
ch_dir <- function(chapter, kind = c("figures", "tables")) {
  match.arg(kind)
  file.path(RESULT_DIR, chapter, kind)
}

int_dir <- function(chapter = NULL, shared = FALSE) {
  if (shared) file.path(RESULT_DIR, "intermediate", "shared")
  else file.path(RESULT_DIR, "intermediate", chapter)
}

rds_dir <- function(chapter = NULL, shared = FALSE) {
  if (shared) file.path(RESULT_DIR, "rds", "shared")
  else file.path(RESULT_DIR, "rds", chapter)
}

log_dir <- function(chapter) file.path(RESULT_DIR, "logs", chapter)

tree_dir <- function(chapter) file.path(RESULT_DIR, chapter, "tree")

ensure_v2_dirs <- function() {
  for (ch in CHAPTERS) {
    for (d in c(ch_dir(ch, "figures"), ch_dir(ch, "tables"),
                int_dir(ch), rds_dir(ch), log_dir(ch), tree_dir(ch))) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    }
  }
  for (d in c(int_dir(shared = TRUE), rds_dir(shared = TRUE), PREPARED_DIR)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  invisible(TRUE)
}

# ---- Input resolution / 输入路径解析 ----------------------------------------
# Prepared tables live either in data/ (shipped with the repository) or in
# result/prepared/ (regenerated when the shipped copy is absent).
# 优先使用随仓库提供的表；缺失时使用重算结果。
first_existing <- function(paths, what) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(paste0(
      what, " not found. Expected one of:\n  ",
      paste(paths, collapse = "\n  "),
      "\nRun 00_preprocessing.Rmd or provide the file.\n",
      what, " 未找到；请运行预处理或补充文件。"), call. = FALSE)
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

# When a step is explicitly recomputed, the regenerated copy is preferred.
# 显式重算时优先使用 result/prepared 下的新副本。
css_merged_path <- function(file) {
  a <- file.path(DATA_DIR, "standard_merge", file)
  b <- file.path(PREPARED_DIR, "standard_merge", file)
  first_existing(if (want_recompute("css")) c(b, a) else c(a, b),
                 paste("merged CSS table", file))
}

css_region_path <- function(region, file) {
  a <- file.path(DATA_DIR, "standard_region", region, file)
  b <- file.path(PREPARED_DIR, "standard_region", region, file)
  first_existing(if (want_recompute("css")) c(b, a) else c(a, b),
                 paste("region CSS table", region, file))
}

host_tree_path <- function(file) {
  first_existing(c(file.path(DATA_DIR, "plant_tree", file),
                   file.path(PREPARED_DIR, "plant_tree", file)),
                 paste("host tree", file))
}

precomp_candidates <- function(file, step) {
  a <- file.path(DATA_DIR, "precomputed", file)
  b <- file.path(PREPARED_DIR, "precomputed", file)
  if (want_recompute(step)) c(b, a) else c(a, b)
}

precomp_path <- function(file, step) {
  first_existing(precomp_candidates(file, step),
                 paste("precomputed file", file))
}

# Beta distances are server results; alpha caches remain local precomputed.
# β 距离为服务器结果；alpha 缓存仍属本地预计算。
bray_path <- function(region) {
  file <- sprintf("bray_dist_%s.rds", region)
  srv <- server_path("ch2", "beta", file)
  dat <- file.path(DATA_DIR, "precomputed", file)
  pre <- file.path(PREPARED_DIR, "precomputed", file)
  first_existing(if (want_recompute("alphabeta")) c(pre, srv, dat)
                 else c(srv, dat, pre),
                 paste("Bray-Curtis distance", region))
}

server_path <- function(...) file.path(SERVER_DIR, ...)

meta_path <- function(dataset) {
  first_existing(file.path(DATA_DIR, "metadata", sprintf("metadata_%s.txt", dataset)),
                 paste("metadata", dataset))
}

require_server <- function(path, script_hint, what = basename(path)) {
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  stop(sprintf(
    "Server result '%s' is missing.\n%s\nProvide the file under server_result/ or run scripts/server/%s with SEED_ALLOW_SERVER_RECOMPUTE=true.\n服务器结果缺失：%s；请补充 server_result/ 文件，或按提示运行服务器脚本。",
    what, path, script_hint, what), call. = FALSE)
}
