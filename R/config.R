# Project-wide configuration / 项目级配置
# Run from the github2 project root, or set SEED_PROJECT_ROOT.

PROJECT_ROOT <- normalizePath(
  Sys.getenv("SEED_PROJECT_ROOT", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
if (!file.exists(file.path(PROJECT_ROOT, "data", "MANIFEST.tsv"))) {
  stop("Not in the github2 project root / 当前目录不是 github2 项目根目录: ", PROJECT_ROOT)
}

DATA_DIR     <- file.path(PROJECT_ROOT, "data")
SERVER_DIR   <- file.path(PROJECT_ROOT, "server_result")
RESULT_DIR   <- file.path(PROJECT_ROOT, "results")
DERIVED_DIR  <- file.path(RESULT_DIR, "derived")
CSS_DIR      <- file.path(DERIVED_DIR, "css")
TREE_DIR     <- file.path(DERIVED_DIR, "trees")
CKPT_DIR     <- file.path(RESULT_DIR, "checkpoints")
FIG_DIR      <- file.path(RESULT_DIR, "figures")
TAB_DIR      <- file.path(RESULT_DIR, "tables")
DIAG_DIR     <- file.path(RESULT_DIR, "diagnostics")
LOG_DIR      <- file.path(RESULT_DIR, "logs")
for (d in c(RESULT_DIR, DERIVED_DIR, CSS_DIR, TREE_DIR, CKPT_DIR, FIG_DIR, TAB_DIR, DIAG_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Threads are capped at 16 for every local or optional server-style step.
detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (is.na(detected) || detected < 1L) detected <- 1L
N_CORES <- max(1L, min(16L, as.integer(Sys.getenv("SEED_THREADS", unset = "16")), detected))
data.table::setDTthreads(N_CORES)
Sys.setenv(OMP_NUM_THREADS = N_CORES, OPENBLAS_NUM_THREADS = N_CORES,
           MKL_NUM_THREADS = N_CORES)
options(mc.cores = N_CORES)

# Recompute flags: empty = reuse checkpoints; "all" or comma-separated step names.
RECOMPUTE <- trimws(strsplit(Sys.getenv("SEED_RECOMPUTE", unset = ""), ",")[[1]])
RECOMPUTE <- RECOMPUTE[nzchar(RECOMPUTE)]
recompute_step <- function(step) {
  identical(tolower(step), "all") || "all" %in% RECOMPUTE || step %in% RECOMPUTE
}
ALLOW_SERVER_RECOMPUTE <- identical(tolower(Sys.getenv("ALLOW_SERVER_RECOMPUTE", "false")), "true")

# Approved analysis constants and seeds / 已批准的分析常数与种子
REGIONS_ALL <- c("v3_v4", "v4", "v5_v6", "v5_v7", "its1", "its2")
REGIONS3    <- c("v4", "v5_v7", "its1")
REGION_LABEL <- c(v4 = "V4", v5_v7 = "V5-V7", its1 = "ITS1", v3_v4 = "V3-V4",
                  v5_v6 = "V5-V6", its2 = "ITS2")
REGION_LAB  <- REGION_LABEL
SEEDS <- list(
  phylo_permutation = 20261106L,
  phylo_bootstrap   = 123L,
  gos                = 20260813L,
  sensitivity        = 20261110L,
  bnti               = 12345L,
  beta_conditional   = 20261110L
)
GOS_ALPHA <- 0.05
GOS_N_ITER <- 1000L
GOS_N_REP <- 1000L
