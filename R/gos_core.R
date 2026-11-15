# ============================================================================
# GOS generalist/specialist core functions / GOS 通才-专才核心函数
#
# 本文件为 GOS 分析的全部函数库（基于旧方法 B：gos_presence 工作流重构）。
# 所有函数通过【参数】接收输入 / 输出文件路径，不在函数内硬编码文件地址；
# 具体文件地址由调用方（Analysis_pipeline.Rmd 第 8 节）传入。
#
# ── 工作流（4 步）────────────────────────────────────────────────────────
#   Step 0  宿主系统发育距离矩阵（两数据集共用）
#             gos_host_dist()
#   Step 1  0-1 矩阵提取与合并
#             gos_presence_merge()
#   Step 2  n_hosts 实际分布 + 零分布 + 通才阈值 + 比较分布图
#             gos_nhost_null() → gos_nhost_threshold() → gos_plot_nhost()
#   Step 3  z 值 + 显著性（基于宿主发育距离的每-n_hosts 零模型 + BH-FDR）
#             gos_mpd_null() → gos_z_scores()
#   Step 4  分类（通才 / 专才 / 中间型 / 缺失）+ 参考序列与物种注释提取
#             gos_classify() → gos_extract_fasta() → gos_extract_tax()
#   Step 5  n_host–相对丰度回归曲线（通才/专才/其他点形 + 线性拟合）
#             gos_nhost_abund_data() → gos_nhost_abund_rank()
#             → gos_nhost_abund_subplot() → gos_nhost_abund_combine_grobs()
#             → gos_nhost_abund_plot()
#   Step 6  通才/专才分类群组成差异（Genus~Class 双向柱状图 + Fisher）
#             gos_composition()
#   Step 7  dispersed/clustered ASV 宿主分布热图 + clustered 物种注释表
#             gos_dispersion_heatmap()
#
# ── 目录约定 ──────────────────────────────────────────────────────────────
#   result/GOS/Temp_data/       中间文件（RDS、日志、调试代码）
#   result/GOS/result_table/    最终结果表（合并矩阵、ASV 信息表、通才/专才表等）
#   result/GOS/result_figure/   最终图形
#   （目录由 gos_init_dirs() 创建）
#
# ── 统计方法要点 ─────────────────────────────────────────────────────────
#   Step 2 零模型：对合并后的 ASV × 宿主物种 0-1 矩阵做 within-column 打乱
#            （每列 0/1 内部重排、保持列序 → 各宿主代表物种数 T 不变），
#            N 次迭代得 n_hosts 零分布；阈值 = 分布图上"实际分布首次大于
#            零分布"的 n_hosts（逐 n_hosts 实际计数 > 零期望均值；零分布随机
#            铺展在低 n_hosts 记录更多、实际长尾在高 n_hosts 保留更多）。
#   Step 3 零模型：对每个出现的 n_hosts = k，随机抽 k 个宿主物种的 MPD
#            （宿主两两遗传距离均值）分布，得 null_mu[k] / null_sd[k]；
#            z = (mean_dist − null_mu[k]) / null_sd[k]；
#            p_disp = P(Z ≥ z)（分散），p_clus = P(Z ≤ z)（聚集），
#            各自 BH-FDR 校正。多检验下不可用固定 |z|>1.96。
#   Step 4 通才 = n_hosts ≥ 阈值（广度达标即为通才；显著聚集不再剥夺通才资格，
#            只作为 dispersion 属性标注）；专才 = n_hosts == 1；
#            中间型 = 2 ≤ n_hosts < 阈值；n_hosts == 0 = 无代表性宿主（NRH，
#            No representative host，非"未检出"）。
#            附加 dispersion 属性列：q_disp < alpha = dispersed（显著分散）、
#            q_clus < alpha = clustered、均不显著 = neutral；显著分散仅为标注属性，
#            不作为通才的硬性门槛。
# ============================================================================

# ---- 0. 通用小工具 ----

# 目录初始化：创建 GOS 输出目录结构，返回命名目录向量
#   dirs["fig"]  —— result_figure（最终图形）
#   dirs["tab"]  —— result_table（最终结果表）
#   dirs["tmp"]  —— Temp_data（中间文件 + 调试代码）
gos_init_dirs <- function(gos_root) {
  dirs <- c(fig = file.path(gos_root, "result_figure"),
            tab = file.path(gos_root, "result_table"),
            tmp = file.path(gos_root, "Temp_data"))
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  invisible(dirs)
}

# 参考序列 fasta 头 → ASV ID 集合（Silva / Unite 映射金标准）
gos_read_ref_ids <- function(fa_path) {
  lines <- readLines(fa_path, warn = FALSE)
  sub("^>", "", lines[grepl("^>", lines)])
}

# ---- Step 0：宿主系统发育距离矩阵（两数据集共用） ----
# tree_file  : 物种树 .nwk
# out_rds    : 输出（list(tree, dist)，dist 为 cophenetic 归一化到 0-1）
gos_host_dist <- function(tree_file, out_rds, force = FALSE) {
  if (file.exists(out_rds) && !force) {
    host_dist <- readRDS(out_rds)
    cat(sprintf("[gos-hostdist] 检查点命中：%s 已存在，跳过计算\n", basename(out_rds)))
  } else {
    tree <- ape::read.tree(tree_file)
    dmat <- stats::cophenetic(tree)          # species × species
    dmat <- dmat / max(dmat)                 # 归一化到 0-1
    host_dist <- list(tree = tree, dist = dmat)
    saveRDS(host_dist, out_rds)
    cat(sprintf("[gos-hostdist] 宿主距离矩阵已保存（%d 个物种，标准化 0-1）→ %s\n",
                nrow(dmat), out_rds))
  }
  invisible(host_dist)
}

# ---- Step 1：0-1 矩阵提取与合并 ----
# 逻辑（与旧 gos_presence_merge 完全一致）：
#   1) 仅保留映射到参考序列（Silva/Unite）的 ASV；
#   2) 丰度 ≠ 0 → 1，得 样本 × ASV 0-1 存在矩阵；
#   3) 对每个宿主物种用 Anderson–Darling 正态性检验决定代表物种数 T
#        （样本数≥10 且正态 p>0.05 → T=均值；≥10 且非正态 → T=中位数；
#          <10 → T=均值；T=as.integer，限幅到 [0, n_ASV]），
#      按该物种样本内各 ASV 流行度降序取前 T 个为代表物种；
#   4) 合并为 ASV × 宿主物种 0-1 矩阵，与宿主树取交集（剔除树外物种）。
# 输出：
#   out_merged_txt —— 合并 0-1 矩阵文本（首列 ASV_ID，result_table/presence_merged_{16s,its}.txt）
#   out_log_csv    —— 每物种 样本数/正态性 p/T/实际代表数/树匹配 日志
#   out_rds        —— 中间 RDS（含 merged 矩阵，供 Step 2/3 使用）
gos_presence_merge <- function(css_file, meta_file, ref_fa, host_dist_rds,
                               out_merged_txt, out_log_csv, out_rds,
                               ds, force = FALSE) {
  if (file.exists(out_rds) && !force) {
    cat(sprintf("[gos-merge:%s] 检查点命中：%s 已存在，跳过合并\n", ds, basename(out_rds)))
    return(readRDS(out_rds))
  }
  dmat <- readRDS(host_dist_rds)$dist

  # 读取 metadata 与 CSS 表
  meta <- as.data.frame(data.table::fread(meta_file, sep = "\t", header = TRUE))
  raw  <- data.table::fread(css_file, sep = "\t", header = TRUE)
  asv_ids <- as.character(raw[[1]])
  if (anyDuplicated(asv_ids) > 0) {
    stop(sprintf("[gos-merge:%s] ASV ID 存在重复，请检查数据\n", ds))
  }
  sample_cols <- colnames(raw)[-1]

  # 样本列 → 宿主物种（metadata species 列；R 列名 X 前缀兜底）
  sp <- meta$species[match(sample_cols, meta$SampleID)]
  if (anyNA(sp)) {
    fix <- which(is.na(sp))
    sp[fix] <- meta$species[match(gsub("^X", "", sample_cols[fix]), meta$SampleID)]
  }
  if (anyNA(sp)) {
    n_un <- sum(is.na(sp))
    cat(sprintf("[gos-merge:%s] 警告：%d 个样本列未匹配到 metadata 的 SampleID，已剔除\n", ds, n_un))
    keep <- !is.na(sp)
    sp   <- sp[keep]
    raw  <- raw[, c(1, which(keep) + 1), with = FALSE]
    sample_cols <- sample_cols[keep]
  }

  # 仅保留映射到参考序列的 ASV
  ref_ids <- gos_read_ref_ids(ref_fa)
  keep_asv <- asv_ids %in% ref_ids
  raw <- raw[keep_asv]
  cat(sprintf("[gos-merge:%s] 参考序列映射：%d/%d 个 ASV 保留\n",
              ds, sum(keep_asv), length(asv_ids)))
  if (nrow(raw) == 0) stop(sprintf("[gos-merge:%s] 没有映射到参考序列的 ASV，无法继续\n", ds))

  # 0-1 存在矩阵（丰度 ≠ 0 → TRUE）
  bin <- as.matrix(raw[, -1, with = FALSE]) > 0
  rownames(bin) <- asv_ids[keep_asv]
  rm(raw); invisible(gc())
  cat(sprintf("[gos-merge:%s] 0-1 存在矩阵：%d ASV × %d 样本\n", ds, nrow(bin), ncol(bin)))

  spp_levels <- unique(sp)
  n_spp <- length(spp_levels)

  # 对每个宿主物种提取代表物种（参考 Extract_representative_species）
  rep_list  <- vector("list", n_spp); names(rep_list)  <- spp_levels
  T_used    <- integer(n_spp);       names(T_used)    <- spp_levels
  pval_used <- rep(NA_real_, n_spp); names(pval_used) <- spp_levels

  for (s in spp_levels) {
    cols <- which(sp == s)
    temp <- bin[, cols, drop = FALSE]
    prevalence <- rowSums(temp)                # 各 ASV 在该宿主样本中的出现次数
    sample_num <- colSums(temp)                # 各样本的 ASV 数（样本丰富度）
    n_samp <- length(sample_num)
    p_val <- NA_real_
    if (n_samp >= 10) {
      # Anderson–Darling 正态性检验（nortest 要求 n ≥ 8）；常数向量时按非正态处理
      p_val <- tryCatch(nortest::ad.test(as.numeric(sample_num))$p.value,
                        warning = function(w) NA_real_,
                        error   = function(e) NA_real_)
      if (is.na(p_val)) p_val <- 0
      T_s <- if (p_val > 0.05) as.integer(mean(sample_num))
             else as.integer(stats::median(sample_num))
    } else {
      T_s <- as.integer(mean(sample_num))
    }
    T_s <- min(max(T_s, 0L), nrow(bin))        # 防御：T 限制在 [0, n_ASV]
    reps <- if (T_s >= 1L) rownames(temp)[order(-prevalence)][seq_len(T_s)]
            else character(0)                  # 流行度降序取前 T 个 ASV
    rep_list[[s]]  <- reps
    T_used[s]      <- T_s
    pval_used[s]   <- p_val
  }

  # 合并为 ASV × 宿主物种 0-1 表（1 = 该 ASV 是此宿主的代表物种）
  merged <- matrix(0L, nrow = nrow(bin), ncol = n_spp,
                   dimnames = list(rownames(bin), spp_levels))
  for (s in spp_levels) merged[rep_list[[s]], s] <- 1L
  rm(bin); invisible(gc())

  # 交集策略：只保留树中物种（其余无法计算遗传距离）
  present <- intersect(spp_levels, colnames(dmat))
  dropped <- setdiff(spp_levels, colnames(dmat))
  merged  <- merged[, present, drop = FALSE]
  cat(sprintf(paste0("[gos-merge:%s] 物种匹配：%d/%d（%.1f%%）在总物种树中；",
                     "未匹配 %d 个已剔除\n"),
              ds, length(present), length(spp_levels),
              100 * length(present) / length(spp_levels), length(dropped)))

  # 输出：合并 0-1 矩阵文本（result_table/presence_merged_{16s,its}.txt）
  out_df <- data.frame(ASV_ID = rownames(merged),
                       as.data.frame(merged, check.names = FALSE),
                       check.names = FALSE)
  data.table::fwrite(out_df, out_merged_txt, sep = "\t")
  cat(sprintf("[gos-merge:%s] 合并 0-1 矩阵已保存：%d ASV × %d 宿主 → %s\n",
              ds, nrow(merged), ncol(merged), out_merged_txt))

  # 输出：合并日志（每物种样本数 + 正态性 p + 代表物种数 T + 树匹配状态）
  spp_n <- table(sp)
  log_df <- data.frame(
    species     = spp_levels,
    n_samples   = as.integer(spp_n[spp_levels]),
    normality_p = round(pval_used, 4),         # Anderson–Darling p（<10 样本或常数向量为 0）
    rep_T       = T_used,                      # 使用的代表物种数
    n_rep       = vapply(rep_list, length, integer(1)),  # 实际提取的代表物种数
    in_tree     = spp_levels %in% colnames(dmat),
    stringsAsFactors = FALSE)
  write.csv(log_df, out_log_csv, row.names = FALSE)

  out <- list(ds = ds, merged = merged, present = present, dropped = dropped,
              spp_n = spp_n, n_asv = nrow(merged), n_samples = length(sp))
  saveRDS(out, out_rds)
  cat(sprintf(paste0("[gos-merge:%s] 完成：%d 个样本 → %d 个宿主物种，%d 个（映射）ASV，",
                     "%d 个 ASV 至少为一个宿主物种的代表\n"),
              ds, length(sp), ncol(merged), nrow(merged), sum(rowSums(merged) > 0)))
  invisible(out)
}

# ---- Step 2：n_hosts 实际分布 + 零分布 ----
# 零模型：within-column 打乱（保持各宿主代表物种数 col_sums 不变）。
#   每迭代对第 j 列随机抽 col_sums[j] 个 ASV 行置 1 → 破坏每个 ASV 的存在向量，
#   同时随机化 n_hosts；列和不变 ⇒ 零分布与实际分布总出现次数相同。
# 高效实现：逐列 sample.int + 累加，再 tabulate 成 n_hosts 计数直方图（0..n_cols）。
# 输出 RDS：actual_count / null_count（n_iter 累积）/ n_iter / n_hosts（逐 ASV）
gos_nhost_null <- function(merged_rds, out_rds, n_iter = 1000,
                           seed = 20260813, force = FALSE) {
  if (file.exists(out_rds) && !force) {
    cat(sprintf("[gos-nhost] 检查点命中：%s 已存在，跳过零分布\n", basename(out_rds)))
    return(readRDS(out_rds))
  }
  mg <- readRDS(merged_rds)
  bin <- mg$merged
  n_asv <- nrow(bin); n_cols <- ncol(bin)
  col_sums <- colSums(bin)                     # 各宿主物种代表物种数 T
  nhosts_actual <- rowSums(bin)
  actual_count <- tabulate(nhosts_actual + 1L, nbins = n_cols + 1L)

  # 零分布循环
  null_count <- integer(n_cols + 1L)
  set.seed(seed)
  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  start_time <- Sys.time()
  for (iter in seq_len(n_iter)) {
    cnt <- integer(n_asv)
    for (j in seq_len(n_cols)) {
      if (col_sums[j] >= 1L) {
        sel <- sample.int(n_asv, col_sums[j])  # 第 j 列随机 T_j-子集
        cnt[sel] <- cnt[sel] + 1L
      }
    }
    null_count <- null_count + tabulate(cnt + 1L, nbins = n_cols + 1L)
    utils::setTxtProgressBar(pb, iter)
    if (iter %% 200 == 0 && iter < n_iter) {
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      eta_min <- round(as.numeric(elapsed) / iter * (n_iter - iter) / 60, 1)
      cat(sprintf("  [%d/%d] 已耗时 %.1f min, 预计剩余 %.1f min\n",
                  iter, n_iter, round(as.numeric(elapsed) / 60, 1), eta_min))
    }
  }
  close(pb)
  elapsed_total <- difftime(Sys.time(), start_time, units = "mins")

  res <- list(ds = mg$ds,
              actual_count = actual_count,     # 长度 n_cols+1：n_hosts = 0..n_cols
              null_count   = null_count,       # n_iter 次迭代累积
              n_iter       = n_iter,
              n_asv        = n_asv,
              n_cols       = n_cols,
              col_sums     = col_sums,
              n_hosts      = nhosts_actual,    # 逐 ASV（阈值候选列出用）
              elapsed_min  = round(as.numeric(elapsed_total), 1))
  saveRDS(res, out_rds)
  cat(sprintf("[gos-nhost:%s] 完成：%d 次迭代（%.1f min），n_hosts 零分布已保存 → %s\n",
              mg$ds, n_iter, as.numeric(elapsed_total), out_rds))
  invisible(res)
}

# n_hosts 通才阈值（实际>零分布 crossing 法）
#   逐 n_hosts 比较：阈值 = 最小的 k（k ≥ 1）使 actual_count[k] > null_count[k]/n_iter，
#   即分布图上实际计数首次超过零分布均值的位置。零分布 within-column 随机铺展在低
#   n_hosts 处记录更多，而实际长尾在高 n_hosts 处保留更多 ASV ⇒ 首个 crossing 点
#   即为通才起始阈值（本数据 16S/ITS 均为 5）。
# 退化处理：若全程无 actual > null（极端数据），退回零分布 95% 分位。
gos_nhost_threshold <- function(null_res, ds = "") {
  act <- null_res$actual_count
  nul <- null_res$null_count
  nit <- null_res$n_iter
  max_t <- length(act) - 1L
  ts <- 0:max_t
  null_avg <- nul / nit                        # 零分布每迭代均值（n_hosts 0..max_t）
  above <- which(act[2:length(act)] > null_avg[2:length(null_avg)])  # k ≥ 1 的 crossing
  if (length(above) == 0) {
    # 退化：全程 actual ≤ null → 退回零分布 95% 分位
    q95 <- as.numeric(quantile(rep(ts, times = nul), probs = 0.95))
    thr <- max(1L, ceiling(q95))
    method <- "null_quantile95_fallback"
    cat(sprintf("[gos-thr:%s] 警告：无 actual>null 的 crossing，退回零分布 95%% 分位：%d\n", ds, thr))
  } else {
    thr <- above[1L]                           # 首个 actual > null_avg 的 n_hosts
    method <- "actual_gt_null"
    cat(sprintf("[gos-thr:%s] 阈值 = 首个 actual>null_avg 的 n_hosts：%d\n", ds, thr))
  }
  A <- rev(cumsum(rev(act)))                   # A(t)：n_hosts ≥ t 的实际 ASV 数
  N <- rev(cumsum(rev(nul))) / nit             # N(t)：零期望 ≥ t 的记录数
  list(threshold = thr,
       method = method,
       n_generalist_candidate = A[thr + 1L],   # n_hosts ≥ 阈值 的实际 ASV 数
       null_q95 = as.numeric(quantile(rep(ts, times = nul), probs = 0.95)),
       df = data.frame(n_hosts = ts, actual = A, null_expected = N, excess = A - N))
}

# n_hosts 实际 vs 零分布对比图（两数据集并排，阈值虚线标注）
gos_plot_nhost <- function(res_16s, res_its, thr_16s, thr_its,
                           y_scale = "log10", x_max = 50,
                           out_pdf, out_jpg) {
  # y_scale: "log10" | "log2" | "identity"（不缩放，原始计数）
  # x_max:   x 轴展示范围上限（n_hosts 超过该值后分布几乎为 0，默认截断到 50）
  stopifnot(y_scale %in% c("log10", "log2", "identity"))

  panel <- function(res, thr, color) {
    act <- res$actual_count
    nul <- res$null_count
    nit <- res$n_iter
    max_t <- min(length(act) - 1L, x_max)     # x 轴截断到 0~x_max
    ts <- 1:max_t                             # 图仅展示 n_hosts ≥ 1（0 主导全轴）
    df <- rbind(
      data.frame(n_hosts = ts, count = as.numeric(act[ts + 1L]),
                 distribution = "Actual", stringsAsFactors = FALSE),
      data.frame(n_hosts = ts, count = nul[ts + 1L] / nit,
                 distribution = "Null (avg/iteration)", stringsAsFactors = FALSE))
    df$distribution <- factor(df$distribution,
                              levels = c("Actual", "Null (avg/iteration)"))

    # y 轴变换：log10 / log2 用 log(count+1) 保 count=0 → y=0 基线；identity 直接用计数
    if (y_scale == "identity") {
      df$y      <- df$count
      y_breaks  <- ggplot2::waiver()
      y_lab     <- "Number of ASVs  [null ÷ iterations]"
    } else if (y_scale == "log2") {
      df$y      <- log2(df$count + 1)
      v         <- c(0, 2^(0:12))             # 0,1,2,4,…,4096
      v         <- v[v <= max(df$count)]
      y_breaks  <- log2(v + 1)
      y_labels  <- as.character(v)
      y_lab     <- "log2(Number of ASVs + 1)  [null ÷ iterations]"
    } else {                                  # "log10"
      df$y      <- log10(df$count + 1)
      v         <- c(0, 10^(0:5))             # 0,1,10,100,1000,10000
      v         <- v[v <= max(df$count)]
      y_breaks  <- log10(v + 1)
      y_labels  <- as.character(v)
      y_lab     <- "log10(Number of ASVs + 1)  [null ÷ iterations]"
    }
    if (y_scale == "identity") y_labels <- ggplot2::waiver()
    n_gt0 <- sum(act[2:length(act)])
    ggplot2::ggplot(df, ggplot2::aes(x = n_hosts, y = y, fill = distribution)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7, alpha = 0.85,
                        color = "black", linewidth = 0.1) +
      ggplot2::geom_vline(xintercept = thr, linetype = "dashed",
                          color = "firebrick", linewidth = 0.6) +
      ggplot2::scale_fill_manual(
        values = c("Actual" = color, "Null (avg/iteration)" = "gray65")) +
      ggplot2::scale_x_continuous(limits = c(0, x_max),
                                  breaks = seq(0, x_max, by = max(1, ceiling(x_max / 10)))) +
      ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels) +
      ggplot2::labs(
        title = sprintf("%s seed microbiota (n_hosts)", toupper(res$ds)),
        subtitle = sprintf("n = %d ASVs (n_hosts ≥ 1); %d iterations; generalist threshold = %d (dashed)",
                           n_gt0, nit, thr),
        x = "n_hosts (number of host species)",
        y = y_lab,
        fill = "") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     legend.position  = "bottom")
  }
  p1 <- panel(res_16s, thr_16s, "#4E79A7")
  p2 <- panel(res_its,  thr_its,  "#81B29A")
  p  <- patchwork::wrap_plots(p1, p2, ncol = 2)

  ggplot2::ggsave(out_pdf, p, width = 16, height = 7, device = grDevices::cairo_pdf)
  ggplot2::ggsave(out_jpg, p, width = 16, height = 7, dpi = 300)
  cat(sprintf("[gos-plot-nhost] n_hosts 分布对比图（y=%s, x≤%d）已保存：%s\n",
              y_scale, x_max, out_pdf))
  invisible(p)
}

# ---- Step 3：z 值 + 显著性 ----
# 每-n_hosts MPD 零模型：对每个出现的 n_hosts = k（≥2），随机抽 k 个宿主物种，
# 计算 MPD = 宿主两两遗传距离均值 = sum(子距离矩阵)/(k²−k)，N 次抽样得
# null_mu[k] / null_sd[k]。宿主集合 = 合并矩阵列（已与树取交）。
# 高效实现：指示矩阵 X（N×n_spp）→ xᵀDx = rowSums((X %*% D) * X)，MPD = xᵀDx/(k(k−1))。
gos_mpd_null <- function(merged_rds, host_dist_rds, out_rds, n_rep = 1000,
                         seed = 20260813, force = FALSE) {
  if (file.exists(out_rds) && !force) {
    cat(sprintf("[gos-mpd] 检查点命中：%s 已存在，跳过零模型\n", basename(out_rds)))
    return(readRDS(out_rds))
  }
  mg  <- readRDS(merged_rds)
  dmat_all <- readRDS(host_dist_rds)$dist
  spp <- intersect(colnames(mg$merged), colnames(dmat_all))
  dmat <- dmat_all[spp, spp, drop = FALSE]
  n_spp <- length(spp)
  nhosts <- rowSums(mg$merged)
  ks <- sort(unique(nhosts[nhosts >= 2]))      # 所有出现的 n_hosts ≥ 2
  cat(sprintf("[gos-mpd:%s] 宿主集合 %d 个物种；为 %d 个不同 n_hosts 值构建零模型\n",
              mg$ds, n_spp, length(ks)))

  set.seed(seed)
  null_mu <- setNames(rep(NA_real_, length(ks)), as.character(ks))
  null_sd <- setNames(rep(NA_real_, length(ks)), as.character(ks))
  pb <- utils::txtProgressBar(min = 0, max = length(ks), style = 3)
  for (j in seq_along(ks)) {
    k <- ks[j]
    S <- t(vapply(seq_len(n_rep), function(m) sample.int(n_spp, k), integer(k)))
    # 指示矩阵 X：N×n_spp，每行一个随机 k-子集
    X <- matrix(0L, n_rep, n_spp)
    X[cbind(rep(seq_len(n_rep), each = k), as.vector(t(S)))] <- 1L
    xDx <- rowSums((X %*% dmat) * X)           # = 2·sum_{a<b} d[S_a,S_b]
    mpd <- xDx / (k * (k - 1))
    null_mu[[j]] <- mean(mpd)
    null_sd[[j]] <- stats::sd(mpd)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)

  res <- list(ds = mg$ds, ks = ks, null_mu = null_mu, null_sd = null_sd,
              n_rep = n_rep, spp = spp, n_spp = n_spp)
  saveRDS(res, out_rds)
  cat(sprintf("[gos-mpd:%s] 完成：%d 个 k 值 × %d 次抽样零模型已保存 → %s\n",
              mg$ds, length(ks), n_rep, out_rds))
  invisible(res)
}

# 每 ASV z 值 + 显著性 + 完整信息表
#   输出表列：ASV_ID, n_hosts, mean_dist, null_mu, null_sd, z,
#            p_disp, p_clus, q_disp, q_clus
#   z = (mean_dist − null_mu[k])/null_sd[k]；仅 n_hosts ≥ 2 有定义（否则 NA）。
#   p_disp = P(Z ≥ z)（分散），p_clus = P(Z ≤ z)（聚集）；
#   q_disp / q_clus = 各自 BH-FDR 校正（在 n_hosts ≥ 2 子集内）。
gos_z_scores <- function(merged_rds, host_dist_rds, mpd_null_rds,
                         out_csv, ds = "") {
  mg   <- readRDS(merged_rds)
  dmat_all <- readRDS(host_dist_rds)$dist
  mpd  <- readRDS(mpd_null_rds)
  bin  <- mg$merged
  spp  <- intersect(colnames(bin), colnames(dmat_all))
  dmat <- dmat_all[spp, spp, drop = FALSE]
  bin  <- bin[, spp, drop = FALSE]

  n_asv <- nrow(bin)
  n_hosts <- rowSums(bin)
  mean_dist <- numeric(n_asv)
  for (i in seq_len(n_asv)) {
    ps <- colnames(bin)[bin[i, ] == 1L]
    n  <- length(ps)
    if (n >= 2) {
      sub <- dmat[ps, ps, drop = FALSE]
      mean_dist[i] <- sum(sub) / (n^2 - n)     # mean_dist：宿主两两遗传距离均值
    }
  }

  has <- n_hosts >= 2
  kvec <- as.character(n_hosts)
  null_mu <- rep(NA_real_, n_asv)
  null_sd <- rep(NA_real_, n_asv)
  if (sum(has) > 0) {
    null_mu[has] <- mpd$null_mu[kvec[has]]
    null_sd[has] <- mpd$null_sd[kvec[has]]
  }
  z <- ifelse(has, (mean_dist - null_mu) / null_sd, NA_real_)
  p_disp <- ifelse(has, stats::pnorm(z, lower.tail = FALSE), NA_real_)
  p_clus <- ifelse(has, stats::pnorm(z, lower.tail = TRUE),  NA_real_)
  q_disp <- q_clus <- rep(NA_real_, n_asv)
  if (sum(has) > 0) {
    q_disp[has] <- stats::p.adjust(p_disp[has], method = "BH")
    q_clus[has] <- stats::p.adjust(p_clus[has], method = "BH")
  }

  out <- data.frame(ASV_ID = rownames(bin), n_hosts = n_hosts,
                    mean_dist = mean_dist, null_mu = null_mu, null_sd = null_sd,
                    z = z, p_disp = p_disp, p_clus = p_clus,
                    q_disp = q_disp, q_clus = q_clus, stringsAsFactors = FALSE)
  data.table::fwrite(out, out_csv, sep = "\t")
  cat(sprintf(paste0("[gos-z:%s] 每 ASV 信息表已保存：%d ASV（n_hosts≥2 有 z 值 %d）→ %s\n"),
              ds, n_asv, sum(has), out_csv))
  invisible(out)
}

# ---- Step 4：分类 + 提取 ----
# 分类：通才 = n_hosts ≥ 阈值（显著聚集的通才仍归通才，聚集性只记录在 dispersion）；
#        专才 = n_hosts == 1；中间型 = 2 ≤ n_hosts < 阈值；
#        无代表性宿主 = n_hosts == 0（NRH，No representative host）。
# 附加 dispersion 属性列（n_hosts ≥ 2 时定义）：
#        显著分散 q_disp < alpha → "dispersed"（宿主系统发育显著分散，作为标注属性）；
#        显著聚集 q_clus < alpha → "clustered"；均不显著 → "neutral"。
gos_classify <- function(full_table, threshold, alpha = 0.05) {
  x <- as.data.frame(full_table)
  cat <- rep("intermediate", nrow(x))
  cat[x$n_hosts == 0] <- "NRH"
  cat[x$n_hosts == 1] <- "specialist"
  cat[x$n_hosts >= threshold] <- "generalist"
  disp <- rep(NA_character_, nrow(x))
  has2 <- !is.na(x$q_disp) & x$n_hosts >= 2
  disp[has2 & x$q_disp < alpha]                        <- "dispersed"
  disp[has2 & x$q_clus < alpha]                        <- "clustered"
  disp[has2 & x$q_disp >= alpha & x$q_clus >= alpha]   <- "neutral"
  x$dispersion <- factor(disp, levels = c("dispersed", "clustered", "neutral"))
  x$category <- factor(cat, levels = c("generalist", "specialist", "intermediate", "NRH"))
  x
}

# 按 ASV ID 集合从参考 fasta 中提取序列（手动解析：跨行序列整块保留）
gos_extract_fasta <- function(asv_ids, ref_fa, out_fa) {
  id_set <- unique(as.character(asv_ids))
  lines <- readLines(ref_fa, warn = FALSE)
  hpos <- which(grepl("^>", lines))
  heads <- sub("^>", "", lines[hpos])
  keep <- heads %in% id_set
  if (sum(keep) == 0) stop("[gos-extract-fa] 没有匹配到参考序列")
  keep_pos <- which(keep)
  out_lines <- character()
  for (i in keep_pos) {
    start <- hpos[i]
    end <- if (i < length(hpos)) hpos[i + 1] - 1L else length(lines)
    out_lines <- c(out_lines, lines[start:end])
  }
  writeLines(out_lines, out_fa)
  cat(sprintf("[gos-extract-fa] 已提取 %d 条参考序列 → %s\n", sum(keep), out_fa))
  invisible(out_fa)
}

# 按 ASV ID 集合从物种注释表提取分类注释
gos_extract_tax <- function(asv_ids, tax_file, out_txt) {
  id_set <- unique(as.character(asv_ids))
  tax <- data.table::fread(tax_file, sep = "\t", header = TRUE, fill = TRUE)  # 部分行 Sub_Species 缺失 → 填充 NA
  keep <- tax$taxid %in% id_set
  if (sum(keep) == 0) stop("[gos-extract-tax] 没有匹配到分类注释")
  data.table::fwrite(tax[keep, ], out_txt, sep = "\t")
  cat(sprintf("[gos-extract-tax] 已提取 %d 条分类注释 → %s\n", sum(keep), out_txt))
  invisible(out_txt)
}

# ============================================================================
# ---- Step 5：n_host - 相对丰度曲线图 ----
# 每个 Phylum（16S）/Class（ITS）一个子图：x = n_hosts，y = log10(相对丰度)，
# 每个点 = 一个 ASV 在一个宿主 species 中的相对丰度（丰度=0 不展示、不计数）。
#   丰度表先按宿主 species 聚合：同 species 的样本列求平均，再转相对丰度
#   （每列 / 列和），最大点数由 ASV×样本 降到 ASV×宿主 species。
#   点色：按相对丰度（高=红、低=蓝，log10 渐变）；alpha=0.25。
#   点形：generalist=三角形、specialist=正方形、其余=圆形。
#   子图外壳 = 上方浅灰背景 Phylum/Class 名称框（banner，不在图内）+ 下方散点图；
#   散点右下标注：ASV 数 / 点数 n / 平均 n_host / 正相关显著性
#           （Spearman 单侧 greater，n_hosts × log10 相对丰度）。
#   y 轴上界截断为 0（相对丰度不可能 > 100%）。
# 数据来源：CSS 丰度表 + metadata（样本→宿主 species）+ gos_classified_* + taxonomy。
# 仅统计 n_hosts >= min_nhost 的 ASV（默认 1）。拼图为 3 行 × 4 列，
# 顺序 = 平均 n_host 降序。
# 拼接全部在 R 中完成（gridExtra，点数经 species 聚合后已足够小，无需 JPG 中转）：
#   组合大图的面板不带图例，颜色（log10 相对丰度，红高蓝低）与形状
#   （generalist△/specialist□/other○）图例共享，统一放在大图底部；
#   颜色标尺上下限在前 top 组内全局固定 → 跨子图颜色含义一致。
# ============================================================================

#' 准备 n_host - 相对丰度长表数据（先按宿主 species 聚合，再转相对丰度）
#'   css 丰度表 → 按宿主 species 对样本列求平均 → 每列除以列和 → 长表。
#'   metadata$SampleID 需与 css 样本列名一致，宿主分组取 species_col 列。
gos_nhost_abund_data <- function(css_file, classified_file, tax_file, meta_file,
                                 group_col, min_nhost = 1, species_col = "species") {
  # 1) CSS 丰度表 + 样本列
  css <- data.table::fread(css_file, sep = "\t", header = TRUE)
  sample_cols <- setdiff(names(css), "ASV_ID")

  # 2) 宿主 species 映射（metadata$SampleID == css 样本列名）
  met <- data.table::fread(meta_file, sep = "\t", header = TRUE)
  met <- met[SampleID %in% sample_cols, c("SampleID", species_col), with = FALSE]
  data.table::setnames(met, c("SampleID", species_col), c("SampleID", "sp"))
  met <- met[!is.na(sp) & sp != "" & sp != "NA"]
  if (nrow(met) == 0) {
    stop("[gos-nhost-abund] metadata 与 css 样本列无交集，无法按宿主 species 聚合")
  }
  sp_names <- unique(met$sp)
  sp2cols  <- split(met$SampleID, met$sp)

  # 3) 分类表：仅 n_hosts >= min_nhost 的 ASV
  cl <- data.table::fread(classified_file, sep = "\t")
  cl <- cl[n_hosts >= min_nhost, .(ASV_ID, n_hosts, category)]

  # 4) 分类注释（fill = TRUE 容忍 Sub_Species 缺失的行）
  tax <- data.table::fread(tax_file, sep = "\t", header = TRUE, fill = TRUE,
                           select = c("taxid", group_col))
  data.table::setnames(tax, c("taxid", group_col), c("ASV_ID", "grp"))

  # 5) 子集 CSS 行 → 数值矩阵，按 species 聚合（同 species 样本列求平均）
  keep <- css$ASV_ID %in% cl$ASV_ID
  rn <- css$ASV_ID[keep]
  m <- as.matrix(css[keep, sample_cols, with = FALSE])
  rm(css); invisible(gc())
  agg <- lapply(sp_names, function(sp) {
    cols <- sp2cols[[sp]]
    if (length(cols) == 1L) m[, cols] else rowMeans(m[, cols, drop = FALSE])
  })
  A <- do.call(cbind, agg)
  colnames(A) <- sp_names
  cs <- colSums(A)          # species 列和 → 相对丰度分母
  rm(m, agg); invisible(gc())

  # 6) 仅取非零格点 → 相对丰度 → 长表
  idx <- which(A > 0, arr.ind = TRUE)
  long <- data.table::data.table(
    ASV_ID = rn[idx[, 1L]],
    rel    = A[idx] / cs[idx[, 2L]]
  )
  rm(A); invisible(gc())
  long[, log10_rel := log10(rel)]

  # 7) 拼接 n_hosts / category / 分类层级
  long[cl,  on = "ASV_ID", c("n_hosts", "category") := .(i.n_hosts, i.category)]
  long[tax, on = "ASV_ID", grp := i.grp]
  long[]
}

#' 按 Phylum/Class 内 ASV 数量降序取前 top 个分类层级（用于展示分组）
#'   排序键 = 每组包含的 ASV 数（n_asv 降序），mean_nhost 仅作并列时的次级排序。
gos_nhost_abund_rank <- function(dt, top = 12L, drop_unassigned = TRUE) {
  asv <- unique(dt[!is.na(grp), .(ASV_ID, n_hosts, grp)])
  if (drop_unassigned) asv <- asv[!grp %in% c("Unassigned", "", "NA")]
  stat <- asv[, .(n_asv = .N, mean_nhost = mean(n_hosts)), by = grp]
  data.table::setorder(stat, -n_asv, -mean_nhost)
  if (nrow(stat) < top) {
    cat(sprintf("[gos-nhost-abund] 警告：只有 %d 个分类层级 < top=%d，全部使用\n",
                nrow(stat), top))
    top <- nrow(stat)
  }
  stat[seq_len(top)]
}

#' 单个子图：n_hosts (x) × log10 相对丰度 (y)。
#' 叠加线性拟合曲线（y = a + b·x）与 95% 置信带（灰色填充 + 虚线上下界），
#' 右下指标框追加拟合方程与 R²。
#' 名称框不在此函数内——由 gos_nhost_abund_banner 加在散点图上方。
#' y 轴上界截断为 0（相对丰度不可能 > 100%，log10(1)=0）。
#' color_limits：共享颜色标尺上下限（log10 相对丰度，默认为 NULL=各子图自适应）。
#' legend：图例位置（"right"/"bottom"/"none"），组合大图的面板用 "none"。
#' 注：颜色直接用 log10_rel（线性渐变），不再用 rel + trans="log10"（等价但更直观，
#'     且 limits/breaks 无需考虑 trans 空间）。
gos_nhost_abund_subplot <- function(dt_grp, grp_name, stat_row, alpha = 0.25,
                                    color_limits = NULL, y_limits = NULL,
                                    legend = "right",
                                    show_axis_titles = TRUE, metric_size = 3.2) {
  n_asv  <- stat_row$n_asv
  n_pts  <- nrow(dt_grp)
  mean_nh <- stat_row$mean_nhost

  p_pos <- tryCatch(
    suppressWarnings(stats::cor.test(dt_grp$n_hosts, dt_grp$log10_rel,
                                    method = "spearman",
                                    alternative = "greater")$p.value),
    error = function(e) NA_real_)

  xr   <- range(dt_grp$n_hosts)
  # y 下界提前算好：也用于把置信带钳制到可见范围；
  # y_limits 非空时使用跨面板共享的 y 轴范围（两行图共享 y 轴）
  y_lo <- if (!is.null(y_limits)) y_limits[1] else min(dt_grp$log10_rel, -0.5)

  # 线性拟合（点少或 x 无变异时不拟合）+ 置信带数据（95% 预测区间）
  has_fit <- n_pts >= 3L && length(unique(dt_grp$n_hosts)) >= 2L
  fit_lab <- character(0)
  band <- NULL
  if (has_fit) {
    fit <- tryCatch(stats::lm(log10_rel ~ n_hosts, data = dt_grp),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      cf <- stats::coef(fit)
      r2 <- summary(fit)$r.squared
      eq <- sprintf("y = %sx %s %.3f",
                    formatC(cf[2], digits = 3, format = "g"),
                    ifelse(cf[1] < 0, "-", "+"),
                    abs(cf[1]))
      fit_lab <- c(eq, sprintf("R² = %.3f", r2))
      # 置信带：在 x 范围内均匀取点 → predict ± t*SE
      xseq <- seq(min(dt_grp$n_hosts), max(dt_grp$n_hosts), length.out = 200)
      pr <- stats::predict(fit, newdata = data.frame(n_hosts = xseq), se.fit = TRUE)
      crit <- stats::qt(0.975, df = fit$df.residual)
      band <- data.frame(x = xseq, yfit = pr$fit,
                         lo = pr$fit - crit * pr$se.fit,
                         hi = pr$fit + crit * pr$se.fit)
      # 钳制到 [y_lo, 0]：退化/微小面板的整段带可能落在 y 轴范围外，
      # 不钳制会被 ggplot 整行删除 → 置信带完全消失
      band$lo   <- pmax(band$lo, y_lo)
      band$hi   <- pmin(band$hi, 0)
      band$yfit <- pmin(pmax(band$yfit, y_lo), 0)
    }
  }

  dt_grp$category <- factor(dt_grp$category,
                            levels = c("generalist", "specialist", "intermediate"))

  # P_pos 用星号标注显著性：*** p<0.001，** p<0.01，* p<0.05，ns 不显著
  p_star <- if (is.na(p_pos)) "NA"
            else if (p_pos < 0.001) "***"
            else if (p_pos < 0.01)  "**"
            else if (p_pos < 0.05)  "*"
            else "ns"
  metric_lab <- paste(c(
    sprintf("ASVs: %d", n_asv),
    sprintf("n: %d", n_pts),
    sprintf("mean n_host: %.2f", mean_nh),
    sprintf("P_pos: %s", p_star),
    fit_lab
  ), collapse = "\n")

  # 坐标轴 title：组合大图的面板里隐去，由组图统一标注（show_axis_titles = FALSE）
  axis_labs <- if (show_axis_titles)
    ggplot2::labs(x = "n_hosts", y = "log10(relative abundance)")
  else
    ggplot2::labs(x = NULL, y = NULL)

  p <- ggplot2::ggplot(dt_grp,
                       ggplot2::aes(x = n_hosts, y = log10_rel,
                                    color = log10_rel, shape = category)) +
    ggplot2::geom_point(size = 1.2, alpha = alpha)
  # 丰度低/中/高三段配色（共享标尺时中色落在 limits 的中点）。
  # 压缩黄色/橙色的占比：蓝、红各保留一段纯色平台，中间过渡带收窄。
  abund_cols <- c("#93BCD0", "#93BCD0", "#FCB366", "#FE816F", "#FE816F")
  abund_vals <- c(0, 0.35, 0.5, 0.65, 1)
  if (!is.null(color_limits)) {
    p <- p + ggplot2::scale_color_gradientn(
      colours = abund_cols, values = abund_vals,
      limits = color_limits,
      breaks = pretty(color_limits, n = 5),
      name = "log10\nrel. abund.")
  } else {
    p <- p + ggplot2::scale_color_gradientn(colours = abund_cols,
                                            values = abund_vals,
                                            name = "log10\nrel. abund.")
  }
  p <- p + ggplot2::scale_shape_manual(
    values = c(generalist = 17, specialist = 15, intermediate = 16),
    labels = c(generalist = "Generalist", specialist = "Specialist",
               intermediate = "Other"),
    name = NULL)
  if (!is.null(band)) {
    p <- p +
      ggplot2::geom_ribbon(data = band, inherit.aes = FALSE,
                           ggplot2::aes(x = x, ymin = lo, ymax = hi),
                           fill = "grey60", alpha = 0.35) +
      ggplot2::geom_line(data = band, inherit.aes = FALSE,
                         ggplot2::aes(x = x, y = yfit),
                         color = "black", linewidth = 0.6) +
      ggplot2::geom_line(data = band, inherit.aes = FALSE,
                         ggplot2::aes(x = x, y = lo),
                         color = "grey30", linewidth = 0.3, linetype = "dashed") +
      ggplot2::geom_line(data = band, inherit.aes = FALSE,
                         ggplot2::aes(x = x, y = hi),
                         color = "grey30", linewidth = 0.3, linetype = "dashed")
  }
  p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(
      limits = if (!is.null(y_limits)) y_limits else c(y_lo, 0),
      expand = ggplot2::expansion(mult = c(0.06, 0))) +
    ggplot2::annotate("label", x = xr[2], y = y_lo,
                      label = metric_lab, size = metric_size, hjust = 1, vjust = 0,
                      fill = grDevices::rgb(1, 1, 1, 0.85),
                      linewidth = 0.25,
                      label.padding = grid::unit(0.35, "lines")) +
    axis_labs +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend,
      legend.box = if (legend == "bottom") "horizontal" else "vertical",
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width  = grid::unit(0.35, "cm"),
      # 字号整体放大 50%（2026-11-11；右下角指标框单独放大 25%，由调用方给出）
      legend.title = ggplot2::element_text(size = 10.5),
      legend.text  = ggplot2::element_text(size = 9.75),
      axis.title   = ggplot2::element_text(size = 13.5),
      axis.text    = ggplot2::element_text(size = 11.25)
    )
}

#' 子图外壳：上方浅灰背景名称框 + 下方散点图（banner 框在最上方，不在图内）
gos_nhost_abund_banner <- function(p, grp_name, banner_h = 0.6, scatter_h = 5,
                                   banner_size = 7.2) {
  banner <- ggplot2::ggplot(data.frame(x = 0, y = 0, lab = grp_name),
                            ggplot2::aes(x = x, y = y, label = lab)) +
    ggplot2::geom_text(fontface = "bold", size = banner_size) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "grey90",
                                                          colour = "grey35",
                                                          linewidth = 0.4),
                   plot.margin = ggplot2::margin(1, 1, 1, 1))
  gridExtra::arrangeGrob(banner, ggplot2::ggplotGrob(p),
                         ncol = 1, heights = c(banner_h, scatter_h))
}

#' 从一张带图例的图中提取 guide-box（颜色标尺 + 形状图例一起），供组合大图共享。
#' 传入的子图需用 legend = "bottom"（水平排布，便于放组合图底部）。
gos_nhost_abund_legend <- function(p) {
  gt <- ggplot2::ggplotGrob(p)
  # ggplot2 4.0.x 图例布局名带位置后缀："guide-box-right"/"guide-box-bottom"/...，
  # 不再是旧的 "guide-box"，用前缀匹配
  idx <- which(grepl("^guide-box", gt$layout$name))
  if (length(idx) == 0L) {
    warning("[gos-nhost-abund] 图中未找到图例（guide-box），请检查 legend 位置参数")
    return(NULL)
  }
  # 只保留实际装有图例的项：剔除 NULL 和 zeroGrob 空占位
  #（ggplot2 4.0.x 会在 5 个位置各放一个 guide-box-* 行，真正的图例只装一个，
  #  其余是 zeroGrob 占位）
  idx <- idx[vapply(idx, function(k) {
    g <- gt$grobs[[k]]
    !is.null(g) && !inherits(g, "zeroGrob")
  }, logical(1))]
  if (length(idx) == 0L) {
    warning("[gos-nhost-abund] 图例 grob 为空，无法提取")
    return(NULL)
  }
  gt$grobs[[idx[1L]]]
}

#' 将子图 grob 拼为 nrow×ncol 大图，可选的共享图例放在底部；同时存 PDF + JPG。
gos_nhost_abund_combine_grobs <- function(grobs, out_base, legend_grob = NULL,
                                          nrow = 3L, ncol = 4L,
                                          panel_w = 5.2, panel_h = 4.6,
                                          dpi = 300) {
  n <- length(grobs)
  if (n > nrow * ncol) {
    stop(sprintf("[gos-nhost-abund] %d 个子图 > %d×%d 网格", n, nrow, ncol))
  }
  if (n < nrow * ncol) {          # 不足时用空白面板补齐
    blank <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
    grobs <- c(grobs, replicate(nrow * ncol - n, blank, simplify = FALSE))
  }
  panels <- do.call(gridExtra::arrangeGrob,
                    c(grobs, list(nrow = nrow, ncol = ncol)))
  H <- nrow * panel_h
  if (!is.null(legend_grob)) {
    full <- gridExtra::arrangeGrob(panels, legend_grob,
                                   ncol = 1, heights = c(H, 0.9))
    H <- H + 0.9
  } else {
    full <- panels
  }
  W <- ncol * panel_w
  ggplot2::ggsave(paste0(out_base, ".pdf"), full, width = W, height = H,
                  device = cairo_pdf, limitsize = FALSE)
  ggplot2::ggsave(paste0(out_base, ".jpg"), full, width = W, height = H,
                  dpi = dpi, device = "jpeg", limitsize = FALSE)
  cat(sprintf("[gos-nhost-abund] 已拼接 %d 个子图（%d×%d）→ %s.{pdf,jpg}\n",
              n, nrow, ncol, out_base))
  invisible(out_base)
}

#' 完整流程：数据 → 排序 → 逐子图 PDF+JPG → 3×4 拼接大图（共享图例）。
#'   子图文件各自带右侧图例（单独查看用）；组合大图面板不带图例，
#'   颜色/形状图例共享、统一放底部。颜色标尺上下限在 top 组内全局固定，
#'   保证跨子图颜色可比（红=高丰度、蓝=低丰度）。
gos_nhost_abund_plot <- function(css_file, classified_file, tax_file, meta_file,
                                 group_col, ds, out_dir, top = 12L, alpha = 0.25,
                                 drop_unassigned = TRUE,
                                 panel_w = 5.2, panel_h = 4.6,
                                 shared_legend = TRUE) {
  sub_dir <- file.path(out_dir, paste0("subplot_", ds))
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  dt   <- gos_nhost_abund_data(css_file, classified_file, tax_file, meta_file, group_col)
  stat <- gos_nhost_abund_rank(dt, top = top, drop_unassigned = drop_unassigned)
  cat(sprintf("[gos-nhost-abund:%s] 前 %d 个 %s（ASV 数量降序）：\n",
              ds, nrow(stat), group_col))
  print(stat[, .(grp, mean_nhost, n_asv)])

  # 共享颜色标尺：top 组内 log10 相对丰度全局范围（上限固定 0）
  color_limits <- c(min(dt[grp %in% stat$grp]$log10_rel), 0)

  # 1) 逐个生成子图（自带右侧图例），保存 PDF+JPG
  for (i in seq_len(nrow(stat))) {
    g  <- stat$grp[i]
    dg <- dt[grp == g]
    p  <- gos_nhost_abund_subplot(dg, g, stat[i], alpha = alpha,
                                  color_limits = color_limits)
    grob <- gos_nhost_abund_banner(p, g)
    tag <- sprintf("%02d_%s", i, gsub("[^A-Za-z0-9._-]", "_", g))
    base <- file.path(sub_dir, sprintf("%s_%s", ds, tag))
    ggplot2::ggsave(paste0(base, ".pdf"), grob, width = panel_w, height = panel_h,
                    device = cairo_pdf, limitsize = FALSE)
    ggplot2::ggsave(paste0(base, ".jpg"), grob, width = panel_w, height = panel_h,
                    dpi = 300, device = "jpeg", limitsize = FALSE)
    cat(sprintf("  已完成 %d/%d: %s (ASV=%d, 点=%d)\n",
                i, nrow(stat), g, stat$n_asv[i], nrow(dg)))
  }

  # 2) 共享图例：从第一张子图（legend="bottom" 水平排布）提取 guide-box
  legend_grob <- NULL
  if (shared_legend) {
    p_leg <- gos_nhost_abund_subplot(dt[grp == stat$grp[1L]], stat$grp[1L],
                                     stat[1L], alpha = alpha,
                                     color_limits = color_limits, legend = "bottom")
    legend_grob <- gos_nhost_abund_legend(p_leg)
  }

  # 3) 拼接 3×4：面板不带图例，颜色/形状图例统一放底部
  panel_grobs <- lapply(seq_len(nrow(stat)), function(i) {
    dg <- dt[grp == stat$grp[i]]
    p  <- gos_nhost_abund_subplot(dg, stat$grp[i], stat[i], alpha = alpha,
                                  color_limits = color_limits, legend = "none")
    gos_nhost_abund_banner(p, stat$grp[i])
  })
  out_base <- file.path(out_dir, sprintf("nhost_abund_curves_%s", ds))
  gos_nhost_abund_combine_grobs(panel_grobs, out_base, legend_grob,
                                nrow = 3L, ncol = 4L,
                                panel_w = panel_w, panel_h = panel_h)
  invisible(list(dt = dt, stat = stat, color_limits = color_limits))
}

# ============================================================================
# ---- Step 6：通才/专才分类群组成差异（Genus~Class 双向柱状图） ----
# 目的：§8.4 已把每个 ASV 分为通才/专才，本函数量化两类在 Genus~Class 各级的
#   组成差异：① 每个分类群在通才/专才中的 ASV 占比（占各自总 ASV 数 %）；
#   ② Fisher 精确检验（通才 vs 专才 2×2）+ BH-FDR → 显著富集方向
#      （三级星号标在 OR 值右侧：* FDR<0.05、** FDR<0.01、*** FDR<0.001）。
#   图 = 双向（diverging）柱状图：x>0 右橙 = 通才占比，x<0 左蓝 = 专才占比；
#   y 轴按 OR 升序（ggplot 底部 = level 1）→ 顶部为 OR 最大（最通才富集）、
#   底部为 OR 最小（仅专才出现 → OR=0 沉底）；通才柱末端标注该分类群的 OR
#   + 三级星号；单侧检出分类群 OR 无意义，退化为类别徽标：仅通才 → OG（橙柱端）、
#   仅专才 → OS（蓝柱左端），显著时 OG/OS 后附星号（FDR 阈值同 OR）。
# 输入：每数据集一张通才表 + 一张专才表（每行 = 一个 ASV，含 Kingdom~Genus 分类列）。
# 输出：组成统计表（全分类群×层级 + fold_GS/OR/p/FDR）+ coverage 表 + 2×2 拼图
#   （行 = 数据集 16S/ITS，列 = Genus/Class）。
# 口径：
#   - "出现频率" = 该分类群包含的 ASV 数（通才+专才之和），与 Step 5 的
#     gos_nhost_abund_rank（ASV 数降序）口径一致。
#   - 排除无效注释：Genus 层 "Unassigned"/"Incertae Sedis"/空值 不进 Top-n 排名，
#     但占比单独记为 coverage 供报告。
#   - 配色沿用 §6.1 PAL1：Generalist=橙 #D35400，Specialist=蓝 #547FB6。
# ============================================================================
gos_composition <- function(gen_tax_files, spe_tax_files, datasets,
                            out_table, out_coverage, out_fig_base,
                            top_class = 12L, top_genus = 20L,
                            pal_gen = "#D35400", pal_spe = "#547FB6",
                            ranks = c("Class", "Order", "Family", "Genus"),
                            plot_ranks = c("Class", "Genus")) {
  # gen_tax_files / spe_tax_files : 命名列表，键 = 数据集键（如 16s/its），
  #                                 值 = 对应通才/专才分类注释表路径
  # datasets                      : 与 gen_tax_files 同键的命名字符向量，值 = 数据集显示名
  # out_table / out_coverage      : 组成统计表 / 无效注释 coverage 表输出路径
  # out_fig_base                  : 拼图输出基名（不含扩展名，自动 .pdf + .jpg）
  # top_class / top_genus         : 各层级 Top-n（通才∪专才各自 Top-n 并集去重，不保证总数）
  # pal_gen / pal_spe             : §6.1 配色（通才橙 / 专才蓝）
  # ranks                         : 表内计算的分类层级（默认 4 层）
  # plot_ranks                    : 参与拼图的层级（默认只画 Class 与 Genus，减少图数）
  stopifnot(length(gen_tax_files) == length(spe_tax_files),
            length(gen_tax_files) == length(datasets),
            all(names(gen_tax_files) %in% names(datasets)))

  drop <- c("Unassigned", "Incertae Sedis", "", "NA", "incertae_sedis")

  # ---- 单数据集 × 全层级：组成 + Fisher ----
  comp_one <- function(gen_file, spe_file, ds) {
    gen <- data.table::fread(gen_file, sep = "\t", header = TRUE, fill = TRUE)
    spe <- data.table::fread(spe_file, sep = "\t", header = TRUE, fill = TRUE)
    Ng <- nrow(gen); Ns <- nrow(spe)
    out <- lapply(ranks, function(rk) {
      g <- gen[[rk]]; s <- spe[[rk]]
      tax <- sort(unique(c(g, s)))
      valid <- tax[!is.na(tax) & !(tax %in% drop)]
      ng <- as.integer(table(factor(g, levels = valid)))
      ns <- as.integer(table(factor(s, levels = valid)))
      total <- ng + ns
      pct_g <- 100 * ng / Ng
      pct_s <- 100 * ns / Ns
      or <- p <- rep(NA_real_, length(valid))
      for (i in seq_along(valid)) {
        ft <- stats::fisher.test(matrix(c(ng[i], Ng - ng[i], ns[i], Ns - ns[i]), nrow = 2))
        or[i] <- as.numeric(ft$estimate)
        p[i]  <- ft$p.value
      }
      data.table::data.table(dataset = ds, rank = rk, taxon = valid,
                             n_generalist = ng, n_specialist = ns, n_total = total,
                             pct_generalist = pct_g, pct_specialist = pct_s,
                             fold_GS = ifelse(pct_s > 0, pct_g / pct_s, Inf),
                             OR = or, p_fisher = p, FDR = stats::p.adjust(p, method = "BH"))
    })
    data.table::rbindlist(out)
  }

  comp_all <- data.table::rbindlist(
    lapply(names(gen_tax_files), function(k)
      comp_one(gen_tax_files[[k]], spe_tax_files[[k]], datasets[[k]])))
  data.table::setorder(comp_all, dataset, rank, -n_total)
  data.table::fwrite(comp_all, out_table, sep = "\t")
  cat(sprintf("[gos-composition] 组成统计表已保存（%d 行）→ %s\n", nrow(comp_all), out_table))

  # ---- coverage：无效注释占比（供报告，不进排名）----
  cover_line <- function(gen_file, spe_file, ds) {
    gen <- data.table::fread(gen_file, sep = "\t", header = TRUE, fill = TRUE)
    spe <- data.table::fread(spe_file, sep = "\t", header = TRUE, fill = TRUE)
    frac <- function(x) 100 * sum(!is.na(x) & (x %in% drop)) / length(x)
    data.frame(dataset = ds,
               generalist_n = nrow(gen), specialist_n = nrow(spe),
               genus_unassigned_pct_generalist = frac(gen$Genus),
               genus_unassigned_pct_specialist = frac(spe$Genus),
               class_unassigned_pct_generalist = frac(gen$Class),
               class_unassigned_pct_specialist = frac(spe$Class))
  }
  cover <- do.call(rbind, lapply(names(gen_tax_files), function(k)
    cover_line(gen_tax_files[[k]], spe_tax_files[[k]], datasets[[k]])))
  data.table::fwrite(cover, out_coverage, sep = "\t")
  cat("\n== annotation coverage (excluded from ranking, % of ASVs) ==\n")
  print(cover, row.names = FALSE)

  # ---- 双向柱状图（单面板）----
  #  x>0 = 通才占比（橙），x<0 = 专才占比（蓝）；通才柱端标注 OR + 三级星号，
  #  单侧检出分类群改标 OG/OS（仅通才/仅专才）。
  #  y 轴按 OR 排序：ggplot 离散 y 轴级别 1 在底部，故 OR 升序 →
  #  最大 OR（最通才富集）落在顶部、最小 OR（最专才富集）落在底部。
  make_div <- function(ds, rk, top_n) {
    d <- comp_all[dataset == ds & rank == rk]
    # 选取：通才按组内 ASV 数取 Top-n、专才按组内 ASV 数取 Top-n，
    # 两者并集去重后全画（不保证总数）；同一分类群两侧同画。
    # 只计入该组内出现过的分类群（另一侧为 0 不影响入选）。
    gen_top <- d[n_generalist > 0][order(-n_generalist)][seq_len(min(top_n, .N))]
    spe_top <- d[n_specialist > 0][order(-n_specialist)][seq_len(min(top_n, .N))]
    d <- unique(data.table::rbindlist(list(gen_top, spe_top)))
    # 排序：OR（Fisher 通才 vs 专才 2×2 优势比），升序排列后 ggplot 底部 = level 1，
    # 顶部 = OR 最大（最通才富集）；仅通才出现 → OR=Inf 自然置顶；
    # 仅专才出现 → OR→0 沉底。NA 不应出现（入选 taxon 至少一侧有 ASV）。
    d <- d[order(OR)]
    d[, taxon := factor(taxon, levels = taxon)]
    # 三级显著性星号（基于 BH-FDR）：* FDR<0.05、** FDR<0.01、*** FDR<0.001。
    d[, star := data.table::fifelse(FDR < 0.001, "***",
              data.table::fifelse(FDR < 0.01, "**",
              data.table::fifelse(FDR < 0.05, "*", "")))]
    # 端部标注（单侧检出分类群 OR 无意义，退化为类别徽标 OG/OS）：
    #   lab_r（通才柱端/右侧）：双侧检出 → "OR(2 位小数) + 星号"；
    #     仅通才（n_specialist==0，OR=Inf）→ "OG"(显著时附星号)。
    #   lab_l（专才柱左端）：仅专才（n_generalist==0，OR→0）→ "OS"(显著时附星号)，
    #     该类无通才柱，故标在蓝柱外侧；双侧检出不在左侧另标（信号已由右侧 OR 表达）。
    d[, star_suf := ifelse(star == "", "", paste0(" ", star))]
    d[, lab_r := data.table::fifelse(
         n_generalist > 0L & n_specialist > 0L,
         data.table::fifelse(is.na(OR), "",
            data.table::fifelse(is.infinite(OR), "∞",
               sprintf("%.2f %s", OR, star))),
         data.table::fifelse(n_generalist > 0L, paste0("OG", star_suf),
                             NA_character_))]
    d[, lab_l := data.table::fifelse(n_specialist > 0L & n_generalist == 0L,
                                     paste0("OS", star_suf), NA_character_)]
    lab_r_dat <- d[pct_generalist > 0 & !is.na(lab_r)]
    lab_l_dat <- d[!is.na(lab_l)]

    long <- rbind(
      d[, .(taxon, group = "Generalist", value = pct_generalist)],
      d[, .(taxon, group = "Specialist", value = -pct_specialist)]
    )
    long[, group := factor(group, levels = c("Generalist", "Specialist"))]
    maxv <- max(d$pct_generalist, d$pct_specialist)

    ggplot2::ggplot(long, ggplot2::aes(x = value, y = taxon, fill = group)) +
      ggplot2::geom_col(width = 0.72, colour = NA) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.35) +
      ggplot2::geom_text(data = lab_r_dat,
                         ggplot2::aes(x = pct_generalist + 0.02 * maxv, y = taxon,
                                      label = lab_r),
                         hjust = 0, vjust = 0.5, size = 2.4, colour = "grey25",
                         inherit.aes = FALSE) +
      ggplot2::geom_text(data = lab_l_dat,
                         ggplot2::aes(x = -pct_specialist - 0.02 * maxv, y = taxon,
                                      label = lab_l),
                         hjust = 1, vjust = 0.5, size = 2.4, colour = "grey25",
                         inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = c(Generalist = pal_gen, Specialist = pal_spe),
                                 name = NULL,
                                 labels = c(Generalist = "Generalist", Specialist = "Specialist")) +
      ggplot2::scale_x_continuous(labels = function(x) sprintf("%.0f%%", abs(x)),
                                  expand = ggplot2::expansion(mult = c(0.05, 0.24))) +
      ggplot2::scale_y_discrete() +
      ggplot2::labs(x = "Proportion of ASVs (%)",
                    title = sprintf("%s  ·  %s", ds, rk)) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold", size = 10.5, hjust = 0.5),
                     axis.title.y = ggplot2::element_blank(),
                     legend.position = "bottom",
                     legend.key.size = grid::unit(0.4, "cm"),
                     legend.text = ggplot2::element_text(size = 9))
  }

  # ---- 2×2 拼图：行 = 数据集（16S 上 / ITS 下），列 = Genus/Class ----
  dsk <- names(datasets)
  p16g <- make_div(datasets[[dsk[1]]], "Genus", top_genus)
  p16c <- make_div(datasets[[dsk[1]]], "Class", top_class)
  pitg <- make_div(datasets[[dsk[2]]], "Genus", top_genus)
  pitc <- make_div(datasets[[dsk[2]]], "Class", top_class)

  combined <- (p16g | p16c) / (pitg | pitc) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Taxonomic composition of generalist vs specialist seed microbes",
      subtitle = paste0(
        "Orange = generalist (right), blue = specialist (left).  ",
        "Labels on orange bars = odds ratio (generalist / specialist);  ",
        "OG / OS = taxon detected only among generalists / specialists.  ",
        "*, **, *** = Fisher exact test, FDR < 0.05 / 0.01 / 0.001."),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
                             plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey35"),
                             legend.position = "bottom"))

  ggplot2::ggsave(paste0(out_fig_base, ".pdf"), combined, width = 12, height = 11,
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(out_fig_base, ".jpg"), combined, width = 12, height = 11, dpi = 300)
  cat(sprintf("[gos-composition] 双向柱状拼图已保存 → %s.{pdf,jpg}\n", out_fig_base))

  # ---- 关键结果摘要打印 ----
  cat("\n== top taxa by ASV frequency (Genus level) ==\n")
  for (k in dsk) {
    print(comp_all[dataset == datasets[[k]] & rank == "Genus"][order(-n_total)][1:10,
          .(taxon, n_generalist, n_specialist, pct_generalist, pct_specialist, FDR)], topn = 10)
  }

  cat("\nGOS COMPOSITION DONE. outputs:", out_table, ";", out_coverage, ";",
      paste0(out_fig_base, ".{pdf,jpg}"), "\n")
  invisible(combined)
}

# ============================================================================
# ---- Step 7：dispersed / clustered ASV 宿主分布热图 + clustered 注释表 ----
# 目的：§8.4 分类输出含 dispersion 属性列（clustered / dispersed / neutral）。
#   本函数用一张热图展示全部 clustered 与 dispersed ASV 在宿主之间的实际存在分布，
#   回答"被取消通才资格的聚集型 ASV 到底集中在哪些宿主"，并输出全部 clustered ASV
#   的物种注释 + 代表宿主属（Repres_genus / Repres_N）。
# 输入（全部参数化，调用方传地址）：
#   classified_files —— 命名列表（键 B/F）→ gos_classified_{16s,its}.txt
#   presence_files  —— 命名列表（键 B/F）→ presence_merged_{16s,its}.txt
#   tax_files       —— 命名列表（键 B/F）→ filtered_merge/{16S,ITS}/Merged_*_taxonomy_filtered.txt
#   meta_files      —— 命名列表（键 B/F）→ standard_merge/metadata_{16S,ITS}.txt
#   tree_file       —— 宿主物种树 .nwk（plant_species_tree_all.nwk）
#   color_map_rds   —— 扁平配色 RDS（含 $order 向量，宿主 Order → 颜色）
#   out_selected    —— 挑选/重命名映射表输出路径（ASV_ID ↔ new_name + grp/db/n_hosts）
#   out_taxonomy    —— clustered ASV 注释表输出路径（Kingdom~Subspecies + Repres_genus/N）
#   out_fig_base    —— 热图输出基名（自动 .pdf + .jpg）
# 口径：
#   - 挑选全部 dispersion ∈ {clustered, dispersed} 且 n_hosts ≥ 1 的 ASV（不按 n_hosts 对齐）。
#   - 重命名 {C|D}_{B|F}_{N}_ASV{seq}（C/D = clustered/dispersed，B/F = 16S/ITS，
#     N = n_hosts，seq 为 (grp, db, n_hosts) 组内顺序号）。
#   - 热图：176 宿主 species（行，按物种树系统发育序）× 全部选中 ASV（列）；
#     C/D×B/F 四组分别聚类（组内存在/不存在二元模式，complete），组间 split；
#     块 = 存在时按 log2(n_hosts) 赋色（高深低浅）、不存在/未测 = 白（NA）；
#     顶部色带 = C/D×B/F 4 色，右侧色带 = 宿主 Order（color_map_flat$order）；
#     列标签隐藏（ASV 数过多）。图内全英文。
#   - Repres_genus / Repres_N：该 ASV 存在的宿主 species 按宿主 Genus 计数，
#     出现次数最多的 Genus（并列 `;` 连接）/ 对应出现次数。
# ============================================================================
gos_dispersion_heatmap <- function(classified_files, presence_files, tax_files,
                                   meta_files, tree_file, color_map_rds,
                                   out_selected, out_taxonomy, out_fig_base,
                                   fig_width = 14, fig_height = 16.2,
                                   jpg_width = 2800, jpg_height = 3240) {
  stopifnot(setequal(names(classified_files), names(presence_files)),
            setequal(names(classified_files), names(tax_files)),
            setequal(names(classified_files), names(meta_files)),
            file.exists(tree_file), file.exists(color_map_rds))
  suppressMessages({ library(ComplexHeatmap); library(circlize) })

  # ---- 1. 挑选全部 clustered / dispersed ASV（n_hosts ≥ 1）----
  sel <- lapply(names(classified_files), function(k) {
    d <- data.table::fread(classified_files[[k]])
    d <- d[n_hosts >= 1 & dispersion %in% c("clustered", "dispersed")]
    d[, grp := ifelse(dispersion == "clustered", "C", "D")]
    d[, db := k]
    cat(sprintf("[gos-disp] %s clustered=%d  dispersed=%d  n_hosts %s\n",
                k, d[grp == "C", .N], d[grp == "D", .N],
                paste(range(d$n_hosts), collapse = "-")))
    d
  })
  names(sel) <- names(classified_files)
  all_sel <- data.table::rbindlist(sel, use.names = TRUE)
  cat("[gos-disp] TOTAL selected ASVs:", nrow(all_sel), "\n")
  print(all_sel[, .N, by = .(grp, db)])

  # ---- 2. 重命名 {C|D}_{B|F}_{N}_ASV{seq}（组内顺序号）----
  all_sel <- all_sel[order(grp, db, n_hosts, ASV_ID)]
  all_sel[, new_name := sprintf("%s_%s_%d_ASV%04d", grp, db, n_hosts, seq_len(.N)),
          by = .(grp, db, n_hosts)]
  id_map <- all_sel[, .(ASV_ID, new_name, grp, db, n_hosts, dispersion)]
  data.table::fwrite(id_map, out_selected, sep = "\t")
  cat("[gos-disp] rename map written:", nrow(id_map), "rows ->", out_selected, "\n")

  # ---- 3. presence 矩阵装配（species × ASV；块存该列 ASV 的 n_hosts）----
  pres <- lapply(presence_files, data.table::fread)
  sp_by_db <- list(); mat_n <- list()
  for (k in names(pres)) {
    sub <- all_sel[db == k]
    p  <- pres[[k]]
    pm <- as.matrix(p[ASV_ID %in% sub$ASV_ID, -"ASV_ID", with = FALSE])
    rownames(pm) <- sub$new_name[match(p[ASV_ID %in% sub$ASV_ID, ASV_ID], sub$ASV_ID)]
    pm <- t(pm)                       # species × ASV
    sp_by_db[[k]] <- rownames(pm)
    m <- pm * matrix(sub$n_hosts[match(colnames(pm), sub$new_name)],
                     nrow = nrow(pm), ncol = ncol(pm), byrow = TRUE)
    m[pm == 0] <- NA                  # 不存在/未测 -> NA（白）
    mat_n[[k]] <- m
    cat(sprintf("[gos-disp] %s presence species=%d  selected ASV cols=%d\n",
                k, nrow(pm), ncol(pm)))
  }

  tr <- ape::read.tree(tree_file)
  tree_sp <- tr$tip.label
  all_sp <- unique(unlist(sp_by_db))
  missing_in_tree <- setdiff(all_sp, tree_sp)
  if (length(missing_in_tree))
    cat("[gos-disp] WARN species in presence but NOT in tree:", length(missing_in_tree), "\n")
  sp_order <- c(tree_sp[tree_sp %in% all_sp], missing_in_tree)  # 树内按树序，树外补末尾

  comb <- matrix(NA, nrow = length(sp_order), ncol = nrow(all_sel),
                 dimnames = list(sp_order, all_sel$new_name))
  for (k in names(mat_n)) {
    idx <- intersect(sp_order, rownames(mat_n[[k]]))
    comb[idx, colnames(mat_n[[k]])] <- mat_n[[k]][idx, ]
  }
  cat("[gos-disp] combined matrix dim:", nrow(comb), "x", ncol(comb), "\n")

  # ---- 4. 元数据：species -> Order / Genus + Order 色带 ----
  meta <- lapply(meta_files, data.table::fread)
  sp_ord <- data.table::rbindlist(
    lapply(meta, function(m) unique(m[, .(species, Order, Genus)])))
  sp_ord <- unique(sp_ord[!is.na(Order) & Order != ""])
  dup <- sp_ord[, .N, by = species][N > 1, species]
  if (length(dup)) cat("[gos-disp] WARN species with >1 Order/Genus rows:", length(dup), "\n")
  sp_ord <- sp_ord[!duplicated(species)]

  cmf <- readRDS(color_map_rds)
  ord_col <- cmf$order
  sp_ord[, Order_col := ord_col[Order]]
  no_ord_col <- sp_ord[is.na(Order_col)]
  if (nrow(no_ord_col)) {
    cat("[gos-disp] WARN species lacking Order color:", nrow(no_ord_col), "\n")
    print(no_ord_col)
  }

  ord_vec <- sp_ord$Order[match(sp_order, sp_ord$species)]
  ord_col_vec <- sp_ord$Order_col[match(sp_order, sp_ord$species)]
  gen_vec <- sp_ord$Genus[match(sp_order, sp_ord$species)]
  names(ord_vec) <- names(ord_col_vec) <- names(gen_vec) <- sp_order

  # ---- 5. 热图配色 + 列聚类 + 注记 ----
  bar4 <- c("C_B" = "#1B6CA8", "D_B" = "#7FBFDF",
            "C_F" = "#238B45", "D_F" = "#8FD4A8")
  asv_grp <- paste(all_sel$grp, all_sel$db, sep = "_")
  names(asv_grp) <- all_sel$new_name

  # n_hosts 块色：log2 变换，低→高 浅黄→深棕红
  brk_l2 <- log2(sort(unique(all_sel$n_hosts)))
  col_n <- circlize::colorRamp2(brk_l2,
    grDevices::colorRampPalette(c("#FDEBD0", "#E67E22", "#641E16"))(length(brk_l2)))
  leg_at <- unique(round(brk_l2[pretty(brk_l2, n = 5)], 2))
  leg_at <- leg_at[leg_at >= min(brk_l2) & leg_at <= max(brk_l2)]

  comb_l2 <- log2(comb)                 # 块存 log2(n_hosts)；NA = 白
  pres_bin <- comb_l2; pres_bin[is.na(pres_bin)] <- 0; pres_bin[pres_bin > 0] <- 1
  grp_lev <- names(bar4)
  col_order_vec <- unlist(lapply(grp_lev, function(g) {
    idx <- which(asv_grp == g)
    if (length(idx) >= 2) {
      hc <- stats::hclust(stats::dist(t(pres_bin[, idx, drop = FALSE])), method = "complete")
      idx[hc$order]
    } else idx
  }))
  names(col_order_vec) <- NULL

  ra <- ComplexHeatmap::rowAnnotation(
    Order = ComplexHeatmap::anno_simple(ord_vec, col = ord_col, na_col = "grey90"),
    show_annotation_name = TRUE,
    annotation_name_rot = 90,
    annotation_name_gp = grid::gpar(fontsize = 8),
    col = list(Order = ord_col))
  ca <- ComplexHeatmap::columnAnnotation(
    Group = asv_grp,
    col = list(Group = bar4),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(3, "mm"))

  ht <- ComplexHeatmap::Heatmap(
    comb_l2,
    name = "n_hosts",
    col = col_n,
    na_col = "white",
    cluster_rows = FALSE,
    row_order = sp_order,
    cluster_columns = FALSE,
    column_order = col_order_vec,
    column_split = factor(asv_grp, levels = names(bar4)),
    column_gap = grid::unit(1.2, "mm"),
    top_annotation = ca,
    right_annotation = ra,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 6),
    show_column_names = FALSE,
    heatmap_legend_param = list(
      title = "n_hosts",
      at = leg_at,
      labels = sprintf("%.0f", 2^leg_at),
      direction = "vertical",
      grid_width = grid::unit(3, "mm")))

  for (fmt in c("pdf", "jpg")) {
    fn <- paste0(out_fig_base, ".", fmt)
    if (fmt == "pdf") grDevices::pdf(fn, width = fig_width, height = fig_height)
    else grDevices::jpeg(fn, width = jpg_width, height = jpg_height, res = 300)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE)
    grDevices::dev.off()
    cat("[gos-disp] figure written:", fn, "\n")
  }

  # ---- 6. clustered ASV 注释表 + Repres_genus / Repres_N ----
  tax <- lapply(tax_files, data.table::fread)
  tax_b <- data.table::rbindlist(list(B = tax$B, F = tax$F), idcol = "db", use.names = TRUE)
  data.table::setnames(tax_b, "taxid", "ASV_ID")
  data.table::setnames(tax_b, "Sub_Species", "Subspecies")

  clus <- all_sel[grp == "C"]          # 只取 clustered ASV（与热图口径一致）
  clus <- merge(clus, tax_b[, .(db, ASV_ID, Kingdom, Phylum, Class, Order, Family,
                                Genus, Species, Subspecies)],
                by = c("db", "ASV_ID"), all.x = TRUE)
  cat("[gos-disp] clustered ASVs:", nrow(clus),
      " taxonomy matched:", sum(!is.na(clus$Genus)), "\n")

  rep_fun <- function(asv_id, db_k) {
    p <- as.matrix(pres[[db_k]][ASV_ID == asv_id, -"ASV_ID", with = FALSE])[1, ]
    sp_pres <- names(p)[p == 1]
    if (!length(sp_pres)) return(list(g = NA_character_, n = 0L))
    g <- sp_ord$Genus[match(sp_pres, sp_ord$species)]
    tb <- table(g[!is.na(g) & g != "" & g != "Unassigned"])
    if (!length(tb)) return(list(g = NA_character_, n = 0L))
    mx <- max(tb)
    list(g = paste(sort(names(tb)[tb == mx]), collapse = ";"), n = as.integer(mx))
  }
  res <- Map(rep_fun, clus$ASV_ID, clus$db)
  clus[, Repres_genus := unlist(lapply(res, `[[`, "g"))]
  clus[, Repres_N := unlist(lapply(res, `[[`, "n"))]

  out_cols <- c("new_name", "dispersion", "n_hosts", "Kingdom", "Phylum", "Class",
                "Order", "Family", "Genus", "Species", "Subspecies", "Repres_genus", "Repres_N")
  data.table::fwrite(clus[, ..out_cols], out_taxonomy, sep = "\t")
  cat("[gos-disp] taxonomy table written:", nrow(clus), "rows ->", out_taxonomy, "\n")
  print(clus[, .N, by = .(db, n_hosts)])

  cat("\nGOS DISPERSION HEATMAP DONE. outputs:", out_selected, ";", out_taxonomy, ";",
      paste0(out_fig_base, ".{pdf,jpg}"), "\n")
  invisible(ht)
}
