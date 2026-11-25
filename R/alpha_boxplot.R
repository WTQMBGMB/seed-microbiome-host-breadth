# =============================================================================
# Alpha-diversity boxplots / α 多样性箱线图
#
# calc_alpha_div() and layered_boxplot() reproduce the established layout
# because Chapter 2 reproduces its exact layout and Kruskal-Wallis annotations.
# 保留原流程的实现，以保证 Ch2 箱线图与 K-W 标注完全一致。
# =============================================================================

.save_alpha_fig <- function(p, prefix, out_dir, width_in, height_in, dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, ".pdf")), p,
                  width = width_in, height = height_in, units = "in",
                  device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, ".jpg")), p,
                  width = width_in, height = height_in, units = "in",
                  device = ragg::agg_jpeg, dpi = dpi, quality = 95, bg = "white")
  invisible(file.path(out_dir, prefix))
}

.format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

calc_alpha_div <- function(otu) {
  # 确保输入是 matrix
  if (!is.matrix(otu)) otu <- as.matrix(otu)

  shannon  <- vegan::diversity(otu, index = "shannon", MARGIN = 1)
  simpson  <- vegan::diversity(otu, index = "simpson", MARGIN = 1)
  observed <- rowSums(otu > 0)
  pielou   <- ifelse(observed > 1, shannon / log(observed), NA_real_)

  df <- data.frame(
    SampleID = rownames(otu),
    Shannon  = shannon,
    Simpson  = simpson,
    Observed = observed,
    Pielou   = pielou,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  return(df)
}

# ============================================================================
# 3. build_color_map — 构建全局 Order→Family→Genus 统一颜色表
#
# 策略：
#   - Order：独立颜色（RColorBrewer Set3 扩展，确保每个 Order 颜色可区分）
#   - Family：所属 Order 颜色的渐变（colorRampPalette 在 Order 色与白色之间）
#   - Genus：继承所属 Family 颜色
#
# 参数：
#   all_meta  合并后的 metadata data.frame，包含 Order/Family/Genus/species 列
#
# 返回：list(order = named_vector, family = named_vector, genus = named_vector,

layered_boxplot <- function(all_data, color_map, alpha_idx = "Shannon",
                            out_dir = ".", prefix = "alpha_boxplot_layered",
                            log10_transform = FALSE) {

  region_names <- names(all_data)
  plot_list <- list()
  kw_list   <- list()
  max_total_taxa <- 0L   # 记录最大的区域总分类行数（用于自适应画布高度）

  # ---- 5.0 共享横轴（α 值）范围，跨区域对齐 ----
  # K-W 检验与中位数排序均在原始值上进行：log10 单调，秩不变，结论一致
  alpha_vals <- unlist(lapply(all_data, function(d) d$alpha[[alpha_idx]]),
                       use.names = FALSE)
  if (log10_transform) {
    # log10 轴：仅用正值取 log 后的范围再回映为原始单位，避免 0 → -Inf
    pos <- alpha_vals[is.finite(alpha_vals) & alpha_vals > 0]
    log_r  <- range(log10(pos))
    lg_pad <- 0.05 * diff(log_r)
    x_range <- 10^c(log_r[1] - lg_pad, log_r[2] + lg_pad)
  } else {
    x_range <- range(alpha_vals, na.rm = TRUE)
    x_pad <- 0.05 * diff(x_range)
    x_range <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  }

  # ---- 5.1 每个区域独立构建长表 + 层级内排序 ----
  for (reg in region_names) {
    dat  <- all_data[[reg]]
    df   <- dat$alpha
    meta <- dat$meta

    # 补齐分类列：alpha_df 已有则直接用，否则从 meta 合并
    tax_cols <- c("Order", "Family", "Genus")
    miss_tax <- setdiff(tax_cols, colnames(df))
    if (length(miss_tax) > 0) {
      df <- merge(df, meta[, c("SampleID", miss_tax), drop = FALSE],
                  by = "SampleID", all.x = TRUE)
    }

    # 三层层级长表（只保留有分类、α 非缺失的样本）
    long <- rbind(
      data.frame(SampleID = df$SampleID, alpha = df[[alpha_idx]],
                 Level = "Genus",  Taxon = df$Genus),
      data.frame(SampleID = df$SampleID, alpha = df[[alpha_idx]],
                 Level = "Family", Taxon = df$Family),
      data.frame(SampleID = df$SampleID, alpha = df[[alpha_idx]],
                 Level = "Order",  Taxon = df$Order)
    )
    long <- long[!is.na(long$Taxon) & long$Taxon != "" &
                 !is.na(long$alpha) & is.finite(long$alpha), , drop = FALSE]
    if (log10_transform) {
      # log10 轴无法表示 0：剔除并记录（K-W/中位数排序仍基于剩余样本）
      n_zero <- sum(long$alpha <= 0)
      if (n_zero > 0) {
        message(sprintf("  [%s] log10 axis: dropping %d rows with %s <= 0",
                        reg, n_zero, alpha_idx))
        long <- long[long$alpha > 0, , drop = FALSE]
      }
    }
    if (nrow(long) == 0) next

    # 每个 Level 内按中位数升序排列（coord_flip 后最低在底部、最高在顶部）
    med <- tapply(long$alpha, list(long$Level, long$Taxon),
                  stats::median, na.rm = TRUE)
    lvl_taxa <- lapply(levels(factor(long$Level)), function(lv) {
      v <- med[lv, , drop = FALSE]
      tx <- colnames(v)[!is.na(v[1, ])]
      if (length(tx) == 0) return(character())
      tx[order(v[1, tx])]                       # 升序：最低→底部，最高→顶部
    })
    long$Taxon <- factor(long$Taxon, levels = unlist(lvl_taxa, use.names = FALSE))
    # Level 加数字前缀：保证字母序 == 期望层级序（Genus→Family→Order）。
    # 若不加前缀，geom_text 层的字符型 Level 会在 facet 时被 factor 成字母序
    # (Family,Genus,Order)，与主数据层因子序冲突，把面板顺序打乱。
    # 前缀字符用户在 PDF 中可自行查找替换清理。
    long$Level <- factor(long$Level,
                         levels = c("Genus", "Family", "Order"),
                         labels = c("1_Genus", "2_Family", "3_Order"))

    # 不限制每个层级的分类数量：显示该区域内全部有样本的分类单元
    keep <- unlist(lvl_taxa, use.names = FALSE)
    long <- long[as.character(long$Taxon) %in% keep, , drop = FALSE]
    max_total_taxa <- max(max_total_taxa, length(unique(long$Taxon)))

    # ---- 5.2 每层级 Kruskal-Wallis 检验 ----
    kw_rows <- lapply(levels(long$Level), function(lv) {
      sub <- long[long$Level == lv, ]
      p <- tryCatch(
        suppressWarnings(stats::kruskal.test(alpha ~ Taxon, data = sub))$p.value,
        error = function(e) NA_real_
      )
      data.frame(Region = reg, Level = lv, kw_p = p,
                 n_taxa = length(unique(sub$Taxon)))
    })
    kw_df <- do.call(rbind, kw_rows)
    kw_df$label <- sprintf("K-W P = %s", sapply(kw_df$kw_p, .format_p))
    # 标注坐标（注意 coord_flip 后：aes x = 分类轴竖直方向，y = α 水平方向）
    # 不能再用 x = -Inf, y = Inf：两者都落在面板外，会被 panel clip 裁掉
    # （PDF 文本流里仍能提取到，但图上看不见）；设置 scale limits 后
    # ±Inf 还会被当作越界值直接删除。
    # 改为面板右下角内部：竖直方向放到该层级最底部分类之下 0.45 个单位
    # （离散轴默认扩展 0.6，位置 0.55 仍在面板内），水平方向放到 α 轴上界
    # （scale 默认再扩展 5%，故上界也在面板内）。
    kw_df$tax_pos   <- 0.55
    kw_df$alpha_pos <- x_range[2]
    # 关键：kw_df$Level 须与主数据层 Level 因子水平一致（含数字前缀），
    # 否则 geom_text 层按字符字母序 factor 化会把 facet 面板顺序打乱。
    kw_df$Level <- factor(kw_df$Level, levels = levels(long$Level))
    kw_list[[reg]] <- kw_df

    # ---- 5.3 颜色：统一 color_map（无匹配则灰色） ----
    tx_colors <- vapply(unique(as.character(long$Taxon)), function(tx) {
      col <- color_map$genus[tx]
      if (is.na(col)) col <- color_map$family[tx]
      if (is.na(col)) col <- color_map$order[tx]
      if (is.na(col)) col <- "#CCCCCC"
      unname(col)
    }, character(1))

    # ---- 5.4 绘图（区域独立：自己的标签与排序，仅颜色统一） ----
    p <- ggplot2::ggplot(long, ggplot2::aes(x = Taxon, y = alpha, fill = Taxon)) +
      ggplot2::geom_boxplot(outlier.size = 0.5, linewidth = 0.3, alpha = 0.85) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.15, size = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(Level ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_manual(values = tx_colors, guide = "none") +
      ggplot2::geom_text(
        data = kw_df,
        ggplot2::aes(x = tax_pos, y = alpha_pos, label = label),
        inherit.aes = FALSE, hjust = 1.05, vjust = 0,
        size = 2.8, color = "black"
      ) +
      ggplot2::labs(title = reg, x = "",
                    y = if (log10_transform) sprintf("log10(%s)", alpha_idx)
                        else alpha_idx) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.y = ggplot2::element_text(size = 7),
        legend.position = "none"
      )
    plot_list[[reg]] <- p
  }

  # ---- 5.5 应用共享横轴（α 值）范围，跨区域对齐（范围在 5.0 已计算） ----
  for (reg in names(plot_list)) {
    if (log10_transform) {
      plot_list[[reg]] <- plot_list[[reg]] +
        ggplot2::scale_y_continuous(trans = "log10", limits = x_range)
    } else {
      plot_list[[reg]] <- plot_list[[reg]] +
        ggplot2::scale_y_continuous(limits = x_range)
    }
  }

  # ---- 5.6 组合并保存 ----
  if (length(plot_list) == 1) {
    p_combined <- plot_list[[1]]
  } else {
    p_combined <- patchwork::wrap_plots(plot_list, ncol = length(plot_list)) +
      patchwork::plot_annotation(
        title = sprintf("Seed Microbial %s Diversity across Host Taxonomic Levels",
                        if (log10_transform) sprintf("log10(%s)", alpha_idx) else alpha_idx),
        subtitle = sprintf(
          "Ranks by median %s within each region x level (high to low); all taxa shown; unified colors; K-W P annotated; %s axis",
          alpha_idx, if (log10_transform) "log10" else "linear")
      )
  }
  # 画布高度自适应：无分类数量限制后，行数多的区域需要更高画布
  h_out <- max(20, ceiling(5 + max_total_taxa * 0.18))
  .save_alpha_fig(p_combined, prefix, out_dir, 18, h_out)

  invisible(list(plot = p_combined,
                 kw_results = do.call(rbind, kw_list)))
}

# ============================================================================
# 6. mean_dist — 将 sample×sample 距离矩阵聚合为 species×species 平均距离
#
# 参考：backup_20260810/obsolete_scripts/reference_code.R 第 299-321 行
#
# 参数：
#   indata       sample×sample 距离矩阵 (matrix 或 dist)
#   classify_df  两列 data.frame: SampleID, groups(物种名)
#
# 返回：group×group 平均距离 matrix（对称，行/列名为 group）
