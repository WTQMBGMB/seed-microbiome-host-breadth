# Shared diversity and host-phylogeny helpers / 共享多样性与宿主系统发育工具
# Adapted from the reviewed host_seed_utils.R; paths point to github2 data/results.
# ============================================================================
# host_seed_utils.R — Host_Seed_microbiota 分析工具函数
#
# 用途：为主分析脚本 host_seed_microbiota.R 提供可复用函数。
# 依赖：ape, vegan, ggtree, ggtreeExtra, phytools, cluster, ecodist,
#       ggplot2, patchwork, dplyr, tidyr, data.table, RColorBrewer
#
# 函数清单：
#   load_region_data       — 读取丰度表 + metadata + 宿主树，对齐后返回 list
#   calc_alpha_div         — 计算 Shannon/Simpson/Observed/Pielou
#   build_color_map        — 构建全局 Order→Family→Genus 统一颜色表
#   phylo_alpha_plot       — 环形宿主树 + 外圈 α 多样性柱状图 + 系统发育信号
#   layered_boxplot        — 分层箱线图（3 区域 × 3 分类层级）
#   mean_dist              — 样本距离矩阵 → 物种平均距离矩阵
#   mantel_mrm_analysis    — Mantel + MRM + dbRDA 完整分析
# ============================================================================


# ---- 辅助：格式化 p 值 ----
.format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# ---- 辅助：统一保存 PDF + JPG ----
save_fig <- function(p, stem, out_dir = FIG_DIR, w = 8, h = 6) {
  save_panel(p, stem, w, h)
}

# ============================================================================
# 1. load_region_data — 读取并对齐一个区域的所有数据
#
# 参数：
#   region    "v4" | "v5_v7" | "its1"
#   base_dir  Downstream_analysis 根目录
#   verbose   是否打印详细信息
#
# 返回：list(
#   otu        = sample×ASV 矩阵 (numeric matrix)
#   meta       = 对齐后的 metadata (data.frame)
#   tree       = 宿主系统发育树 (phylo)
#   region     = 区域名
#   n_samples  = 最终样本数
#   n_asv      = 最终 ASV 数
#   n_species  = 最终宿主物种数
#   dropped_tips = 被裁剪的 tip 名称 (character)
# )
# ============================================================================
load_region_data <- function(region, base_dir = ".", verbose = TRUE) {

  # ---- 1.1 路径映射 ----
  path_map <- list(
    v4    = list(abund = file.path(CSS_DIR, "region", "v4", "css_v4_asvtab.txt"),
                 meta  = file.path(DATA_DIR, "metadata", "metadata_16S.txt"),
                 tree  = file.path(TREE_DIR, "plant_species_tree_v4.nwk"),
                 region_filter = "V4", meta_type = "16S"),
    v5_v7 = list(abund = file.path(CSS_DIR, "region", "v5_v7", "css_v5_v7_asvtab.txt"),
                 meta  = file.path(DATA_DIR, "metadata", "metadata_16S.txt"),
                 tree  = file.path(TREE_DIR, "plant_species_tree_v5_v7.nwk"),
                 region_filter = "V5_V7", meta_type = "16S"),
    its1  = list(abund = file.path(CSS_DIR, "region", "its1", "css_its1_asvtab.txt"),
                 meta  = file.path(DATA_DIR, "metadata", "metadata_ITS.txt"),
                 tree  = file.path(TREE_DIR, "plant_species_tree_its1.nwk"),
                 region_filter = "ITS1", meta_type = "ITS")
  )
  cfg <- path_map[[region]]
  if (is.null(cfg)) stop("Unknown region: ", region)

  # ---- 1.2 读取丰度表 ----
  abund_path <- cfg$abund
  if (verbose) message(sprintf("[load_region_data] %s: reading %s", region, abund_path))
  raw <- data.table::fread(abund_path, sep = "\t", header = TRUE)
  asv_ids <- as.character(raw[[1]])
  otu_mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(otu_mat) <- asv_ids
  rm(raw); invisible(gc())

  # 转置为 sample × ASV
  otu <- t(otu_mat)
  rm(otu_mat); invisible(gc())

  if (verbose) message(sprintf("  abundance matrix: %d samples × %d ASVs", nrow(otu), ncol(otu)))

  # ---- 1.3 读取 metadata ----
  meta_path <- cfg$meta
  meta_full <- data.table::fread(meta_path, sep = "\t", header = TRUE)
  # 按 Region 筛选
  meta <- meta_full[meta_full$Region == cfg$region_filter, ]
  if (verbose) message(sprintf("  metadata: %d samples (filtered from %d by Region=%s)",
                               nrow(meta), nrow(meta_full), cfg$region_filter))
  rm(meta_full)

  # ---- 1.4 读取宿主树 ----
  tree_path <- cfg$tree
  if (verbose) message(sprintf("  reading tree: %s", tree_path))
  tree <- ape::read.tree(tree_path)
  if (verbose) message(sprintf("  tree: %d tips", length(tree$tip.label)))

  # ---- 1.5 对齐丰度表与 metadata ----
  # 统一样本名（丰度表列名可能含 X 前缀）
  meta_samples <- as.character(meta$SampleID)
  otu_samples  <- rownames(otu)

  common_samples <- intersect(otu_samples, meta_samples)
  if (length(common_samples) == 0) {
    # 尝试修复：检查是否是 X 前缀问题
    otu_fixed <- gsub("^X", "", otu_samples)
    common_samples <- intersect(otu_fixed, meta_samples)
    if (length(common_samples) > 0) {
      rownames(otu) <- otu_fixed
      otu_samples <- otu_fixed
      common_samples <- intersect(otu_samples, meta_samples)
    }
  }
  otu  <- otu[common_samples, , drop = FALSE]
  meta <- as.data.frame(meta[match(common_samples, meta$SampleID), ])
  rownames(meta) <- meta$SampleID

  # 清理零丰度 ASV
  asv_totals <- colSums(otu)
  zero_asv <- which(asv_totals <= 0)
  if (length(zero_asv) > 0) {
    if (verbose) message(sprintf("  removing %d zero-abundance ASVs", length(zero_asv)))
    otu <- otu[, -zero_asv, drop = FALSE]
  }

  if (verbose) message(sprintf("  after alignment: %d samples, %d ASVs", nrow(otu), ncol(otu)))

  # ---- 1.6 对齐宿主树 tip 与 metadata species ----
  tree_tips  <- tree$tip.label
  meta_species <- unique(meta$species[!is.na(meta$species) & meta$species != ""])

  # 匹配 tips 与 species
  common_species <- intersect(tree_tips, meta_species)
  dropped_tips   <- setdiff(tree_tips, meta_species)
  missing_species <- setdiff(meta_species, tree_tips)

  if (verbose) {
    message(sprintf("  tree tips: %d, meta species: %d, intersect: %d",
                    length(tree_tips), length(meta_species), length(common_species)))
    if (length(dropped_tips) > 0)
      message(sprintf("  dropped tips (%d): %s", length(dropped_tips),
                      paste(head(dropped_tips, 10), collapse = ", "),
                      if (length(dropped_tips) > 10) " ..." else ""))
    if (length(missing_species) > 0)
      message(sprintf("  species not in tree (%d): %s", length(missing_species),
                      paste(head(missing_species, 10), collapse = ", "),
                      if (length(missing_species) > 10) " ..." else ""))
  }

  # 裁剪树
  if (length(dropped_tips) > 0) {
    tree <- ape::drop.tip(tree, dropped_tips)
  }

  # 交集覆盖率检查
  coverage <- length(common_species) / length(tree_tips)
  if (coverage < 0.8) {
    warning(sprintf("[load_region_data] %s: tree tip coverage %.1f%% < 80%%", region, coverage * 100))
  }

  list(
    otu          = otu,
    meta         = meta,
    tree         = tree,
    region       = region,
    n_samples    = nrow(otu),
    n_asv        = ncol(otu),
    n_species    = length(common_species),
    dropped_tips = dropped_tips
  )
}

# ============================================================================
# 2. calc_alpha_div — 计算 α 多样性指数
#
# 参数：
#   otu  sample×ASV 矩阵（每行一个样本）
#
# 返回：data.frame(SampleID, Shannon, Simpson, Observed, Pielou)
# ============================================================================
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
#            palette_order = named_vector)
# ============================================================================
build_color_map <- function(all_meta) {

  # 提取所有唯一 Order
  orders <- sort(unique(all_meta$Order[!is.na(all_meta$Order) & all_meta$Order != ""]))
  n_orders <- length(orders)

  # 用 RColorBrewer Set3 (12色) + Paired (12色) + 自定义扩展
  brew_set3 <- RColorBrewer::brewer.pal(12, "Set3")
  brew_paired <- RColorBrewer::brewer.pal(12, "Paired")
  brew_set1 <- RColorBrewer::brewer.pal(9, "Set1")
  # 手动扩展色（来自 PAL2 经验色板）
  extra_colors <- c(
    "#855B67", "#B6BCCA", "#78A495", "#8A926B", "#5D6C57",
    "#995B70", "#846992", "#D4C0C2", "#F7E09B", "#099396",
    "#E7D6C2", "#C5D8E7", "#641E16", "#BFEFEE", "#D35400",
    "#1D3557", "#5B2C6F", "#E63946", "#355070", "#CB997E",
    "#2A9D8F", "#FFE8D6", "#6D597A", "#F9C74F", "#2E4053",
    "#BA4A00", "#7DCEA0", "#424949", "#E5989B", "#E07A5F"
  )
  all_palette <- unique(c(brew_set3, brew_paired, brew_set1, extra_colors))

  # 分配 Order 颜色
  order_colors <- setNames(all_palette[1:n_orders], orders)

  # 为每个 Order 下的 Family 生成渐变色
  family_colors <- c()
  genus_colors  <- c()

  for (ord in orders) {
    fams <- sort(unique(all_meta$Family[all_meta$Order == ord &
                                          !is.na(all_meta$Family) &
                                          all_meta$Family != ""]))
    n_fam <- length(fams)
    if (n_fam == 0) next

    base_color <- order_colors[ord]
    if (n_fam == 1) {
      fam_col <- setNames(base_color, fams)
    } else {
      # 从 Order 颜色到浅色渐变
      fam_pal <- colorRampPalette(c(base_color, "#F5F5F5"))(n_fam + 2)
      fam_col <- setNames(fam_pal[1:n_fam], fams)
    }
    family_colors <- c(family_colors, fam_col)

    # Genus 继承 Family 颜色
    for (fam in fams) {
      gen <- sort(unique(all_meta$Genus[all_meta$Family == fam &
                                          !is.na(all_meta$Genus) &
                                          all_meta$Genus != ""]))
      if (length(gen) > 0) {
        genus_colors <- c(genus_colors, setNames(rep(fam_col[fam], length(gen)), gen))
      }
    }
  }

  list(
    order  = order_colors,
    family = family_colors,
    genus  = genus_colors,
    palette_order = order_colors
  )
}

# ============================================================================
# 3b. build_color_map_flat — 新建：扁平配色（同 Order 的 Family/Genus 共用 Order 颜色）
#
# 与 build_color_map 的关系：原函数保留不动；本函数是另建的新方案，返回结构
# 完全一致（order/family/genus/palette_order），可直接互换使用。
#
# 策略：
#   - Order：使用用户推荐的高对比度初始色板（30 色，默认），减少相邻色混淆
#   - Family：一律使用所属 Order 的颜色（不再做白色渐变）
#   - Genus ：一律使用所属 Order 的颜色（继承 Family = Order 色）
#
# 参数：
#   all_meta      合并后的 metadata data.frame，包含 Order/Family/Genus 列
#   order_palette 自定义 Order 色板向量；NULL 时使用内置推荐 30 色
#
# 返回：list(order = named_vector, family = named_vector, genus = named_vector,
#            palette_order = named_vector)
# ============================================================================
build_color_map_flat <- function(all_meta, order_palette = NULL) {

  orders <- sort(unique(all_meta$Order[!is.na(all_meta$Order) & all_meta$Order != ""]))
  n_orders <- length(orders)

  # 用户推荐的高对比度初始 Order 色板（30 色）
  if (is.null(order_palette)) {
    order_palette <- c(
      "#3B528B", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#377EB8",
      "#FFFF33", "#A65628", "#F781BF", "#00A6A6", "#D73027", "#4575B4",
      "#91CF60", "#762A83", "#FEE08B", "#FDAE61", "#1B7837", "#C51B7D",
      "#5E4FA2", "#3288BD", "#66C2A5", "#FC8D59", "#E78AC3", "#A6D854",
      "#FFD92F", "#8DA0CB", "#E5C494", "#B3B3B3", "#1F78B4", "#B15928"
    )
  }
  # 色板不足时用 RColorBrewer 扩展（避免 NA 色）
  if (n_orders > length(order_palette)) {
    extra <- unique(c(RColorBrewer::brewer.pal(12, "Set3"),
                      RColorBrewer::brewer.pal(12, "Paired"),
                      RColorBrewer::brewer.pal(9, "Set1")))
    extra <- extra[!(extra %in% order_palette)]
    order_palette <- c(order_palette, extra)
  }
  order_colors <- setNames(order_palette[1:n_orders], orders)

  family_colors <- c()
  genus_colors  <- c()
  for (ord in orders) {
    base <- order_colors[ord]
    fams <- sort(unique(all_meta$Family[all_meta$Order == ord &
                                          !is.na(all_meta$Family) &
                                          all_meta$Family != ""]))
    if (length(fams) > 0) {
      family_colors <- c(family_colors, setNames(rep(base, length(fams)), fams))
    }
    for (fam in fams) {
      gen <- sort(unique(all_meta$Genus[all_meta$Family == fam &
                                          !is.na(all_meta$Genus) &
                                          all_meta$Genus != ""]))
      if (length(gen) > 0) {
        genus_colors <- c(genus_colors, setNames(rep(base, length(gen)), gen))
      }
    }
  }

  list(
    order  = order_colors,
    family = family_colors,
    genus  = genus_colors,
    palette_order = order_colors
  )
}

# ============================================================================
# 4. phylo_alpha_plot — 环形宿主树 + 外圈 α 多样性柱状图 + 系统发育信号
#
# 设计约定：外圈仅画纯柱状图（不画误差棒）。误差棒在 ggtree 下易错位，
# 难以调试，按用户要求一律取消；物种均值之间的差异由柱高直接体现。
#
# 参数：
#   data            load_region_data 的输出 list
#   alpha_idx       要展示的 α 多样性指数名称 ("Shannon" | "Simpson" | "Observed" | "Pielou")
#   color_map       build_color_map 的输出 list
#   out_dir         输出目录
#   prefix          输出文件名前缀（如 "phylo_alpha_v4"）
#   show_species    是否显示物种标签（默认 FALSE：物种多时标签拥挤，仅保留 Order 颜色点）
#   legend_pos      图例位置；传 "none" 则隐藏全部图例（供多面板拼接时单独绘制共享图例）
#   fill_limits     渐变填充的数值范围 c(low, high)【作图尺度，即 log10_transform 之后的值域】；
#                   传 NULL 时由数据自动决定
#   log10_transform 是否对柱状图 x 值（及填充色）做 log10 变换。适合跨度达数量级的
#                   计数型指数（如 Observed）：均值 ≤ 0 的物种无法取对数，会从柱环剔除
#                   （树上点仍保留），并打印提示。系统发育信号检验仍基于原始值。
#   save            是否在本函数内落盘（FALSE 时仅返回 plot，供外部拼图）
#
# 返回：invisible(list(plot, phylosig_K, phylosig_L, alpha_species))
# ============================================================================
phylo_alpha_plot <- function(data, alpha_idx = "Shannon", color_map = NULL,
                             out_dir = ".", prefix = "phylo_alpha",
                             show_species = FALSE, legend_pos = "right",
                             fill_limits = NULL, log10_transform = FALSE,
                             save = TRUE) {

  if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
    stop("需要安装 ggtreeExtra 包：BiocManager::install('ggtreeExtra')")
  }

  otu  <- data$otu
  meta <- data$meta
  tree <- data$tree

  # ---- 4.1 计算 α 多样性 ----
  alpha_df <- calc_alpha_div(otu)
  alpha_df <- merge(alpha_df, meta[, c("SampleID", "species", "Order", "Family", "Genus")],
                    by = "SampleID", all.x = TRUE)

  # ---- 4.2 过滤：仅保留树中存在的物种 ----
  alpha_df <- alpha_df[alpha_df$species %in% tree$tip.label, ]

  # ---- 4.3 物种水平聚合（均值 ± SE） ----
  alpha_sp <- alpha_df %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      mean_val = mean(.data[[alpha_idx]], na.rm = TRUE),
      se_val   = if (dplyr::n() > 1) sd(.data[[alpha_idx]], na.rm = TRUE) / sqrt(dplyr::n()) else NA_real_,
      n_sample = dplyr::n(),
      Order    = dplyr::first(Order),
      Family   = dplyr::first(Family),
      .groups  = "drop"
    )
  alpha_sp <- as.data.frame(alpha_sp)
  rownames(alpha_sp) <- alpha_sp$species

  # ---- 4.4 提取性状向量（按树 tip 顺序） ----
  tips_in_tree <- tree$tip.label
  trait <- setNames(alpha_sp[tips_in_tree, "mean_val"], tips_in_tree)
  trait <- trait[!is.na(trait)]

  # 裁剪树以匹配有数据的 tips
  valid_tips <- names(trait)
  if (length(valid_tips) < length(tree$tip.label)) {
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, valid_tips))
  }

  # ---- 4.5 系统发育信号检验 ----
  # log10 尺度与柱环保持一致 / use the same log10 scale as the bar ring.
  sig_trait <- trait
  tree_sig <- tree
  if (log10_transform) {
    sig_trait <- log10(sig_trait)
    keep_sig <- is.finite(sig_trait)
    if (!all(keep_sig)) {
      sig_trait <- sig_trait[keep_sig]
      tree_sig <- ape::drop.tip(tree, setdiff(tree$tip.label, names(sig_trait)))
    }
  }
  ps_K <- tryCatch(
    phytools::phylosig(tree_sig, sig_trait[tree_sig$tip.label], method = "K", test = TRUE),
    error = function(e) NULL
  )
  ps_L <- tryCatch(
    phytools::phylosig(tree_sig, sig_trait[tree_sig$tip.label], method = "lambda", test = TRUE),
    error = function(e) NULL
  )

  # 构建 subtitle
  sub_lines <- c()
  if (!is.null(ps_K)) {
    sub_lines <- c(sub_lines,
                   sprintf("Blomberg's K = %.3f, P = %s", ps_K$K, .format_p(ps_K$P)))
  }
  if (!is.null(ps_L)) {
    sub_lines <- c(sub_lines,
                   sprintf("Pagel's λ = %.3f, P = %s", ps_L$lambda, .format_p(ps_L$P)))
  }
  phylosig_subtitle <- paste(sub_lines, collapse = "; ")

  message(sprintf("[phylo_alpha_plot] %s: %s", prefix, phylosig_subtitle))

  # ---- 4.6 准备树注释数据 ----
  tip_order <- alpha_sp[tree$tip.label, "Order"]
  tip_order[is.na(tip_order)] <- "Unknown"
  tip_data <- data.frame(label = tree$tip.label, Order = tip_order,
                         stringsAsFactors = FALSE)

  # ---- 4.7 准备柱状图数据 ----
  # 仅保留画柱所需列，避免与树数据列名冲突（Order 等已通过 tip_data 传入）
  # 说明：se_val 仅保留在 alpha_sp 返回值中；柱状图不画误差棒，无需带入
  bar_data <- alpha_sp[tree$tip.label, c("species", "mean_val")]
  bar_data$species <- factor(bar_data$species, levels = tree$tip.label)

  # 可选 log10 变换（适合 Observed 等数量级跨度大的计数型指数）：
  # log10(0) = -Inf，故剔除均值 <= 0 的物种柱（树上的点仍保留），并打印提示。
  # 注意：fill_limits 须按变换后的尺度给出；系统发育信号检验仍用原始值。
  bar_axis_title <- alpha_idx
  if (log10_transform) {
    n_nonpos <- sum(bar_data$mean_val <= 0, na.rm = TRUE)
    if (n_nonpos > 0) {
      message(sprintf("  [%s] log10 transform: %d species with %s <= 0 dropped from bar ring",
                      prefix, n_nonpos, alpha_idx))
      bar_data <- bar_data[bar_data$mean_val > 0, , drop = FALSE]
    }
    bar_data$mean_val <- log10(bar_data$mean_val)
    bar_axis_title <- sprintf("log10(%s)", alpha_idx)
  }

  # ---- 4.8 绘制环形树 ----
  # 基础环形树 + Order 颜色映射
  p_tree <- ggtree::ggtree(tree, layout = "circular", size = 0.4)
  p_tree <- p_tree %<+% tip_data

  # 添加分支颜色（按 Order）
  p_tree <- p_tree + ggtree::geom_tippoint(ggplot2::aes(color = Order), size = 1.5)
  # species 标签可选：物种多时标签会拥挤无法辨认，默认关闭（仅保留 Order 颜色点）
  if (show_species) {
    p_tree <- p_tree +
      ggtree::geom_tiplab(ggplot2::aes(label = gsub("_", " ", label)),
                          size = 1.8, offset = 0.15, fontface = "italic")
  }

  if (!is.null(color_map) && !is.null(color_map$order)) {
    p_tree <- p_tree +
      ggplot2::scale_color_manual(values = color_map$order, na.value = "grey70")
  }

  # ---- 4.9 添加外圈 α 多样性柱状图 ----
  # 用 ggtreeExtra::geom_fruit 添加柱状图环
  alpha_range <- range(bar_data$mean_val, na.rm = TRUE)
  alpha_range[1] <- max(0, alpha_range[1] - 0.1 * diff(alpha_range))
  alpha_range[2] <- alpha_range[2] + 0.1 * diff(alpha_range)

  p <- p_tree +
    ggtreeExtra::geom_fruit(
      data = bar_data,
      geom = geom_col,
      mapping = ggplot2::aes(y = species, x = mean_val, fill = mean_val),
      orientation = "y",
      width = 0.7,
      offset = 0.08,
      pwidth = 0.25,
      axis.params = list(axis = "x", title = bar_axis_title,
                         text.size = 1.5, title.size = 2.5)
    ) +
    ggplot2::scale_fill_gradient(low = "#F7E09B", high = "#3A5A40",
                                 name = bar_axis_title, limits = fill_limits,
                                 guide = ggplot2::guide_colourbar(
                                   title.position = "top", barwidth = 8, barheight = 0.8
                                 ))

  # 注意：不添加误差棒。按用户要求（ggtree 调试困难），外圈仅保留纯柱状图。

  # ---- 4.10 添加标题 ----
  region_label <- data$region
  p <- p +
    ggplot2::labs(
      title = sprintf("Host Phylogeny & %s Diversity (%s)", alpha_idx, region_label),
      subtitle = phylosig_subtitle
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey40"),
      legend.position = legend_pos
    )

  # ---- 4.11 保存 ----
  if (save) save_fig(p, prefix, out_dir, w = 16, h = 14)

  # ---- 4.12 返回 ----
  invisible(list(
    plot = p,
    phylosig_K = ps_K,
    phylosig_L = ps_L,
    alpha_species = alpha_sp
  ))
}

# ============================================================================
# 5. layered_boxplot — 分层箱线图（3 区域 × 3 分类层级）
#
# 使用统一的 Order→Family→Genus 颜色映射，展示 α 多样性在各分类层级的分布。
# 每个区域独立成图（并排），区域内每个层级只保留有样本的分类单元：
#   - 每个区域 × 层级 内部按 α 中位数从高到低排序（不跨区域排序 / 不共享分类标签）
#   - 无样本的分类单元不绘制（不会出现空箱线图）
#   - 仅颜色使用统一的 color_map
#
# 参数：
#   all_data    list(region_name = list(alpha = alpha_df, meta = metadata_df))
#   color_map   build_color_map 的输出
#   alpha_idx   展示的 α 指数（默认 "Shannon"）
#   out_dir     输出目录
#   log10_transform  TRUE 时对 α 值轴做 log10 变换（适合跨度大的计数型指数如 Observed）
#
# 返回：invisible(list(plot = combined_plot, kw_results))
# ============================================================================
layered_boxplot <- function(all_data, color_map, alpha_idx = "Shannon",
                            out_dir = ".", prefix = "alpha_boxplot_layered",
                            log10_transform = FALSE) {

  region_names <- names(all_data)
  plot_list <- list()
  kw_list   <- list()
  max_total_taxa <- 0L   # 记录最大的区域总分类行数（用于自适应画布高度）

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
        ggplot2::aes(x = -Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = -0.1, vjust = 1.5,
        size = 2.8, color = "grey30"
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

  # ---- 5.5 共享横轴（α 值）范围，跨区域对齐 ----
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
  save_fig(p_combined, prefix, out_dir, w = 18, h = h_out)

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
# ============================================================================
mean_dist <- function(indata, classify_df) {
  if (inherits(indata, "dist")) {
    indata <- as.matrix(indata)
  }
  colnames(classify_df) <- c("SampleID", "groups")
  common_list <- intersect(rownames(indata), classify_df$SampleID)
  classify_df <- classify_df[classify_df$SampleID %in% common_list, , drop = FALSE]
  indata <- indata[common_list, common_list, drop = FALSE]

  group_list <- unique(classify_df$groups)
  group_list_num <- length(group_list)
  group_matrix <- matrix(data = 1, ncol = group_list_num, nrow = group_list_num)
  colnames(group_matrix) <- group_list
  rownames(group_matrix) <- group_list

  for (i in seq_len(group_list_num)) {
    row_group <- group_list[i]
    row_id <- classify_df$SampleID[classify_df$groups == row_group]
    for (j in seq_len(group_list_num)) {
      col_group <- group_list[j]
      col_id <- classify_df$SampleID[classify_df$groups == col_group]
      group_matrix[i, j] <- mean(indata[row_id, col_id, drop = FALSE])
    }
  }
  return(group_matrix)
}

# ============================================================================
# 7. mantel_mrm_analysis — Mantel + MRM + dbRDA 完整分析
#
# 参数：
#   data       load_region_data 的输出 list
#   out_dir    输出目录
#   bray_dist  可选：预计算的 Bray-Curtis 距离矩阵（避免重复计算）
#
# 返回：invisible(list(mantel, mrm, dbrda, plots, bray_dist))
# ============================================================================
mantel_mrm_analysis <- function(data, out_dir = ".", bray_dist = NULL) {

  otu  <- data$otu
  meta <- data$meta
  tree <- data$tree
  region_name <- data$region

  message(sprintf("[mantel_mrm_analysis] Processing %s ...", region_name))

  # ---- 7.1 计算 Bray-Curtis 距离 ----
  if (is.null(bray_dist)) {
    message("  computing Bray-Curtis distance (sample level)...")
    bray_dist <- vegan::vegdist(otu, method = "bray")
  } else {
    message("  using precomputed Bray-Curtis distance...")
  }
  bray_matrix <- as.matrix(bray_dist)

  # ---- 7.2 宿主遗传距离 ----
  tree_matrix <- stats::cophenetic(tree)
  message(sprintf("  cophenetic matrix: %d × %d species", nrow(tree_matrix), ncol(tree_matrix)))

  # ---- 7.3 物种水平平均距离 ----
  message("  aggregating to species-level mean distances...")
  classify_df <- meta[, c("SampleID", "species")]
  species_matrix <- mean_dist(bray_matrix, classify_df)

  # 对齐两个矩阵
  common_sp <- intersect(rownames(tree_matrix), rownames(species_matrix))
  message(sprintf("  common species: %d", length(common_sp)))
  if (length(common_sp) < 10) {
    warning(sprintf("[mantel_mrm_analysis] %s: only %d common species, skipping", region_name, length(common_sp)))
    return(invisible(NULL))
  }
  tree_mat_aligned    <- tree_matrix[common_sp, common_sp]
  species_mat_aligned <- species_matrix[common_sp, common_sp]

  # ---- 7.4 Mantel 检验 ----
  message("  running Mantel test...")
  mantel_res <- vegan::mantel(tree_mat_aligned, species_mat_aligned,
                              method = "pearson", permutations = 999)
  mantel_r <- mantel_res$statistic
  mantel_p <- mantel_res$signif

  # ---- 7.5 MRM ----
  message("  running MRM...")
  mrm_res <- ecodist::MRM(as.dist(species_mat_aligned) ~ as.dist(tree_mat_aligned),
                          nperm = 999)
  mrm_r2  <- mrm_res$r.squared[1]
  mrm_r2p <- mrm_res$r.squared[2]
  mrm_F   <- mrm_res$F.test[1]
  mrm_Fp  <- mrm_res$F.test[2]

  # ---- 7.6 dbRDA ----
  message("  running dbRDA...")
  pcoa_tree <- cmdscale(tree_mat_aligned, k = min(5, nrow(tree_mat_aligned) - 1), eig = TRUE)
  tree_vectors <- as.data.frame(pcoa_tree$points)
  dbrda_res <- vegan::capscale(species_mat_aligned ~ ., data = tree_vectors)
  dbrda_anova_all <- anova(dbrda_res, permutations = 999)
  dbrda_r2 <- sum(dbrda_res$CCA$eig) / dbrda_res$tot.chi
  dbrda_p  <- dbrda_anova_all$`Pr(>F)`[1]

  # ---- 7.7 可视化 ----
  # 散点图数据
  get_df <- function(tree_mat, species_mat) {
    idx <- upper.tri(tree_mat)
    data.frame(
      phylo = tree_mat[idx],
      micro = species_mat[idx]
    )
  }
  df_plot <- get_df(tree_mat_aligned, species_mat_aligned)

  # 颜色：16S 用 #4E79A7，ITS 用 #81B29A
  is_its <- grepl("its", region_name, ignore.case = TRUE)
  pt_color <- if (is_its) "#81B29A" else "#4E79A7"

  # Mantel+MRM 标注
  mantel_label <- sprintf(
    "Mantel r = %.3f\nP = %s\n\nMRM R² = %.3f\nP = %s",
    mantel_r, .format_p(mantel_p), mrm_r2, .format_p(mrm_r2p)
  )

  # Mantel+MRM 散点图
  p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = phylo, y = micro)) +
    ggplot2::geom_point(color = pt_color, alpha = 0.2, size = 1.2) +
    ggplot2::geom_smooth(method = "lm", color = pt_color, se = FALSE, linewidth = 1.2) +
    ggplot2::annotate("text", x = Inf, y = -Inf, label = mantel_label,
                      hjust = 1.05, vjust = -0.2, size = 3.5) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 13)) +
    ggplot2::labs(title = sprintf("Mantel & MRM (%s)", region_name),
                  x = "Host phylogenetic distance",
                  y = "Microbial community dissimilarity (Bray-Curtis)")

  # dbRDA 排序图
  # 注意：必须用 vegan::scores。主脚本 library() 顺序里 phytools 在 vegan 之后 attach，
  # 其导出的 scores 泛型会遮蔽 vegan 的 scores，导致 scores(capscale) 走 scores.default
  # 返回 character 向量（表现为 "object 'CAP1' not found"）。
  dbrda_label <- sprintf("dbRDA R² = %.3f\nP = %s", dbrda_r2, .format_p(dbrda_p))
  sc <- vegan::scores(dbrda_res, display = "sites")
  if (is.list(sc)) sc <- sc[[grep("site", names(sc), value = TRUE)[1]]]
  scores_mat <- as.matrix(sc)
  n_cap <- length(dbrda_res$CCA$eig)
  if (n_cap >= 2 && ncol(scores_mat) >= 2) {
    colnames(scores_mat)[1:2] <- c("CAP1", "CAP2")   # 兜底重命名，防轴命名差异
    scores_df <- as.data.frame(scores_mat[, 1:2, drop = FALSE])
    p2 <- ggplot2::ggplot(scores_df, ggplot2::aes(x = CAP1, y = CAP2)) +
      ggplot2::geom_point(color = pt_color, alpha = 0.7, size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
      ggplot2::annotate("text", x = Inf, y = Inf, label = dbrda_label,
                        hjust = 1.1, vjust = 1.5, size = 3.5) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold", size = 13)) +
      ggplot2::labs(title = sprintf("dbRDA (%s)", region_name),
                    x = "CAP1", y = "CAP2")
  } else {
    # 约束轴不足 2 个：无法做二维排序图，用占位图并提示
    warning(sprintf(
      "[mantel_mrm %s] dbRDA has only %d constrained axis (R²=%.3f) — ordination plot skipped",
      region_name, n_cap, dbrda_r2))
    p2 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = sprintf("dbRDA ordination not drawn\n(%d constrained axis)\nR² = %.3f",
                                        n_cap, dbrda_r2),
                        size = 4, color = "grey50") +
      ggplot2::labs(title = sprintf("dbRDA (%s)", region_name)) +
      ggplot2::theme_void()
  }

  # ---- 7.8 返回 ----
  results <- list(
    region       = region_name,
    mantel_r     = mantel_r,
    mantel_p     = mantel_p,
    mrm_r2       = mrm_r2,
    mrm_r2p      = mrm_r2p,
    mrm_F        = mrm_F,
    mrm_Fp       = mrm_Fp,
    dbrda_r2     = dbrda_r2,
    dbrda_p      = dbrda_p,
    n_species    = length(common_sp),
    plot_mantel  = p1,
    plot_dbrda   = p2,
    bray_dist    = bray_dist
  )

  invisible(results)
}
