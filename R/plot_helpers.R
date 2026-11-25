# =============================================================================
# Plotting helper / 绘图辅助函数
#
# create_color_legend() reproduces the established legend layout because
# Chapter 1 and Chapter 3 use its exact layout. / 该函数按旧流程原样保留，供 Ch1/Ch3 图例使用。
# =============================================================================

create_color_legend <- function(data, save_path = "legend.pdf", shape = "circle", 
                                layout = "col", n = 2, title = "title", 
                                label_size = 4, title_size = 6, ifsave = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包: install.packages('ggplot2')")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("请先安装grid包: install.packages('grid')")
  }
  
  library(ggplot2)
  library(grid)
  
  # 数据预处理与检查
  unique_data <- unique(data)
  if (any(duplicated(unique_data[,1]))) {
    warning("去重后标签列仍有重复元素，请检查数据")
  }
  if (any(duplicated(unique_data[,2]))) {
    warning("去重后颜色列仍有重复元素，请检查数据")
  }
  colnames(unique_data) <- c("category", "color")
  unique_data$color <- as.character(unique_data$color)
  color_check <- grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", unique_data$color)
  if (!all(color_check)) {
    stop("颜色列包含非16进制格式的颜色值")
  }
  title_size = title_size*4
  n_items <- nrow(unique_data)
  
  # 关键增强：计算最长标签字符数，并动态调整色块间距
  max_label_length <- max(nchar(as.character(unique_data$category)))
  # 根据规则计算间距调整量：每10个字符增加1单位间距
  width_adjust_factor <- floor(max_label_length / 5) 
  # 基础间距 + 动态调整量
  base_width_per_col <- 1.5
  adjusted_width_per_col <- base_width_per_col + width_adjust_factor
  
  # 根据布局模式计算参数
  if (layout == "col") {
    ncol <- n
    nrow <- ceiling(n_items / ncol)
    # 使用动态调整后的间距计算图形宽度
    width <- ncol * (2+width_adjust_factor*0.5)
    height <- nrow * 1 + 1.0
    legend_direction <- "vertical"
  } else if (layout == "row") {
    nrow <- n
    ncol <- ceiling(n_items / nrow)
    width <- ncol * (1.5+width_adjust_factor*0.5)
    height <- nrow * 1.5 + 1.0
    legend_direction <- "horizontal"
  } else {
    stop("layout参数必须是'col'或'row'")
  }
  
  # 设置图形符号
  shapes <- list(
    circle = 21,
    square = 22,
    diamond = 23,
    triangle = 24
  )
  
  if (!shape %in% names(shapes)) {
    stop("不支持的形状，请选择: ", paste(names(shapes), collapse = ", "))
  }
  shape_code <- shapes[[shape]]
  
  # 创建图例数据框
  legend_df <- data.frame(
    category = unique_data$category,
    color = unique_data$color,
    x = rep(1:ncol, length.out = n_items),
    y = rep(nrow:1, each = ncol)[1:n_items]
  )
  
  # 关键增强：在col模式下，根据调整后的间距重新计算x轴位置
  if (layout == "col") {
    legend_df <- data.frame(
      category = unique_data$category,
      color = unique_data$color,
      x = rep(1:ncol, each = nrow)[1:n_items],
      y = rep(nrow:1, length.out = n_items)
    )
    # 计算每个色块的新x坐标，考虑动态间距
    x_positions <- numeric(ncol)
    for (i in 1:ncol) {
      x_positions[i] <- (i - 1) * adjusted_width_per_col + 1
    }
    x_value = c()
    timess = 0
    for (cols in c(1:ncol)){
      timess = timess+nrow
      if (timess<n_items){
        x_value = c(x_value,rep(x_positions[[cols]],nrow))
      }else{
        rest = n_items + nrow - timess
        x_value = c(x_value,rep(x_positions[[cols]],rest))
      }
    }
    legend_df$x <- x_value
  }
  
  # 创建图形
  p <- ggplot(legend_df, aes(x = x, y = y, fill = category)) +
    geom_point(shape = shape_code, size = 8, color = "black", aes(fill = category)) +
    scale_fill_manual(values = setNames(unique_data$color, unique_data$category)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = title_size, hjust = 0.5, face = "bold")
    ) +
    ggtitle(title) +
    coord_fixed(ratio = 1)
  
  # 添加标签
  if (layout == "col") {
    p <- p + geom_text(aes(label = category, x = x+0.3*ncol, y = y), 
                       size = label_size, hjust = 0, vjust = 0.5)+
      theme(plot.margin = margin(t = 10,  # 顶部边距
                                 r = 10+adjusted_width_per_col*5,  # 右侧边距
                                 b = 10,  # 底部边距
                                 l = 10,  # 左侧边距
                                 unit = "pt")  # 单位，如"pt"（磅）、"cm"（厘米）
      )+coord_cartesian(xlim =c(0.7, max(legend_df$x)+adjusted_width_per_col*0.8))
  } else if (layout == "row") {
    p <- p + geom_text(aes(label = category, x = x, y = y - 0.3), 
                       size = label_size, hjust = 0.5, vjust = 1)+
      coord_cartesian(xlim =c(0.5, max(legend_df$x)+adjusted_width_per_col*0.1),ylim =c(0.5, max(legend_df$y)))+
      theme(plot.margin = margin(t = 10,  # 顶部边距
                                 r = 10+adjusted_width_per_col*7,  # 右侧边距
                                 b = 10+adjusted_width_per_col*2,  # 底部边距
                                 l = 10+adjusted_width_per_col*7,  # 左侧边距
                                 unit = "pt"))
  }
  p
  # 保存图形
  width = width/1.4
  height = height/2
  if (ifsave){
    ggsave(save_path, p, width = width, height = height, units = "in", dpi = 300)
    message("图例已保存至: ", normalizePath(save_path))
    message("图例尺寸: ", round(width, 2), " x ", round(height, 2), " 英寸")
  }
  message("包含 ", n_items, " 个图例项")
  message("最长标签字符数: ", max_label_length, "，间距调整量: +", width_adjust_factor)
  
  return(p)
}

