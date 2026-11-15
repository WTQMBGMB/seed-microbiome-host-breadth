# CSS normalization, taxonomy aggregation and host trees / CSS、分类聚合与宿主树

css_normalize <- function(abund_path, out_path, norm_factor_path = NULL) {
  tab <- data.table::fread(abund_path, sep = "\t", header = TRUE)
  ids <- as.character(tab[[1]])
  otu <- data.matrix(tab[, -1, with = FALSE])
  rownames(otu) <- ids
  rm(tab)
  bad <- which(!vapply(seq_len(ncol(otu)), function(j) {
    x <- otu[, j]
    all(is.finite(x)) && all(x >= 0) && all(x == floor(x))
  }, logical(1)))
  if (length(bad)) stop("CSS requires integer counts / CSS 需要整数计数: ",
                        paste(colnames(otu)[head(bad, 10)], collapse = ", "))
  mre <- metagenomeSeq::newMRexperiment(otu)
  p <- metagenomeSeq::cumNormStatFast(mre)
  mre <- metagenomeSeq::cumNorm(mre, p = p)
  css <- metagenomeSeq::MRcounts(mre, norm = TRUE, log = FALSE)
  out <- data.table::as.data.table(css, keep.rownames = TRUE)
  data.table::setnames(out, "rn", "ASV_ID")
  data.table::fwrite(out, out_path, sep = "\t", quote = FALSE)
  if (!is.null(norm_factor_path)) {
    nf <- metagenomeSeq::normFactors(mre)
    data.table::fwrite(data.frame(sample = names(nf), css_quantile = as.numeric(nf)),
                       norm_factor_path, sep = "\t", quote = FALSE)
  }
  invisible(out)
}

aggregate_taxa_levels <- function(css_path, tax_path, out_dir, dataset,
                                  levels = c("Phylum", "Class", "Order", "Family", "Genus")) {
  css <- data.table::fread(css_path, sep = "\t", header = TRUE)
  tax <- data.table::fread(tax_path, sep = "\t", header = TRUE)
  ids <- as.character(css[[1]])
  if (!setequal(ids, as.character(tax[[1]]))) {
    stop("CSS and taxonomy IDs differ / CSS 与 taxonomy 的 ASV ID 不一致")
  }
  tax <- tax[match(ids, as.character(tax[[1]])), ]
  mat <- as.matrix(css[, -1, with = FALSE])
  rownames(mat) <- ids
  files <- character(length(levels))
  for (i in seq_along(levels)) {
    lv <- levels[[i]]
    grp <- as.character(tax[[lv]])
    grp[is.na(grp) | !nzchar(grp)] <- "Unassigned"
    sums <- rowsum(mat, group = grp, reorder = TRUE)
    out <- data.table::as.data.table(sums, keep.rownames = TRUE)
    data.table::setnames(out, "rn", lv)
    f <- file.path(out_dir, paste0(tolower(lv), "_", dataset, "_asvtab.txt"))
    data.table::fwrite(out, f, sep = "\t", quote = FALSE)
    files[[i]] <- f
  }
  invisible(files)
}

standardize_sample_names <- function(css_dir) {
  files <- list.files(css_dir, pattern = "_asvtab\\.txt$", recursive = TRUE, full.names = TRUE)
  for (f in files) {
    x <- data.table::fread(f, sep = "\t", header = TRUE)
    old <- names(x)[-1]
    new <- gsub("-", ".", old, fixed = TRUE)
    new <- ifelse(grepl("^[0-9]", new), paste0("X", new), new)
    data.table::setnames(x, old, new)
    data.table::fwrite(x, f, sep = "\t", quote = FALSE)
  }
  invisible(files)
}

build_host_trees <- function(tree_in, metadata_16s, metadata_its, out_dir) {
  m16 <- data.table::fread(metadata_16s, sep = "\t", header = TRUE)
  mits <- data.table::fread(metadata_its, sep = "\t", header = TRUE)
  tax <- unique(rbind(m16[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "species", "timetree")],
                      mits[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "species", "timetree")]))
  tree <- ape::read.tree(tree_in)
  tree$tip.label <- tax$species[match(tree$tip.label, tax$timetree)]
  ape::write.tree(tree, file.path(out_dir, "plant_species_tree_all.nwk"))
  cut <- function(meta, region, name) {
    tips <- unique(meta$species[meta$Region == region])
    t <- ape::keep.tip(tree, tips)
    ape::write.tree(t, file.path(out_dir, name))
  }
  cut(m16, "V4", "plant_species_tree_v4.nwk")
  cut(m16, "V5_V7", "plant_species_tree_v5_v7.nwk")
  cut(mits, "ITS1", "plant_species_tree_its1.nwk")
  one <- tax[!duplicated(tax$Genus), ]
  g <- ape::keep.tip(tree, one$species)
  g$tip.label <- one$Genus[match(g$tip.label, one$species)]
  ape::write.tree(g, file.path(out_dir, "plant_Genus_tree_all.nwk"))
  invisible(TRUE)
}
draw_pie_chart <- function(dataframe, 
                           data_col,
                           label = TRUE,
                           pre = TRUE, 
                           legend_t = "legend title", 
                           limit_n = 26,
                           color_list = c("#F7E09B","#93CEDF","#9EB395", "#E9BDA0", "#D5ABBF", "#846992", "#D98267", "#73BEA0", "#E4A644", "#A23F25",
                                          "#3E4069","#517D7E","#CA820A","#C1D2E2","#808DAA","#E5A984","#552379","#B24057","#A49BC6","#D35A21",
                                          "#F1D20B","#547FB6","#2D659E", "#c89dc8","#7acdc4", "#c78980",
                                          "#6A994E","#BC6C25","#3F37C9","#4895EF","#F72585","#4CC9F0","#F8961E","#90A959","#F9844A","#9970AB",
                                          "#F9C74F","#90BE6D","#43AA8B","#577590","#F8961E","#277DA1","#F94144","#4D908E","#F9844A","#9067BE"),
                           other_name = c("unknown","Unassigned")
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Missing ggplot2; run the first-time environment setup in Analysis_pipeline.Rmd.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Missing dplyr; run the first-time environment setup in Analysis_pipeline.Rmd.")
  }
  
  library(ggplot2)
  library(dplyr)
  
  if (!is.numeric(limit_n) || limit_n <= 0 || limit_n != as.integer(limit_n)) {
    stop("limit_n必须是正整数")
  }
  limit_n <- as.integer(limit_n)
  freq_table <- dataframe %>%
    count(!!sym(data_col)) %>%
    mutate(percentage = n / sum(n) * 100,
           label_percent = paste0(!!sym(data_col)," (",round(percentage, 2), "%)")) %>% 
    arrange(desc(percentage))
  
  names(freq_table)[1] = "names"
  other_rows <- freq_table[freq_table$names %in% other_name, ]
  keep_rows <- freq_table[!freq_table$names %in% other_name, ]
  freq_table <- rbind(keep_rows, other_rows)
  if (nrow(freq_table) > limit_n) {
    # 提取前n个类别
    top_n <- freq_table[1:limit_n, ]
    # 计算剩余类别的总和
    other_sum <- freq_table[(limit_n + 1):nrow(freq_table), ] %>%
      summarise(
        n = sum(n),
        percentage = sum(percentage),
        label_percent = paste0("other (", round(sum(percentage), 2), "%)"),
        names = "other"
      )
    # 合并数据
    freq_table <- bind_rows(top_n, other_sum)
  }
  
  freq_table$names <- factor(freq_table$names, levels = c(as.character(freq_table$names), "other")[1:nrow(freq_table)])
  freq_table$label_percent <- factor(freq_table$label_percent, levels = c(as.character(freq_table$label_percent))[1:nrow(freq_table)])
  p <- ggplot(freq_table, aes(x = "", y = n, fill = names)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(fill = legend_t) +
    theme(legend.position = "none")
  
  if (label) {
    if (pre) {
      p <- p + 
        geom_text(aes(label = label_percent), 
                  position = position_stack(vjust = 0.5), size = 5, colour = "black")
    } else {
      p <- p + 
        geom_text(aes(label = names), 
                  position = position_stack(vjust = 0.5), size = 5, colour = "black")
    }
  } else {
    if (pre) {
      p <- ggplot(freq_table, aes(x = "", y = n, fill = label_percent)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void() +
        labs(fill = legend_t) +
        theme(legend.position = "right")
    } else {
      p <- ggplot(freq_table, aes(x = "", y = n, fill = names)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void() +
        labs(fill = legend_t) +
        theme(legend.position = "right")
    }
  }
  
  # 颜色列表（扩展颜色数量，避免颜色不足）
  color_list = color_list
  
  # 确保颜色数量足够（循环使用颜色列表）
  color_values <- rep(color_list, length.out = nrow(freq_table))
  if (nrow(freq_table)>=limit_n){
    color_values[nrow(freq_table)]="#eeeeee"
  }
  p <- p + scale_fill_manual(values = color_values)

  return(p)
}

