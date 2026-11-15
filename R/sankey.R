# Sankey rendering restored from the archived plotting scripts
# 按归档绘图脚本恢复桑基图：GOS 五列图与 KO 两组图

# 平滑色带 / Smooth ribbon between two full node blocks.
sankey_ribbon <- function(x1, x2, y1min, y1max, y2min, y2max,
                          id, taxon, type, n = 80L) {
  tt <- seq(0, 1, length.out = n)
  ss <- tt * tt * (3 - 2 * tt)
  data.table::data.table(
    x = x1 + (x2 - x1) * tt,
    ymin = y1min + (y2min - y1min) * ss,
    ymax = y1max + (y2max - y1max) * ss,
    id = id, Taxon = taxon, type = type)
}

# GOS 五列桑基：节点高度、列顺序和色带均按归档表绘制。
# GOS five-column Sankey: node heights, column order and ribbons follow
# the archived tables.
draw_gos_sankey <- function(nodes, links, taxon_col, node_colours,
                            title = NULL, subtitle = NULL,
                            show_legend = TRUE, legend_title = NULL,
                            node_border = "white", node_width = 0.32,
                            panel_margin = c(6, 10, 6, 10)) {
  nodes <- data.table::copy(data.table::as.data.table(nodes))
  links <- data.table::copy(data.table::as.data.table(links))
  if (taxon_col %in% names(links) && taxon_col != "Taxon") {
    data.table::setnames(links, taxon_col, "Taxon")
  }
  stopifnot(all(c("Category", "Node", "n", "frac", "ymin", "ymax", "x") %in%
                  names(nodes)))
  stopifnot(all(c("Pair", "From", "To", "Taxon", "type") %in% names(links)))

  categories <- nodes[order(x), unique(Category)]
  nodes[, key := paste(Category, Node, sep = "||")]
  links[, `:=`(from_cat = sub("->.*$", "", Pair), to_cat = sub("^.*->", "", Pair))]
  li <- match(paste(links$from_cat, links$From, sep = "||"), nodes$key)
  ri <- match(paste(links$to_cat, links$To, sep = "||"), nodes$key)
  stopifnot(!anyNA(li), !anyNA(ri))

  ribbons <- data.table::rbindlist(lapply(seq_len(nrow(links)), function(i) {
    sankey_ribbon(nodes$x[li[i]] + node_width / 2,
                  nodes$x[ri[i]] - node_width / 2,
                  nodes$ymin[li[i]], nodes$ymax[li[i]],
                  nodes$ymin[ri[i]], nodes$ymax[ri[i]],
                  i, links$Taxon[i], links$type[i])
  }))
  other_fill <- unname(node_colours["Other"])
  if (is.na(other_fill)) other_fill <- "#D9D9D9"
  col_tot <- nodes[, .(n_total = sum(n)), by = Category]
  col_tot[, x := match(Category, categories)]

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ribbons[type %in% c("dropout", "emerge")],
      ggplot2::aes(x = x, ymin = ymin, ymax = ymax, group = id, fill = Taxon),
      alpha = 0.18, colour = NA) +
    ggplot2::geom_ribbon(
      data = ribbons[type == "shared"],
      ggplot2::aes(x = x, ymin = ymin, ymax = ymax, group = id, fill = Taxon),
      alpha = 0.55, colour = NA) +
    ggplot2::geom_ribbon(
      data = ribbons[type == "other"],
      ggplot2::aes(x = x, ymin = ymin, ymax = ymax, group = id),
      fill = other_fill, alpha = 0.45, colour = NA) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = x - node_width / 2, xmax = x + node_width / 2,
                   ymin = ymin, ymax = ymax, fill = Node),
      colour = node_border, linewidth = 0.35) +
    ggplot2::geom_text(
      data = nodes[frac >= 0.02],
      ggplot2::aes(x = x, y = (ymin + ymax) / 2, label = Node),
      size = 2.5, fontface = "bold", colour = "grey10") +
    ggplot2::geom_text(
      data = col_tot,
      ggplot2::aes(x = x, y = 1.025,
                   label = sprintf("%s\n(n=%s)", Category,
                                   format(n_total, big.mark = ","))),
      inherit.aes = FALSE, size = 3.0, fontface = "bold",
      lineheight = 0.95, vjust = 0) +
    ggplot2::scale_fill_manual(values = node_colours, name = legend_title,
                               guide = if (show_legend) "legend" else "none") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(categories), labels = categories,
      limits = c(0.55, length(categories) + 0.45), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.12), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      legend.key.size = grid::unit(0.32, "cm"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.5, hjust = 0),
      plot.margin = ggplot2::margin(panel_margin[1], panel_margin[2],
                                    panel_margin[3], panel_margin[4]),
      axis.text.x = ggplot2::element_text(face = "bold", size = 9,
                                          margin = ggplot2::margin(t = 4)))
}

# KO 两组桑基：六列分组重叠；shared / dropout / emerge 分别按归档透明度绘制。
# KO grouped Sankey: overlapping six groups; shared / dropout / emerge use archived alphas.
draw_ko_sankey <- function(nodes, links, absolute = FALSE,
                           show_legend = TRUE, legend_title = "KEGG level2",
                           pgpg = FALSE) {
  nodes <- data.table::copy(data.table::as.data.table(nodes))
  links <- data.table::copy(data.table::as.data.table(links))
  stopifnot(all(c("Group", "Category", "n", "frac", "ymin", "ymax", "x") %in% names(nodes)))
  stopifnot(all(c("Pair", "From", "To", "Category", "type") %in% names(links)))

  groups <- nodes[order(x), unique(Group)]
  cat_row <- unique(as.character(nodes$Category))
  top_cats <- rev(setdiff(cat_row, c("Other", "Unknown")))
  cat_levels <- c(top_cats, "Other", if ("Unknown" %in% cat_row) "Unknown")
  nodes[, Category := factor(Category, levels = cat_levels)]
  links[, Category := factor(Category, levels = cat_levels)]

  ko_pal <- c("#FCB366", "#FE816F", "#90D7C7", "#93BCD0", "#AF3F4F",
              "#E5C28C", "#CA8BB6", "#1A6C92", "#EAB46E", "#8CC3BD",
              "#8F9EC5", "#6C77B5", "#138A8E", "#4B9D77", "#8E60A8")
  cols <- stats::setNames(ko_pal[seq_along(top_cats)], top_cats)
  cols["Other"] <- "#D9D9D9"
  if ("Unknown" %in% cat_levels) cols["Unknown"] <- "#BDBDBD"

  pgpg_map <- c(
    "COLONIZATION-PLANT_DERIVED_SUBSTRATE_USAGE" = "Plant-derived substrate utilization",
    "CE-BACTERIAL_FITNESS"                       = "Bacterial fitness",
    "NEUTRALIZING_ABIOTIC_STRESS"                = "Abiotic stress alleviation",
    "NEUTRALIZING_BIOTIC_STRESS"                 = "Biotic stress alleviation",
    "XENOBIOTICS_BIODEGRADATION"                 = "Xenobiotic degradation",
    "PHOSPHATE_SOLUBILIZATION"                   = "Phosphate solubilization",
    "NITROGEN_ACQUISITION"                       = "Nitrogen acquisition",
    "CE-QUORUM_SENSING_RESPONSE|BIOFILM_FORMATION" = "Quorum sensing and biofilm formation",
    "COLONIZATION-SURFACE_ATTACHMENT"            = "Surface attachment",
    "CE-BACTERIAL_SECRETION"                     = "Bacterial secretion",
    "PLANT_VITAMIN_PRODUCTION"                   = "Plant-beneficial vitamin production",
    "COLONIZATION-MOTILITY|CHEMOTAXIS"           = "Motility and chemotaxis",
    "HEAVY_METAL_DETOXIFICATION"                 = "Heavy metal detoxification",
    "CE-CELL_ENVELOPE_REMODELLING"               = "Cell envelope remodeling",
    "IRON_ACQUISITION"                           = "Iron acquisition")
  wrap_lab <- function(x, width = 44) {
    vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
           character(1))
  }
  legend_labels <- cat_levels
  if (pgpg) {
    mapped <- unname(pgpg_map[cat_levels])
    legend_labels <- unname(ifelse(is.na(mapped), cat_levels, mapped))
  }

  nodes[, key := paste(Group, Category, sep = "||")]
  links[, `:=`(from_cat = sub("->.*$", "", Pair), to_cat = sub("^.*->", "", Pair))]
  li <- match(paste(links$from_cat, links$From, sep = "||"), nodes$key)
  ri <- match(paste(links$to_cat, links$To, sep = "||"), nodes$key)
  stopifnot(!anyNA(li), !anyNA(ri))

  ribbons <- data.table::rbindlist(lapply(seq_len(nrow(links)), function(i) {
    sankey_ribbon(nodes$x[li[i]] + 0.36 / 2,
                  nodes$x[ri[i]] - 0.36 / 2,
                  nodes$ymin[li[i]], nodes$ymax[li[i]],
                  nodes$ymin[ri[i]], nodes$ymax[ri[i]],
                  i, links$Category[i], links$type[i])
  }))
  ribbons[, Category := factor(Taxon, levels = cat_levels)]

  col_tot <- nodes[, .(n_total = sum(n), x = x[1]), by = Group]
  max_n <- max(col_tot$n_total)
  if (absolute) {
    col_tot[, y_head := n_total + max_n * 0.03]
    y_lim <- c(0, max_n * 1.12)
  } else {
    col_tot[, y_head := 1.03]
    y_lim <- c(0, 1.10)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ribbons[type == "shared"],
      ggplot2::aes(x = x, ymin = ymin, ymax = ymax, group = id, fill = Category),
      alpha = 0.55, colour = NA) +
    ggplot2::geom_ribbon(
      data = ribbons[type != "shared"],
      ggplot2::aes(x = x, ymin = ymin, ymax = ymax, group = id, fill = Category),
      alpha = 0.18, colour = NA) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = x - 0.36 / 2, xmax = x + 0.36 / 2,
                   ymin = ymin, ymax = ymax, fill = Category),
      colour = NA) +
    ggplot2::geom_text(
      data = col_tot,
      ggplot2::aes(x = x, y = y_head,
                   label = sprintf("n=%s", format(n_total, big.mark = ","))),
      inherit.aes = FALSE, size = 3.0, fontface = "bold", vjust = 0) +
    ggplot2::scale_fill_manual(values = cols, name = legend_title,
                               labels = unname(wrap_lab(legend_labels)),
                               guide = if (show_legend) "legend" else "none") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(groups), labels = groups,
      limits = c(0.55, length(groups) + 0.45), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      legend.key.size = grid::unit(0.32, "cm"),
      legend.text = ggplot2::element_text(size = 7, lineheight = 0.9),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      axis.text.x = ggplot2::element_text(face = "bold", size = 8,
                                          margin = ggplot2::margin(t = 4)))
  if (!absolute) {
    p <- p + ggplot2::geom_text(
      data = nodes[frac >= 0.02 & n > 0],
      ggplot2::aes(x = x, y = (ymin + ymax) / 2,
                   label = sprintf("%.1f%%", 100 * frac)),
      size = 2.8, colour = "grey10", fontface = "bold")
  }
  p
}
