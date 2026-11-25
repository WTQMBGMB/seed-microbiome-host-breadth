# =============================================================================
# Shared helpers / 公共函数
#
# Metadata cleaning, fixed palettes, plotting theme, and output writers.
# Chapter-specific analysis stays in the chapter Rmd files.
# 元数据清洗、固定配色、绘图主题与输出函数；具体分析保留在章节 Rmd 中。
# =============================================================================

META_COLS <- c("SampleID", "Kingdom", "Phylum", "Class", "Order", "Family",
               "Genus", "species", "timetree", "bio_project", "file",
               "source", "Cite", "extract", "Country", "Region", "IF_SD",
               "Platform")

# Platform labels are harmonised before they enter any model.
# 平台标签先统一，再进入模型。
normalize_platform <- function(x) {
  x <- trimws(x)
  out <- rep("Unspecified", length(x))
  out[grepl("miseq",   x, ignore.case = TRUE)] <- "MiSeq"
  out[grepl("hiseq",   x, ignore.case = TRUE)] <- "HiSeq"
  out[grepl("novaseq", x, ignore.case = TRUE)] <- "NovaSeq"
  out[grepl("nextseq", x, ignore.case = TRUE)] <- "NextSeq"
  out[grepl("aviti",   x, ignore.case = TRUE)] <- "AVITI"
  out[grepl("454",     x)]                      <- "454"
  out[is.na(x) | x == ""]                       <- "Unspecified"
  out
}

# Country names follow the Natural Earth spelling used in the map.
# 国家名统一为 Natural Earth 拼写，供世界地图匹配。
CNTRY_FIX <- c(
  "Cote d Ivoire"                    = "Côte d'Ivoire",
  "Italia"                           = "Italy",
  "Korea"                            = "South Korea",
  "USA"                              = "United States of America",
  "United States"                    = "United States of America",
  "Czech Republic"                   = "Czechia",
  "England"                          = "United Kingdom",
  "Republic of the Congo"            = "Congo",
  "Democratic Republic of the Congo" = "Dem. Rep. Congo"
)

fix_country <- function(x) {
  x <- trimws(x)
  x[is.na(x) | x == ""] <- "Unknown"
  hit <- x %in% names(CNTRY_FIX)
  x[hit] <- unname(CNTRY_FIX[x[hit]])
  x
}

# Read one kingdom's metadata and enforce the 18-column contract.
# 读取单个 kingdom 的 metadata，并校验 18 列结构。
read_meta <- function(ds) {
  stopifnot(ds %in% c("16S", "ITS"))
  f <- meta_path(ds)
  m <- data.table::fread(f, sep = "\t", na.strings = c("", "NA"),
                         quote = "", data.table = TRUE)
  if (ncol(m) != 18L) {
    stop(sprintf("%s: expected 18 columns, found %d (check trailing tabs) /\n%s：应为 18 列，实际 %d 列，请检查行尾制表符",
                 f, ncol(m), f, ncol(m)), call. = FALSE)
  }
  miss <- setdiff(META_COLS, names(m))
  if (length(miss)) stop("metadata is missing columns: ", paste(miss, collapse = ", "))
  m[, `:=`(
    dataset    = ds,
    kingdom    = if (ds == "16S") "Bacteria" else "Fungi",
    PlatformN  = normalize_platform(Platform),
    CountryRaw = Country
  )]
  m[, Country := fix_country(Country)]
  m[, extract := trimws(extract)]
  m[]
}

read_all_meta <- function() {
  data.table::rbindlist(list(read_meta("16S"), read_meta("ITS")), fill = TRUE)
}

# ---- Fixed palettes / 固定配色 ----------------------------------------------
KING_COL <- c(Bacteria = "#3A5A40", Fungi = "#D97B29")

REG_COL <- c(
  V3_V4 = "#3A5A40", V4    = "#F2CC8F",
  V5_V6 = "#81B29A", V5_V7 = "#3D405B",
  ITS1  = "#F4A261", ITS2  = "#099396"
)

REG_LEVELS <- c("V3_V4", "V4", "V5_V6", "V5_V7", "ITS1", "ITS2")

SRC_LEVELS <- c("S1-57 Collected", "S58 Simonin 2022", "S59 In house")

COL_POS <- "#E67E22"
COL_NEG <- "#2E86C1"
COL_NS  <- "#9E9E9E"

PAL_QUAL <- c(
  "#3A5A40", "#D97B29", "#81B29A", "#3D405B", "#F4A261",
  "#099396", "#E7D6C2", "#547FB6", "#C5D8E7", "#BF6173",
  "#E07A5F", "#90DBF4", "#641E16", "#A7C957", "#D35400",
  "#8C6D31", "#C58B3C", "#2E86C1", "#9E9E9E", "#F2CC8F"
)

# Shared multi-level palette; colours are assigned by factor order.
# 通用多水平色盘，按因子顺序取色。
PAL_MULTI <- c(
  "#424949", "#3A5A40", "#8A926B", "#E4E2C9", "#78A495",
  "#ADD565", "#CFE1D3", "#7DCEA0", "#F1FAEE", "#1D3557",
  "#355070", "#613767", "#547FB6", "#099396", "#C5D8E7",
  "#A8DADC", "#90DBF4", "#0FBFBC", "#4E56AB", "#855B67",
  "#BF6173", "#FDE0E4", "#D4C0C2", "#995B70", "#846992",
  "#552379", "#E56B6F", "#CB4047", "#E76F51", "#BA4A00",
  "#F4A261", "#F9C74F", "#F2CC8F", "#E97DAE"
)

pal_multi <- function(n) {
  if (n > length(PAL_MULTI)) {
    stop(sprintf("PAL_MULTI has %d colours, requested %d", length(PAL_MULTI), n))
  }
  PAL_MULTI[seq_len(n)]
}

# ---- Font and theme / 字体与主题 --------------------------------------------
if (.Platform$OS.type == "windows" &&
    !"Arial" %in% names(grDevices::windowsFonts())) {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

theme_seed <- function(base_size = 8) {
  ggplot2::theme_bw(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      text              = ggplot2::element_text(face = "plain", size = base_size),
      axis.text         = ggplot2::element_text(face = "plain", size = base_size),
      axis.title        = ggplot2::element_text(face = "plain", size = base_size),
      legend.text       = ggplot2::element_text(face = "plain", size = base_size),
      legend.title      = ggplot2::element_text(face = "plain", size = base_size),
      strip.text        = ggplot2::element_text(face = "plain", size = base_size),
      plot.title        = ggplot2::element_blank(),
      plot.subtitle     = ggplot2::element_blank(),
      plot.tag          = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank()
    )
}

# Append a timestamped progress line to the chapter log.
# 向章节日志追加带时间戳的进度行。
make_progress <- function(chapter) {
  force(chapter)
  function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), ..., "\n",
        file = file.path(log_dir(chapter), paste0(chapter, "_progress.log")),
        append = TRUE)
  }
}

# ---- Output writers / 输出函数 ----------------------------------------------
# PDF + JPG in millimetres; the JPG device keeps UTF-8 symbols intact.
# 同时输出 PDF 与 JPG；JPG 使用 ragg 以保证 Unicode 符号正确。
save_fig <- function(p, chapter, stem, width_mm, height_mm, dpi = 300,
                     stem_check = TRUE) {
  if (stem_check) stopifnot(grepl("^Fig_[0-9]+_", stem))
  d <- ch_dir(chapter, "figures")
  ggplot2::ggsave(file.path(d, paste0(stem, ".pdf")), p,
                  width = width_mm, height = height_mm, units = "mm",
                  device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(d, paste0(stem, ".jpg")), p,
                  width = width_mm, height = height_mm, units = "mm",
                  device = ragg::agg_jpeg, dpi = dpi, quality = 95,
                  bg = "white")
  invisible(file.path(d, stem))
}

# Embed saved figures in the rendered document, in place of the code chunk.
# The JPG copy is used so the same call works for HTML and PDF output.
# 把已保存的图件嵌入渲染结果；使用 JPG，HTML 与 PDF 输出均适用。
# Paths are made relative to the rendered document so the HTML stays movable.
# 路径相对渲染文档，生成的 HTML 可整体移动。
embed_fig <- function(chapter, stems, width = "100%") {
  paths <- file.path(ch_dir(chapter, "figures"), paste0(stems, ".jpg"))
  miss <- !file.exists(paths)
  if (any(miss)) {
    stop("figure not found: ", paste(paths[miss], collapse = ", "), call. = FALSE)
  }
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  # render_one.R sets seed.html_dir; interactive knitting writes next to the Rmd.
  out_dir <- getOption("seed.html_dir", PROJECT_ROOT)
  paths <- xfun::relative_path(paths, normalizePath(out_dir, winslash = "/"),
                               error = FALSE)
  tags <- if (knitr::is_latex_output()) {
    sprintf("\\includegraphics[width=%s]{%s}", width, paths)
  } else {
    sprintf('<p><img src="%s" style="width:%s; max-width:100%%;" /></p>',
            paths, width)
  }
  knitr::asis_output(paste(tags, collapse = "\n"))
}

# Write a table atomically; byte-identical files are left untouched.
# 原子写出表格；内容未变的文件不重写。
write_tab <- function(dt, chapter, name) {
  stopifnot(grepl("^Tab_[0-9]+_", name))
  f <- file.path(ch_dir(chapter, "tables"), paste0(name, ".tsv"))
  tmp <- paste0(f, ".tmp-", Sys.getpid())
  data.table::fwrite(dt, tmp, sep = "\t", quote = FALSE, na = "NA")
  on.exit(unlink(tmp), add = TRUE)
  same <- file.exists(f) && file.size(f) == file.size(tmp) &&
    identical(readBin(f, "raw", file.size(f)),
              readBin(tmp, "raw", file.size(tmp)))
  if (same) return(invisible(f))
  ok <- isTRUE(file.rename(tmp, f)) || isTRUE(file.copy(tmp, f, overwrite = TRUE))
  if (!ok) stop(sprintf("could not overwrite %s (locked by another process?)", f))
  invisible(f)
}

write_int <- function(dt, chapter, name, shared = FALSE) {
  data.table::fwrite(dt, file.path(int_dir(chapter, shared), paste0(name, ".tsv")),
                     sep = "\t", quote = FALSE, na = "NA")
}

save_rds <- function(obj, chapter, name, shared = FALSE) {
  saveRDS(obj, file.path(rds_dir(chapter, shared), paste0(name, ".rds")))
}

read_int <- function(chapter, name, shared = FALSE) {
  data.table::fread(file.path(int_dir(chapter, shared), paste0(name, ".tsv")),
                    sep = "\t", na.strings = c("", "NA"), quote = "")
}

write_log <- function(chapter, name, lines) {
  writeLines(as.character(lines), file.path(log_dir(chapter), paste0(name, ".log")))
}

n_distinct_nonmissing <- function(x) {
  data.table::uniqueN(x[!is.na(x) & x != "" & x != "Unknown"])
}
