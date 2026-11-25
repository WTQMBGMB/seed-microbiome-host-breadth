# =============================================================================
# Preprocessing helpers / 预处理函数
#
# CSS normalisation, host-tree derivation, alpha/beta caches, trait-table
# integration, and materialisation of read-only server results.
# CSS 标准化、宿主树派生、alpha/β 缓存、性状表整合，以及服务器结果落位。
# =============================================================================

# ---- CSS normalisation / CSS 标准化 -----------------------------------------
# Cumulative sum scaling (metagenomeSeq); input must be integer counts.
# 输入必须为非负整数计数；输出为样本×ASV 的标准化丰度表。
css_normalize <- function(abund_path, out_path, first_col = "ASV_ID",
                          norm_factor_path = NULL, verbose = TRUE) {
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("Package metagenomeSeq is required / 需要 R 包 metagenomeSeq")
  }
  if (verbose) message("[CSS] reading ", basename(abund_path))
  tab <- data.table::fread(abund_path, sep = "\t", header = TRUE)
  ids <- as.character(tab[[1]])
  otu <- data.matrix(tab[, -1, with = FALSE])
  rownames(otu) <- ids
  rm(tab)

  bad <- which(!vapply(seq_len(ncol(otu)), function(j) {
    x <- otu[, j]
    all(is.finite(x)) && all(x >= 0) && all(x == floor(x))
  }, logical(1)))
  if (length(bad)) {
    stop("CSS requires integer counts; problematic columns: ",
         paste(colnames(otu)[head(bad, 10)], collapse = ", "))
  }

  mre <- metagenomeSeq::newMRexperiment(otu)
  p <- metagenomeSeq::cumNormStatFast(mre)
  if (verbose) message("[CSS] quantile p = ", format(p, digits = 4))
  mre <- metagenomeSeq::cumNorm(mre, p = p)
  css <- metagenomeSeq::MRcounts(mre, norm = TRUE, log = FALSE)

  out <- data.table::as.data.table(css, keep.rownames = TRUE)
  data.table::setnames(out, "rn", first_col)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out, out_path, sep = "\t", quote = FALSE)

  if (!is.null(norm_factor_path)) {
    nf <- metagenomeSeq::normFactors(mre)
    nf_tab <- data.frame(sample = names(nf), css_quantile = as.numeric(nf))
    dir.create(dirname(norm_factor_path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(nf_tab, norm_factor_path, sep = "\t", quote = FALSE)
  }
  invisible(out)
}

# Build every CSS table that is not already available.
# 生成缺失的 CSS 表；默认随仓库提供的表存在时不重算。
ensure_css <- function(force = FALSE) {
  merged <- list(
    `16s` = c("filtered_merge/16S/Merged_16S_table_filtered.txt", "16s"),
    its   = c("filtered_merge/ITS/Merged_ITS_table_filtered.txt", "its")
  )
  for (ds in names(merged)) {
    out <- file.path(PREPARED_DIR, "standard_merge", paste0("css_", ds, "_asvtab.txt"))
    nf  <- file.path(PREPARED_DIR, "standard_merge", paste0("css_", ds, "_normFactors.txt"))
    shipped <- file.path(DATA_DIR, "standard_merge", paste0("css_", ds, "_asvtab.txt"))
    if (force || !file.exists(shipped)) {
      css_normalize(file.path(DATA_DIR, merged[[ds]][1]), out, norm_factor_path = nf)
    }
  }

  regions <- c("its1", "its2", "v3_v4", "v4", "v5_v6", "v5_v7")
  for (reg in regions) {
    upper <- toupper(reg)
    inp <- file.path(DATA_DIR, "filtered_region", upper,
                     paste0(upper, "_merged_asvtab_filtered.txt"))
    out <- file.path(PREPARED_DIR, "standard_region", reg,
                     paste0("css_", reg, "_asvtab.txt"))
    nf  <- file.path(PREPARED_DIR, "standard_region", reg,
                     paste0("css_", reg, "_normFactors.txt"))
    shipped <- file.path(DATA_DIR, "standard_region", reg,
                         paste0("css_", reg, "_asvtab.txt"))
    if (force || !file.exists(shipped)) {
      css_normalize(inp, out, norm_factor_path = nf)
    }
  }
  invisible(TRUE)
}

# ---- Host trees / 宿主树 ----------------------------------------------------
# Derive region and genus trees from the TimeTree export and metadata.
# 从 TimeTree 导出树与 metadata 派生区段树和属级树。
build_host_trees <- function(out_dir = file.path(PREPARED_DIR, "plant_tree")) {
  metadata_16s <- data.table::fread(meta_path("16S"), sep = "\t")
  metadata_its <- data.table::fread(meta_path("ITS"), sep = "\t")
  tax <- unique(rbind(metadata_16s[, 2:9], metadata_its[, 2:9]))
  plant <- unique(tax$timetree)
  stopifnot(!anyNA(plant), length(plant) > 0L)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(plant, file.path(out_dir, "species_all.txt"))

  tree <- ape::read.tree(file.path(DATA_DIR, "plant_tree", "species_all.nwk"))
  tree <- ape::keep.tip(tree, intersect(tree$tip.label, plant))
  stopifnot(setequal(tree$tip.label, plant))
  tree$tip.label <- tax$species[match(tree$tip.label, tax$timetree)]
  ape::write.tree(tree, file.path(out_dir, "plant_species_tree_all.nwk"))

  v4 <- ape::keep.tip(tree, unique(metadata_16s$species[metadata_16s$Region == "V4"]))
  v57 <- ape::keep.tip(tree, unique(metadata_16s$species[metadata_16s$Region == "V5_V7"]))
  its1 <- ape::keep.tip(tree, unique(metadata_its$species[metadata_its$Region == "ITS1"]))
  ape::write.tree(v4, file.path(out_dir, "plant_species_tree_v4.nwk"))
  ape::write.tree(v57, file.path(out_dir, "plant_species_tree_v5_v7.nwk"))
  ape::write.tree(its1, file.path(out_dir, "plant_species_tree_its1.nwk"))

  genus <- tax[!duplicated(tax$Genus)]
  gt <- ape::keep.tip(tree, intersect(tree$tip.label, genus$species))
  gt$tip.label <- genus$Genus[match(gt$tip.label, genus$species)]
  stopifnot(!anyNA(gt$tip.label))
  ape::write.tree(gt, file.path(out_dir, "plant_Genus_tree_all.nwk"))
  invisible(out_dir)
}

ensure_host_trees <- function() {
  need <- c("plant_species_tree_all.nwk", "plant_species_tree_v4.nwk",
            "plant_species_tree_v5_v7.nwk", "plant_species_tree_its1.nwk",
            "plant_Genus_tree_all.nwk")
  have <- file.exists(file.path(DATA_DIR, "plant_tree", need))
  if (all(have)) return(invisible(TRUE))
  build_host_trees()
  invisible(TRUE)
}

# ---- Per-region tables / 区段表 ---------------------------------------------
REGION_CFG <- list(
  v3_v4 = list(dataset = "16S", filt = "V3_V4"),
  v4    = list(dataset = "16S", filt = "V4"),
  v5_v6 = list(dataset = "16S", filt = "V5_V6"),
  v5_v7 = list(dataset = "16S", filt = "V5_V7"),
  its1  = list(dataset = "ITS", filt = "ITS1"),
  its2  = list(dataset = "ITS", filt = "ITS2")
)

# Load a region's CSS table and align it to metadata SampleID.
# 读取区段 CSS 表并按 metadata SampleID 对齐。
load_region_otu <- function(reg) {
  cfg <- REGION_CFG[[reg]]
  raw <- data.table::fread(
    css_region_path(reg, paste0("css_", reg, "_asvtab.txt")),
    sep = "\t", header = TRUE)
  mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(mat) <- as.character(raw[[1]])
  otu <- t(mat)
  rm(raw, mat)
  meta_all <- data.table::fread(meta_path(cfg$dataset), sep = "\t")
  meta <- as.data.frame(meta_all[meta_all$Region == cfg$filt, ])
  common <- intersect(rownames(otu), meta$SampleID)
  if (!length(common)) {
    fixed <- gsub("^X", "", rownames(otu))
    common <- intersect(fixed, meta$SampleID)
    if (length(common)) rownames(otu) <- fixed
  }
  otu <- otu[common, , drop = FALSE]
  meta <- meta[match(common, meta$SampleID), , drop = FALSE]
  zero <- colSums(otu) <= 0
  if (any(zero)) otu <- otu[, !zero, drop = FALSE]
  list(otu = otu, meta = meta)
}

calc_alpha_div <- function(otu) {
  if (!is.matrix(otu)) otu <- as.matrix(otu)
  shannon <- vegan::diversity(otu, index = "shannon", MARGIN = 1)
  simpson <- vegan::diversity(otu, index = "simpson", MARGIN = 1)
  observed <- rowSums(otu > 0)
  pielou <- ifelse(observed > 1, shannon / log(observed), NA_real_)
  data.frame(SampleID = rownames(otu), Shannon = shannon, Simpson = simpson,
             Observed = observed, Pielou = pielou, stringsAsFactors = FALSE)
}

# Bray-Curtis via L1 distance / row sums (identical to vegan::vegdist).
# 用 L1 距离与行和计算 Bray-Curtis，与 vegan::vegdist 等价。
compute_bray <- function(otu) {
  rs <- rowSums(otu)
  dist_mat <- Rfast::Dist(otu, method = "manhattan", parallel = TRUE,
                          result = "matrix")
  bray <- dist_mat / outer(rs, rs, "+")
  rownames(bray) <- colnames(bray) <- rownames(otu)
  diag(bray) <- 0
  stopifnot(all(is.finite(bray)))
  as.dist(bray)
}

# Build alpha and Bray caches for the three manuscript regions.
# 为三个稿件区段生成 alpha 与 Bray 缓存。
ensure_alpha_beta <- function(force = FALSE) {
  out_dir <- file.path(PREPARED_DIR, "precomputed")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (reg in c("v4", "v5_v7", "its1")) {
    alpha_f <- file.path(out_dir, paste0("alpha_div_", reg, ".rds"))
    bray_f  <- file.path(out_dir, paste0("bray_dist_", reg, ".rds"))
    alpha_in <- file.path(DATA_DIR, "precomputed", paste0("alpha_div_", reg, ".rds"))
    bray_in  <- server_path("ch2", "beta", paste0("bray_dist_", reg, ".rds"))
    if (!force && file.exists(alpha_in) && file.exists(bray_in)) next
    d <- load_region_otu(reg)
    if (force || !file.exists(alpha_f)) {
      ad <- calc_alpha_div(d$otu)
      ad <- merge(ad, d$meta[, c("SampleID", "species", "Order", "Family", "Genus")],
                  by = "SampleID", all.x = TRUE)
      saveRDS(ad, alpha_f)
    }
    if (force || !file.exists(bray_f)) {
      saveRDS(compute_bray(d$otu), bray_f)
    }
  }
  invisible(TRUE)
}

# ---- Trait table / 性状表整合 -----------------------------------------------
# Reproduce the final integration used by Chapter 5: BugBase, FAPROTAX2 and
# metaTraits annotations are matched to the Chapter 4 ASV universe.
# 按 Ch5 最终口径整合 BugBase、FAPROTAX2 与 metaTraits 注释；
# GOS 的 n_hosts/mean_dist 等列不被 Ch5 使用，故不重建。
build_trait_table <- function(asv_ids) {
  tax <- data.table::fread(
    file.path(DATA_DIR, "filtered_merge", "16S", "Merged_16S_taxonomy_filtered.txt"),
    sep = "\t", check.names = FALSE)
  tax <- tax[taxid %in% asv_ids]

  # BugBase: derive the two categorical summaries used downstream.
  bug <- data.table::fread(
    file.path(DATA_DIR, "function_prediction", "bugbase_result.txt"),
    sep = "\t", check.names = FALSE)
  bug[, Aerobic_type := data.table::fifelse(
    Aerobic == "TRUE", "Aerobic",
    data.table::fifelse(Anaerobic == "TRUE", "Anaerobic",
      data.table::fifelse(Facultatively_Anaerobic == "TRUE",
                          "Facultatively_Anaerobic", "Unknown")))]
  bug[, Gram_type := data.table::fifelse(
    Gram_Positive == "TRUE", "Positive",
    data.table::fifelse(Gram_Negative == "TRUE", "Negative", "Unknown"))]
  bug_ids <- as.character(bug$taxid)
  bug <- bug[, setdiff(names(bug), "taxid"), with = FALSE]
  data.table::setnames(bug, paste0("bugbase_", names(bug)))

  # FAPROTAX2: screen the full prediction table, then attach target ASVs.
  fp <- data.table::fread(
    file.path(DATA_DIR, "function_prediction", "faprotax_result.txt"),
    sep = "\t", check.names = FALSE, na.strings = "NA")
  fp_ids <- as.character(fp[[1]])
  fp <- fp[, -1, with = FALSE]
  keep <- vapply(fp, function(v) {
    known <- sum(!is.na(v))
    sum(is.na(v)) / length(v) < .5 && known >= 100 &&
      sum(v == TRUE, na.rm = TRUE) / known > .01 &&
      length(unique(stats::na.omit(v))) == 2
  }, logical(1))
  fp <- fp[, keep, with = FALSE]
  data.table::setnames(fp, paste0("faprotax_", names(fp)))

  # metaTraits: boolean and continuous matrices are screened separately.
  mb <- data.table::fread(
    file.path(DATA_DIR, "function_prediction", "asv_trait_matrix_clean_boolean.tsv"),
    sep = "\t", check.names = FALSE)
  mb_ids <- as.character(mb[[1]])
  mb <- mb[, -1, with = FALSE]
  keep_b <- vapply(mb, function(x) {
    mean(is.na(x)) < 0.5 && mean(x == 1, na.rm = TRUE) > 0.05
  }, logical(1))
  mb <- mb[, keep_b, with = FALSE]
  mb <- as.data.frame(lapply(mb, function(x) as.logical(x)))
  names(mb) <- paste0("metatraits_", names(mb))

  mc <- data.table::fread(
    file.path(DATA_DIR, "function_prediction", "asv_trait_matrix_clean_continuous.tsv"),
    sep = "\t", check.names = FALSE)
  mc_ids <- as.character(mc[[1]])
  mc <- mc[, -1, with = FALSE]
  keep_c <- vapply(mc, function(x) mean(is.na(x)) < 0.5, logical(1))
  mc <- mc[, keep_c, with = FALSE]
  data.table::setnames(mc, paste0("metatraits_", names(mc)))

  # Preserve the taxonomy row order and fill annotations by match().
  idx_bug <- match(tax$taxid, bug_ids)
  idx_fp  <- match(tax$taxid, fp_ids)
  idx_mb  <- match(tax$taxid, mb_ids)
  idx_mc  <- match(tax$taxid, mc_ids)

  bind <- function(base, mat, idx) {
    if (!ncol(mat)) return(base)
    out <- mat[idx, , drop = FALSE]
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    cbind(as.data.frame(base), out)
  }
  out <- bind(tax, bug, idx_bug)
  out <- bind(out, fp, idx_fp)
  out <- bind(out, mb, idx_mb)
  out <- bind(out, mc, idx_mc)
  data.table::as.data.table(out)
}

# ---- Server assets / 服务器结果落位 -----------------------------------------
# Copy read-only server results into the result layout expected by chapters.
# 将只读服务器结果复制到各章读取的 result 路径。
materialize_server_assets <- function(strict = TRUE) {
  copy_one <- function(src, dst, required = TRUE) {
    if (!file.exists(src)) {
      if (required && strict) {
        stop("Missing server result: ", src,
             "\nProvide it under server_result/ or rerun the documented server script.")
      }
      return(invisible(FALSE))
    }
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dst, overwrite = TRUE)
    invisible(TRUE)
  }
  copy_tree <- function(src, dst, required = FALSE) {
    if (!dir.exists(src)) {
      if (required && strict) stop("Missing server directory: ", src)
      return(invisible(FALSE))
    }
    files <- list.files(src, recursive = TRUE, full.names = TRUE)
    for (f in files) {
      rel <- sub(paste0("^", gsub("\\\\", "/", src), "/?"), "", gsub("\\\\", "/", f))
      copy_one(f, file.path(dst, rel), required = required)
    }
    invisible(TRUE)
  }

  copy_one(server_path("ch1", "host_overview_tree.pdf"),
           file.path(ch_dir("ch1_dataset", "figures"), "host_overview_tree.pdf"))
  copy_one(server_path("ch2", "mixed_pool_adonis.rds"),
           file.path(rds_dir("ch2_host_effect"), "mixed_pool_adonis.rds"))
  for (f in c("partA_units.rds", "partB_pools.rds")) {
    copy_one(server_path("ch2", "adonis", f), file.path(rds_dir("ch2_host_effect"), f))
  }
  copy_tree(server_path("ch2", "mixed_pool_raw"),
            file.path(int_dir("ch2_host_effect"), "mixed_pool"))
  copy_one(server_path("ch3", "asv_ref_silva.nwk"),
           file.path(tree_dir("ch3_phylogeny"), "asv_ref_silva.nwk"))
  copy_one(server_path("ch3", "boot_lambda_samples.rds"),
           file.path(rds_dir("ch3_phylogeny"), "boot_lambda_samples.rds"))
  copy_one(server_path("ch4", "gos_hostn_pos_16s_midpoint.nwk"),
           file.path(tree_dir("ch4_generality"), "gos_hostn_pos_16s_midpoint.nwk"))
  for (f in c("Tab_4_tree_phylosig.tsv", "Tab_4_tree_generalist_D.tsv")) {
    copy_one(server_path("ch4", f), file.path(ch_dir("ch4_generality", "tables"), f))
  }
  copy_tree(server_path("ch4", "server_raw"),
            file.path(tree_dir("ch4_generality"), "server_raw"))
  copy_tree(server_path("ch4", "itol_annotation"),
            file.path(tree_dir("ch4_generality"), "itol_annotation"))
  invisible(TRUE)
}
