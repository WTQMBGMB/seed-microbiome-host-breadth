# ============================================================================
# Mixed-pool PCoA + adonis2 (3 effects: Host Genus, source, Other)
#
# Other = interaction(Country, IF_SD, PlatformN, drop = TRUE)
# Model: bray ~ Genus + source + Other, by = "margin"
#
# Run on server (Linux, vegan 2.7.x):
#   Rscript _mixed_pool_adonis_server.R <work_dir> <nperm> <ncore>
#
# Inputs in <work_dir>:
#   bray_dist_{v4,v5_v7,its1}.rds
#   metadata_16S.txt, metadata_ITS.txt
# Outputs in <work_dir>/out:
#   mixed_pool_adonis_<reg>.rds   (adonis object + design/labels)
#   mixed_pool_adonis_3regions.tsv
#   mixed_pool_design_<reg>.tsv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(vegan)
})
options(stringsAsFactors = FALSE)

a <- commandArgs(trailingOnly = TRUE)
work <- a[1]
nperm <- as.integer(a[2])
ncore <- as.integer(a[3])
outdir <- file.path(work, "out")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
CNTRY_FIX <- c(
  "Cote d Ivoire" = "Côte d'Ivoire", "Italia" = "Italy",
  "Korea" = "South Korea", "USA" = "United States of America",
  "United States" = "United States of America", "Czech Republic" = "Czechia",
  "England" = "United Kingdom", "Republic of the Congo" = "Congo",
  "Democratic Republic of the Congo" = "Dem. Rep. Congo")
fix_country <- function(x) {
  x <- trimws(x); x[is.na(x) | x == ""] <- "Unknown"
  hit <- x %in% names(CNTRY_FIX); x[hit] <- unname(CNTRY_FIX[x[hit]]); x
}
read_meta <- function(ds) {
  m <- fread(file.path(work, sprintf("metadata_%s.txt", ds)),
             sep = "\t", na.strings = c("", "NA"), quote = "")
  m[, `:=`(dataset = ds, PlatformN = normalize_platform(Platform))]
  m[, Country := fix_country(Country)]
  m[, extract := trimws(extract)]
  m
}

REG <- c(v4 = "16S", v5_v7 = "16S", its1 = "ITS")
all_tab <- list()
for (reg in names(REG)) {
  cat(sprintf("\n===== [%s] %s =====\n", reg, format(Sys.time())))
  bd  <- readRDS(file.path(work, sprintf("bray_dist_%s.rds", reg)))
  labs <- attr(bd, "Labels")
  m <- read_meta(REG[[reg]])[match(labs, SampleID)]
  stopifnot(nrow(m) == length(labs), all(!is.na(m$SampleID)))

  m[, Other := interaction(Country, IF_SD, PlatformN, drop = TRUE, sep = "|")]
  keep <- complete.cases(m[, .(Genus, source, Country, IF_SD, PlatformN)])
  cat(sprintf("samples=%d complete=%d dropped=%d\n",
              length(labs), sum(keep), sum(!keep)))
  m <- m[keep]; labs <- labs[keep]

  dm <- as.matrix(bd)[keep, keep, drop = FALSE]
  D  <- as.dist(dm)
  design <- data.frame(
    Genus  = droplevels(factor(m$Genus)),
    source = droplevels(factor(m$source)),
    Other  = droplevels(factor(m$Other))
  )
  rownames(design) <- labs
  cat(sprintf("levels: Genus=%d source=%d Other=%d | n=%d\n",
              nlevels(design$Genus), nlevels(design$source),
              nlevels(design$Other), nrow(design)))
  cat(sprintf("rank checks: full=%d genus=%d source=%d other=%d\n",
              qr(model.matrix(~ Genus + source + Other, design))$rank,
              qr(model.matrix(~ Genus, design))$rank,
              qr(model.matrix(~ source, design))$rank,
              qr(model.matrix(~ Other, design))$rank))

  set.seed(20261120)
  t0 <- Sys.time()
  ad <- adonis2(D ~ Genus + source + Other, data = design,
                by = "margin", permutations = nperm, parallel = ncore)
  cat(sprintf("adonis2 done in %.1f min with %d cores, %d permutations\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              ncore, nperm))
  print(ad)

  saveRDS(list(ad = ad, design = design, labels = labs, region = reg,
               nperm = nperm, ncore = ncore),
          file.path(outdir, sprintf("mixed_pool_adonis_%s.rds", reg)))
  fwrite(cbind(labels = labs, design),
         file.path(outdir, sprintf("mixed_pool_design_%s.tsv", reg)),
         sep = "\t")

  at <- as.data.frame(ad); at$term <- rownames(at)
  at <- at[!at$term %in% c("Residual", "Total"), ]
  all_tab[[reg]] <- data.table(
    Region = reg, term = at$term, Df = at$Df, SumOfSqs = at$SumOfSqs,
    R2 = at$R2, F = at$F, P = at$`Pr(>F)`, n = nrow(design),
    n_genus = nlevels(design$Genus), n_source = nlevels(design$source),
    n_other = nlevels(design$Other), nperm = nperm, ncore = ncore)
}
tab <- rbindlist(all_tab)
fwrite(tab, file.path(outdir, "mixed_pool_adonis_3regions.tsv"), sep = "\t")
cat("\nSaved outputs to", outdir, "\n")
print(tab)
