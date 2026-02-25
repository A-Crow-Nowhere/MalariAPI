#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(ExomeDepth)
  library(GenomicRanges)
  library(Rsamtools)
}))

usage <- function() {
  cat("
Usage:
  run_exomedepth.R --sample <NAME> --test-bam <test.bam> --ref-bams <a.bam,b.bam,...>
                   --exons-bed <exons.bed> --out-tsv <calls.tsv> --out-bed <calls.bed>

Notes:
  - BAMs must be coordinate-sorted and indexed.
  - Exons BED must be: chrom  start0  end  [name]
  - This helper intentionally avoids rtracklayer to prevent XML/RCurl dependency issues.
\n", file=stderr())
}

arg_val <- function(args, key, default=NA_character_) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing value for", key))
  args[[i + 1]]
}

read_bed4 <- function(path) {
  # Try data.table for speed/robustness if present; else base R.
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(
      path,
      sep = "\t",
      header = FALSE,
      data.table = FALSE,
      showProgress = FALSE,
      fill = TRUE
    )
  } else {
    dt <- utils::read.table(
      path,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    )
  }

  if (ncol(dt) < 3) stop("BED must have >=3 columns: chrom start end")

  chrom  <- as.character(dt[[1]])
  start0 <- as.integer(dt[[2]])
  end    <- as.integer(dt[[3]])

  name <- rep(".", length(chrom))
  if (ncol(dt) >= 4) {
    name <- as.character(dt[[4]])
    name[is.na(name) | name == ""] <- "."
  }

  # Convert BED start0 -> 1-based start for ExomeDepth bed.frame
  start1 <- pmax(1L, start0 + 1L)

  data.frame(
    chromosome = chrom,
    start = start1,
    end = end,
    name = name,
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h","--help"))) {
  usage()
  quit(status = 0)
}

sample    <- arg_val(args, "--sample")
test_bam  <- arg_val(args, "--test-bam")
ref_csv   <- arg_val(args, "--ref-bams")
exons_bed <- arg_val(args, "--exons-bed")
out_tsv   <- arg_val(args, "--out-tsv")
out_bed   <- arg_val(args, "--out-bed")

if (any(is.na(c(sample, test_bam, ref_csv, exons_bed, out_tsv, out_bed)))) {
  usage()
  stop("Missing required arguments.")
}

ref_bams <- strsplit(ref_csv, ",", fixed = TRUE)[[1]]
ref_bams <- ref_bams[nzchar(ref_bams)]
if (length(ref_bams) < 1) stop("No reference BAMs parsed from --ref-bams")

if (!file.exists(test_bam)) stop(paste("test BAM not found:", test_bam))
if (!file.exists(exons_bed)) stop(paste("exons BED not found:", exons_bed))
for (b in ref_bams) if (!file.exists(b)) stop(paste("ref BAM not found:", b))

bed.frame <- read_bed4(exons_bed)
include_chr <- unique(bed.frame$chromosome)

message("[exomedepth] test=", test_bam)
message("[exomedepth] refs=", length(ref_bams))
message("[exomedepth] exons=", nrow(bed.frame))

message("[exomedepth] counting test BAM...")
test_counts <- ExomeDepth::getBamCounts(
  bed.frame = bed.frame[, c("chromosome","start","end")],
  bam.files = test_bam,
  include.chr = include_chr,
  referenceFasta = NULL
)$counts

message("[exomedepth] counting reference BAM(s)...")
ref_mat <- NULL
for (b in ref_bams) {
  message("[exomedepth]  ref: ", b)
  cts <- ExomeDepth::getBamCounts(
    bed.frame = bed.frame[, c("chromosome","start","end")],
    bam.files = b,
    include.chr = include_chr,
    referenceFasta = NULL
  )$counts
  ref_mat <- cbind(ref_mat, cts)
}

ref_counts <- rowSums(ref_mat)

ed <- new(
  "ExomeDepth",
  test = as.numeric(test_counts),
  reference = as.numeric(ref_counts),
  formula = "cbind(test, reference) ~ 1"
)

message("[exomedepth] calling CNVs...")
ed <- CallCNVs(
  x = ed,
  transition.probability = 1e-4,
  chromosome = bed.frame$chromosome,
  start = bed.frame$start,
  end = bed.frame$end,
  name = bed.frame$name
)

cnvs <- ed@CNV.calls

if (is.null(cnvs) || nrow(cnvs) == 0) {
  message("[exomedepth] no calls")
  write.table(data.frame(), file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("", file = out_bed)
  quit(status = 0)
}

write.table(cnvs, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

bed_out <- data.frame(
  chrom = cnvs$chromosome,
  start = pmax(0L, as.integer(cnvs$start) - 1L),
  end   = as.integer(cnvs$end),
  name  = paste0(
    cnvs$type,
    ";BF=", signif(cnvs$BF, 3),
    ";ratio=", signif(cnvs$read.ratio, 3),
    ";sample=", sample
  ),
  stringsAsFactors = FALSE
)

write.table(bed_out, file = out_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

message("[exomedepth] done")
