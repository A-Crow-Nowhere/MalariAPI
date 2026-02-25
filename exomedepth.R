#!/usr/bin/env Rscript
# ============================================================
# ExomeDepth batch CNV calling (cohort-as-reference)
# Compatible with archived ExomeDepth 1.1.16 (CRAN archive)
# - getBamCounts() returns a data.frame (NOT SummarizedExperiment)
# - PfDd2 contigs (NO "chr" prefix): include.chr = FALSE
#
# Output:
#   <out_dir>/<sample>.exomedepth.cnv_calls.tsv
#   <out_dir>/<sample>.exomedepth.object.rds
#   <out_dir>/exomedepth_manifest.tsv
# ============================================================

suppressPackageStartupMessages({
  library(VGAM)
  library(aod)
  library(ExomeDepth)
  library(Rsamtools)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
bam_dir  <- "/scratch/njb8sg/MalariAPI/scratch/clover/cannon/mimic_data/final_bams"
bed_file <- "~/MalariAPI/genomes/exons.bed"
out_dir  <- "/scratch/njb8sg/MalariAPI/scratch/clover/cannon/mimic_data/exomedepth"

# Optional: restrict chromosomes (NULL keeps all)
include_chr <- NULL

# Reference building
top_k_refs <- 6          # number of top-correlated refs to aggregate
min_ref_samples <- 2    # minimum other samples required
transition_prob <- 1e-4  # HMM transition probability (ExomeDepth default-ish)

# Filtering bins (removes dead exons where both test+ref are 0)
drop_all_zero_bins <- TRUE

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Logging helpers
# -----------------------------
say <- function(...) cat("[exomedepth] ", sprintf(...), "\n", sep = "")
die <- function(...) stop(sprintf(...), call. = FALSE)

# -----------------------------
# Helpers
# -----------------------------
read_targets <- function(bed_path) {
  # robust to tabs/spaces, avoids Windows CR
  bed <- read.table(
    bed_path,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  if (ncol(bed) < 3) die("BED must have >= 3 columns: chrom start end [name]")
  
  bed <- bed[, 1:min(4, ncol(bed)), drop = FALSE]
  if (ncol(bed) < 4) bed[[4]] <- paste0("target_", seq_len(nrow(bed)))
  
  colnames(bed)[1:4] <- c("chromosome", "start", "end", "name")
  
  bed$chromosome <- trimws(sub("\r$", "", bed$chromosome))
  bed$start <- as.integer(bed$start)
  bed$end   <- as.integer(bed$end)
  bed$name  <- trimws(sub("\r$", "", bed$name))
  
  bed
}

find_bam_index <- function(bam) {
  bai1 <- paste0(bam, ".bai")
  bai2 <- sub("\\.bam$", ".bai", bam)
  if (file.exists(bai1)) return(bai1)
  if (file.exists(bai2)) return(bai2)
  NA_character_
}

get_bam_seqs <- function(bam) {
  names(scanBamHeader(bam)[[1]]$targets)
}

make_reference_by_correlation <- function(test_counts, ref_counts_mat, k = 6) {
  # ref_counts_mat: bins x refs
  corrs <- apply(ref_counts_mat, 2, function(r) suppressWarnings(cor(r, test_counts)))
  corrs[is.na(corrs)] <- -Inf
  
  ord <- order(corrs, decreasing = TRUE)
  k <- min(k, length(ord))
  chosen <- ord[seq_len(k)]
  
  list(chosen = chosen, corrs = corrs)
}

# -----------------------------
# Discover BAMs
# -----------------------------
bam_files <- sort(list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE))
if (length(bam_files) < 1) die("No .bam files found in bam_dir: %s", bam_dir)
say("Found %d BAMs", length(bam_files))

# Require BAM indexes
idx <- vapply(bam_files, find_bam_index, character(1))
if (any(is.na(idx))) {
  missing <- bam_files[is.na(idx)]
  die("Missing BAM index for:\n%s\n\nRun: samtools index <bam>", paste(missing, collapse = "\n"))
}

# Read targets
targets <- read_targets(bed_file)
say("Read %d targets from BED", nrow(targets))
say("BED contigs (unique): %s", paste(unique(targets$chromosome), collapse = ", "))

# BAM contigs
bam_seqs <- get_bam_seqs(bam_files[1])
say("BAM contigs: %s", paste(bam_seqs, collapse = ", "))

# Verify contigs match (no chr prefix expected)
missing_in_bam <- setdiff(unique(targets$chromosome), bam_seqs)
if (length(missing_in_bam) > 0) {
  say("ERROR: These BED contigs are missing from BAM header: %s",
      paste(head(missing_in_bam, 25), collapse = ", "))
  die("Contig mismatch between BED and BAM. (Are you sure include.chr=FALSE and same reference?)")
}

# Optional chromosome filtering
if (!is.null(include_chr)) {
  keep <- targets$chromosome %in% include_chr
  say("Filtering targets by include_chr: keeping %d / %d", sum(keep), nrow(targets))
  targets <- targets[keep, , drop = FALSE]
}

# -----------------------------
# Count reads over targets (expensive; happens ONCE)
# -----------------------------
say("Counting reads over targets for all BAMs (include.chr=FALSE)...")
counts_gr <- getBamCounts(
  bed.frame   = targets,
  bam.files   = bam_files,
  include.chr = FALSE
)

# In ExomeDepth 1.1.16 this is a data.frame
if (!is.data.frame(counts_gr)) {
  say("NOTE: counts_gr is not a data.frame (class=%s). This script assumes ExomeDepth archived behavior.",
      paste(class(counts_gr), collapse = ", "))
}

say("counts_gr class: %s", paste(class(counts_gr), collapse = ", "))
say("counts_gr columns (first 12): %s", paste(names(counts_gr)[1:min(12, ncol(counts_gr))], collapse = ", "))

# Ensure expected annotation columns exist
need_cols <- c("chromosome", "start", "end")
if (!all(need_cols %in% names(counts_gr))) {
  die("counts_gr missing required columns. Has: %s", paste(names(counts_gr), collapse = ", "))
}
# Determine the name column (could be 'name' if provided)
name_col <- if ("name" %in% names(counts_gr)) "name" else NA_character_
if (is.na(name_col)) {
  # fallback: create a name column if absent
  say("WARNING: counts_gr has no 'name' column. Creating synthetic names.")
  counts_gr$name <- paste0("target_", seq_len(nrow(counts_gr)))
  name_col <- "name"
}

# Extract count matrix (all columns after the first 4 are usually samples)
# But don't assume: identify sample columns as those that are NOT annotation columns.
anno_cols <- c("chromosome", "start", "end", name_col)
sample_cols <- setdiff(names(counts_gr), anno_cols)

# Some ExomeDepth outputs include extra annotation columns; keep only numeric-ish sample cols:
is_numish <- vapply(counts_gr[, sample_cols, drop = FALSE], function(x) is.numeric(x) || is.integer(x), logical(1))
sample_cols <- sample_cols[is_numish]

if (length(sample_cols) < 2) {
  die("Could not identify sample count columns. sample_cols=%s", paste(sample_cols, collapse = ", "))
}

count_mat <- as.matrix(counts_gr[, sample_cols, drop = FALSE])
say("Counts matrix: %d targets x %d samples", nrow(count_mat), ncol(count_mat))

# Name samples nicely
colnames(count_mat) <- sub("\\.bam$", "", basename(colnames(count_mat)))
sample_names <- colnames(count_mat)

# Bin lengths (inclusive coords)
bin_len_all <- as.integer(counts_gr$end - counts_gr$start + 1L)
if (any(bin_len_all <= 0, na.rm = TRUE)) {
  say("WARNING: Found %d bins with non-positive length. Will be dropped by filters.",
      sum(bin_len_all <= 0, na.rm = TRUE))
}

# Prepare manifest
results_index <- data.frame(
  sample = sample_names,
  calls_tsv = NA_character_,
  exome_obj_rds = NA_character_,
  stringsAsFactors = FALSE
)

# -----------------------------
# Per-sample CNV calling
# -----------------------------
for (i in seq_along(sample_names)) {
  test_name <- sample_names[i]
  say("=== %s (%d/%d) ===", test_name, i, length(sample_names))
  
  test_counts_all <- as.numeric(count_mat[, i])
  
  # Reference matrix = all other samples
  ref_cols <- setdiff(seq_along(sample_names), i)
  if (length(ref_cols) < min_ref_samples) {
    say("Skipping %s: not enough other samples for reference (have %d, need >= %d)",
        test_name, length(ref_cols), min_ref_samples)
    next
  }
  ref_counts_all <- count_mat[, ref_cols, drop = FALSE]
  
  # Filtering bins
  ref_mean <- rowMeans(ref_counts_all)
  if (drop_all_zero_bins) {
    keep_idx <- which((bin_len_all > 0L) & ((test_counts_all + ref_mean) > 0))
  } else {
    keep_idx <- which(bin_len_all > 0L)
  }
  
  say("Keeping %d / %d bins after filter", length(keep_idx), length(bin_len_all))
  
  # Slice everything by keep_idx using a single annotation frame (prevents length mismatch)
  ann <- counts_gr[keep_idx, c("chromosome", "start", "end", name_col)]
  colnames(ann)[colnames(ann) == name_col] <- "name"
  
  test_counts <- as.numeric(test_counts_all[keep_idx])
  ref_counts  <- ref_counts_all[keep_idx, , drop = FALSE]
  bin_len     <- bin_len_all[keep_idx]
  
  chrom_keep <- as.character(ann$chromosome)
  start_keep <- as.integer(ann$start)
  end_keep   <- as.integer(ann$end)
  name_keep  <- as.character(ann$name)
  
  # Hard checks
  if (!(length(test_counts) == nrow(ann) &&
        length(test_counts) == length(chrom_keep) &&
        length(test_counts) == length(start_keep) &&
        length(test_counts) == length(end_keep) &&
        length(test_counts) == length(name_keep))) {
    
    say("DEBUG lengths: test=%d ann=%d chrom=%d start=%d end=%d name=%d",
        length(test_counts), nrow(ann), length(chrom_keep), length(start_keep),
        length(end_keep), length(name_keep))
    
    die("Annotation length mismatch (this should not happen with ann slicing).")
  }
  
  # Choose top-K correlated reference samples
  ref_fit <- make_reference_by_correlation(test_counts, ref_counts, k = top_k_refs)
  chosen <- ref_fit$chosen
  
  # Integer aggregate reference (IMPORTANT for beta-binomial)
  ref_agg <- as.numeric(rowSums(ref_counts[, chosen, drop = FALSE]))
  
  # Debug reference selection
  top_cor <- max(ref_fit$corrs, na.rm = TRUE)
  say("Chose %d refs (top cor = %.3f). Ref columns: %s",
      length(chosen), top_cor,
      paste(sample_names[ref_cols][chosen], collapse = ", "))
  
  # Sanity: integer-ish?
  if (any(abs(ref_agg - round(ref_agg)) > 1e-6)) {
    say("WARNING: reference has non-integer values (unexpected with rowSums).")
  }
  
  # Construct ExomeDepth object
  exome_obj <- new(
    "ExomeDepth",
    test      = test_counts,
    reference = ref_agg,
    formula   = "cbind(test, reference) ~ 1"
  )
  
  # Call CNVs
  say("Calling CNVs...")
  exome_obj <- CallCNVs(
    x = exome_obj,
    transition.probability = transition_prob,
    chromosome = chrom_keep,
    start      = start_keep,
    end        = end_keep,
    name       = name_keep
  )
  
  calls <- exome_obj@CNV.calls
  say("CNV calls for %s: %d rows", test_name, nrow(calls))
  
  # Write outputs
  calls_file <- file.path(out_dir, paste0(test_name, ".exomedepth.cnv_calls.tsv"))
  rds_file   <- file.path(out_dir, paste0(test_name, ".exomedepth.object.rds"))
  
  write.table(calls, calls_file, sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(exome_obj, rds_file)
  
  results_index$calls_tsv[i] <- calls_file
  results_index$exome_obj_rds[i] <- rds_file
  
  say("Wrote: %s", calls_file)
}

# Manifest
manifest_file <- file.path(out_dir, "exomedepth_manifest.tsv")
write.table(results_index, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)
say("Done. Manifest: %s", manifest_file)
