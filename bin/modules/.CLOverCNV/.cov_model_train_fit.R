#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  # NOTE: No optparse. We parse args manually for portability.
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Missing package: mgcv")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing package: data.table")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing package: jsonlite")
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing value after", flag))
  args[[i + 1]]
}

# Allow repeated flags: --exclude-contigs A --exclude-contigs B
get_args_multi <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(character())
  vals <- character()
  for (i in idx) {
    if (i == length(args)) stop(paste("Missing value after", flag))
    vals <- c(vals, args[[i + 1]])
  }
  vals
}

# Parse exclude contigs string(s) into vector
parse_contigs <- function(x) {
  if (length(x) == 0) return(character())
  # x may be multiple strings; split on commas/whitespace
  parts <- unlist(strsplit(paste(x, collapse = " "), "[,\\s]+"))
  parts <- parts[nzchar(parts)]
  unique(parts)
}

depth_path <- get_arg("--depth")
nuc_path   <- get_arg("--nuc")
map_path   <- get_arg("--map")
window     <- as.integer(get_arg("--window", "100"))
seed       <- as.integer(get_arg("--seed", "42"))

# Optional filters
min_window_bp <- as.integer(get_arg("--min-window-bp", "0"))  # 0 = off
exclude_contigs <- parse_contigs(get_args_multi("--exclude-contigs"))

out_model  <- get_arg("--out-model")
out_meta   <- get_arg("--out-meta")
out_table  <- get_arg("--out-table")
out_plot   <- get_arg("--out-plot")

stopifnot(!is.null(depth_path), !is.null(nuc_path), !is.null(map_path),
          !is.null(out_model), !is.null(out_meta), !is.null(out_table), !is.null(out_plot))

library(data.table)
library(mgcv)
library(jsonlite)

set.seed(seed)

# -----------------------------
# Read inputs
# -----------------------------
depth <- fread(depth_path, header = FALSE)
setnames(depth, c("chrom", "start", "end", "mean_depth"))

# bedtools nuc: header present, many columns
nuc <- fread(nuc_path, header = TRUE)

# Prefer pct_at if available; else fall back to pct_gc
pct_at_col <- NULL
for (nm in c("pct_at", "pct.AT", "pct_at_1", "X.pct_at")) {
  if (nm %in% names(nuc)) { pct_at_col <- nm; break }
}
if (is.null(pct_at_col)) {
  hits <- grep("pct.*at", names(nuc), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) pct_at_col <- hits[[1]]
}

pct_gc_col <- NULL
for (nm in c("pct_gc", "pct.GC", "pct_gc_1", "X.pct_gc")) {
  if (nm %in% names(nuc)) { pct_gc_col <- nm; break }
}
if (is.null(pct_gc_col)) {
  hits <- grep("pct.*gc", names(nuc), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) pct_gc_col <- hits[[1]]
}

if (is.null(pct_at_col) && is.null(pct_gc_col)) {
  stop("Could not find pct_at or pct_gc column in bedtools nuc output. Columns: ",
       paste(names(nuc), collapse = ", "))
}

# Keep only essential columns; bedtools nuc chrom/start/end are usually first three
# Keep pct_at if present, else pct_gc
if (!is.null(pct_at_col)) {
  nuc_sub <- nuc[, .(
    chrom = get(names(nuc)[1]),
    start = as.integer(get(names(nuc)[2])),
    end   = as.integer(get(names(nuc)[3])),
    pct_at_raw = as.numeric(get(pct_at_col))
  )]
} else {
  nuc_sub <- nuc[, .(
    chrom = get(names(nuc)[1]),
    start = as.integer(get(names(nuc)[2])),
    end   = as.integer(get(names(nuc)[3])),
    pct_gc_raw = as.numeric(get(pct_gc_col))
  )]
}

# Map file: chrom start end map_mean (NA allowed)
map <- fread(map_path, header = FALSE)
setnames(map, c("chrom", "start", "end", "map_mean"))
map[, map_mean := suppressWarnings(as.numeric(map_mean))]

# -----------------------------
# Join
# -----------------------------
setkey(depth, chrom, start, end)
setkey(nuc_sub, chrom, start, end)
setkey(map, chrom, start, end)

dt <- depth[nuc_sub][map]

# -----------------------------
# Optional filters
# -----------------------------
# Exclude contigs
if (length(exclude_contigs) > 0) {
  dt <- dt[!(chrom %in% exclude_contigs)]
}

# Size filter (off by default)
# Useful if windows are variable-sized; for fixed makewindows it's mostly redundant.
dt[, win_bp := as.integer(end - start)]
if (!is.na(min_window_bp) && min_window_bp > 0) {
  dt <- dt[win_bp >= min_window_bp]
}

# -----------------------------
# Features: pct_at scaling
# -----------------------------
# bedtools nuc may output pct_at/pct_gc either as:
#   - fractions in [0,1]  (as you observed)
#   - percents in [0,100]
#
# We normalize to AT fraction in [0,1].
if ("pct_at_raw" %in% names(dt)) {
  x <- dt$pct_at_raw
  if (all(is.na(x))) stop("pct_at column is all NA after join.")
  # Detect scaling: if max > 1.5 assume 0-100 percent
  if (max(x, na.rm = TRUE) > 1.5) x <- x / 100
  dt[, pct_at := pmax(0, pmin(1, x))]
} else {
  x <- dt$pct_gc_raw
  if (all(is.na(x))) stop("pct_gc column is all NA after join.")
  if (max(x, na.rm = TRUE) > 1.5) x <- x / 100
  dt[, pct_at := pmax(0, pmin(1, 1 - x))]
}

dt[, map_mean := ifelse(is.na(map_mean), NA_real_, map_mean)]

# Response: mean_depth; stabilize with log1p
dt <- dt[!is.na(mean_depth)]
dt[, y := log1p(mean_depth)]

# -----------------------------
# Fit model (GAM)
# -----------------------------
use_map <- any(!is.na(dt$map_mean))

if (use_map) {
  fit <- gam(y ~ s(pct_at, k = 20) + s(map_mean, k = 20), data = dt, method = "REML")
} else {
  fit <- gam(y ~ s(pct_at, k = 20), data = dt, method = "REML")
}

# -----------------------------
# Compute adjustment factor
# -----------------------------
pred_y <- predict(fit, newdata = dt)
pred_depth <- pmax(1e-9, expm1(pred_y))
baseline <- median(pred_depth, na.rm = TRUE)
dt[, pred_depth := pred_depth]
dt[, adj := baseline / pred_depth]

# -----------------------------
# Write outputs
# -----------------------------
out_dt <- dt[, .(chrom, start, end, mean_depth, pct_at, map_mean, pred_depth, adj)]
fwrite(out_dt, file = out_table, sep = "\t", compress = "gzip")

saveRDS(fit, file = out_model)

meta <- list(
  window_bp = window,
  seed = seed,
  used_mappability = use_map,
  exclude_contigs = exclude_contigs,
  min_window_bp = min_window_bp,
  nuc_pct_scale = "auto-detected (fraction 0-1 vs percent 0-100)",
  response = "log1p(mean_depth)",
  adjustment = "adj = median(pred_depth) / pred_depth; corrected = raw * adj"
)
writeLines(toJSON(meta, pretty = TRUE, auto_unbox = TRUE), con = out_meta)

# Diagnostics plot
png(out_plot, width = 1400, height = 900)
par(mfrow = c(1, 2))
plot(dt$pct_at, dt$mean_depth, pch = 16, cex = 0.3,
     xlab = "AT fraction (per window)", ylab = "Mean depth (raw)",
     main = "Raw depth vs AT")
plot(dt$pct_at, dt$adj, pch = 16, cex = 0.3,
     xlab = "AT fraction (per window)", ylab = "Adjustment factor",
     main = "Adjustment vs AT")
dev.off()

cat("OK\n")
