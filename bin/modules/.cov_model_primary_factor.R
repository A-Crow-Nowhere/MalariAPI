#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(NA_character_)
  args[[i + 1]]
}

features_path <- get_arg("--features")
model_path    <- get_arg("--model")
out_expected  <- get_arg("--out-expected")
out_cfactor   <- get_arg("--out-cfactor")

if (is.na(features_path) || is.na(model_path) || is.na(out_expected) || is.na(out_cfactor)) {
  stop("Usage: Rscript .cov_model_primary_factor.R --features X --model Y --out-expected A --out-cfactor B")
}

saye <- function(...) cat("[cov_model_primary_factor.R] ", ..., "\n", sep="", file=stderr())

dt <- read.delim(features_path, header=TRUE, sep="\t",
                 stringsAsFactors=FALSE, check.names=FALSE)

to_num <- function(x) suppressWarnings(as.numeric(x))
for (nm in intersect(c("start","end","depth","AT","M","pct_at"), names(dt))) {
  dt[[nm]] <- to_num(dt[[nm]])
}

model <- readRDS(model_path)

saye("model class: ", paste(class(model), collapse=","))
saye("model formula: ", paste(deparse(formula(model)), collapse=" "))

tt <- tryCatch(terms(model), error=function(e) NULL)
if (!is.null(tt)) {
  needed_terms <- all.vars(delete.response(tt))
  saye("needed vars (terms): ", paste(needed_terms, collapse=", "))
}

# Ensure pct_at exists
if (!("pct_at" %in% names(dt))) {
  if ("AT" %in% names(dt)) {
    x <- dt[["AT"]]
    # AT usually 0..1 => convert to percent
    if (all(is.na(x) | (x >= 0 & x <= 1.2))) {
      dt[["pct_at"]] <- x * 100.0
    } else {
      dt[["pct_at"]] <- x
    }
    saye("created pct_at from AT")
  } else {
    stop("Neither pct_at nor AT found in features table.")
  }
}

saye("pct_at range (as provided): ", paste(range(dt$pct_at, na.rm=TRUE), collapse=" .. "))

# Validate model.frame on predictors only (don't require response)
if (!is.null(tt)) {
  ttx <- delete.response(tt)
  mf_try <- tryCatch(model.frame(ttx, data=dt, na.action=na.pass), error=function(e) e)
  if (inherits(mf_try, "error")) {
    saye("ERROR building model.frame(): ", conditionMessage(mf_try))
    stop("model.frame() failed; predictors in features do not match what the model expects.")
  }
}

# Pull training-scale hints from mgcv
train_rng <- NULL
if (!is.null(model$var.summary) && "pct_at" %in% names(model$var.summary)) {
  train_rng <- range(model$var.summary[["pct_at"]], na.rm=TRUE)
} else if (!is.null(model$model) && "pct_at" %in% names(model$model)) {
  train_rng <- range(model$model[["pct_at"]], na.rm=TRUE)
}
if (!is.null(train_rng)) {
  saye("pct_at training range (from model): ", paste(train_rng, collapse=" .. "))
} else {
  saye("pct_at training range: <unknown> (model lacks var.summary/model$pct_at)")
}

predict_expected <- function(newdata) {
  pr <- as.numeric(mgcv::predict.gam(model, newdata=newdata, type="response"))
  pr[!is.finite(pr)] <- NA_real_
  pr[pr <= 0] <- NA_real_
  pr
}

# Try multiple scalings; choose the first that works
candidates <- list(
  list(label="as_is",  fn=function(z) z),
  list(label="div100", fn=function(z) z/100.0),
  list(label="mul100", fn=function(z) z*100.0)
)

# If we know training range, reorder tries to favor the likely correct scaling
if (!is.null(train_rng)) {
  # If training looks like fraction-scale (<= ~1.2), prefer /100 first
  if (train_rng[2] <= 1.2) {
    candidates <- list(
      list(label="div100", fn=function(z) z/100.0),
      list(label="as_is",  fn=function(z) z),
      list(label="mul100", fn=function(z) z*100.0)
    )
  }
  # If training looks like percent-scale (>= ~10), prefer as-is first (already)
}

pred <- rep(NA_real_, nrow(dt))
used <- NA_character_

for (cand in candidates) {
  dt_try <- dt
  dt_try$pct_at <- cand$fn(dt$pct_at)

  rng_try <- range(dt_try$pct_at, na.rm=TRUE)
  saye("trying pct_at scaling '", cand$label, "' range=", paste(rng_try, collapse=" .. "))

  pr <- tryCatch(predict_expected(dt_try), error=function(e) e)
  if (inherits(pr, "error")) {
    saye("predict ERROR under '", cand$label, "': ", conditionMessage(pr))
    next
  }
  if (!all(is.na(pr))) {
    pred <- pr
    used <- cand$label
    break
  }
}

if (is.na(used) || all(is.na(pred))) {
  stop("predict() returned all NA/nonpositive for expected depth under all tested pct_at scalings. Model object may be unusable for prediction.")
}

saye("prediction succeeded using pct_at scaling: ", used)

# Correction factor
scale <- median(pred, na.rm=TRUE)
cf <- scale / pred
cf[!is.finite(cf)] <- NA_real_

# Clamp extremes (tweak/remove later)
cf[cf < 0.2] <- 0.2
cf[cf > 5.0] <- 5.0

# Enforce integer coords + valid intervals
dt[["start"]] <- as.integer(round(dt[["start"]]))
dt[["end"]]   <- as.integer(round(dt[["end"]]))
ok <- !is.na(dt[["start"]]) & !is.na(dt[["end"]]) & dt[["start"]] >= 0L & dt[["end"]] > dt[["start"]]

dt <- dt[ok, , drop=FALSE]
pred <- pred[ok]
cf   <- cf[ok]

outE <- data.frame(chrom=dt[["chrom"]], start=dt[["start"]], end=dt[["end"]], val=pred)
outC <- data.frame(chrom=dt[["chrom"]], start=dt[["start"]], end=dt[["end"]], val=cf)

write.table(outE, file=out_expected, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(outC, file=out_cfactor,  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
