#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  library(tidyverse)
  library(purrr)
}))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
this_dir <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile)), error=function(e) NA)
if (is.na(this_dir) || !nzchar(this_dir)) {
  # when running via Rscript <file>, sys.frame trick may fail; fall back
  args0 <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args0, value = TRUE)
  if (length(file_arg) == 1) this_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
}
if (is.na(this_dir) || !dir.exists(this_dir)) {
  stop("Cannot resolve helper directory for SVCROWS module")
}

source(file.path(this_dir, "Scavenge.R"))
source(file.path(this_dir, "Hunt.R"))
source(file.path(this_dir, "Tools.R"))

usage <- function() {
  cat(paste0(
"SVCROWS (MAPI module helper)\n\n",
"Usage:\n",
"  svcrows <command> [args...]\n\n",
"Commands:\n",
"  scavenge   Run Scavenge mode (directory in -> directory out)\n",
"  hunt       Run Hunt mode (directory in + feature list -> directory out)\n",
"  summarize  Summarize all FCL output in an output directory\n",
"  consensus-to-query  Convert consensus output to clean query list\n",
"  concat-samples      Concatenate many per-sample files into one file\n\n",
"Scavenge positional form (legacy):\n",
"  svcrows scavenge INPUT_QUERY_LIST_DIR OUTPUT_DIR [ExpandRORegion] [BPfactor] [DefaultSizes] [xs] [xl] [y1s] [y1l] [y2s] [y2l]\n\n",
"Hunt positional form (legacy):\n",
"  svcrows hunt INPUT_QUERY_LIST_DIR FEATURELIST_TSV OUTPUT_DIR [BPfactor] [DefaultSizes] [xs] [xl] [y1s] [y1l] [y2s] [y2l]\n\n",
"Flag form (recommended):\n",
"  svcrows scavenge --in DIR --out DIR [--expand T/F] [--bpfactor T/F] [--defaultsizes T/F]\n",
"                 [--xs N --xl N --y1s N --y1l N --y2s N --y2l N]\n",
"  svcrows hunt --in DIR --features FILE --out DIR [--bpfactor T/F] [--defaultsizes T/F]\n",
"              [--xs N --xl N --y1s N --y1l N --y2s N --y2l N]\n\n",
"Notes:\n",
"  - Booleans accept: TRUE/FALSE, T/F, 1/0, yes/no\n",
"  - Numeric params are optional unless DefaultSizes=FALSE\n"
), file=stderr())
}

b <- function(x){
  if (is.null(x) || length(x)==0) return(NULL)
  s <- tolower(as.character(x))
  if (s %in% c("true","t","1","yes","y")) return(TRUE)
  if (s %in% c("false","f","0","no","n")) return(FALSE)
  if (s %in% c("na","")) return(NA)
  as.logical(x)
}

n <- function(x){
  if (is.null(x) || length(x)==0) return(NULL)
  s <- tolower(as.character(x))
  if (s %in% c("na","")) return(NA_real_)
  as.numeric(x)
}

# very small flag parser: --key value OR --key=value
parse_flags <- function(args) {
  flags <- list()
  pos <- character()

  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]

    # stop parsing flags at "--"
    if (identical(a, "--")) {
      if (i < length(args)) pos <- c(pos, args[(i+1):length(args)])
      break
    }

    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)

      # key=value form: --foo=bar
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1]]
        k <- parts[[1]]
        v <- paste(parts[-1], collapse = "=")
        flags[[k]] <- v
      } else {
        # if next token exists and is not another flag, treat as value
        if (i < length(args) && !startsWith(args[[i+1]], "--")) {
          flags[[key]] <- args[[i+1]]
          i <- i + 1
        } else {
          # boolean flag: --foo (no value) => TRUE
          flags[[key]] <- "TRUE"
        }
      }
    } else {
      pos <- c(pos, a)
    }

    i <- i + 1
  }

  list(flags = flags, pos = pos)
}

argv <- commandArgs(trailingOnly=TRUE)
if (length(argv)==0 || argv[[1]] %in% c("-h","--help","help")) {
  usage(); quit(status=0)
}

cmd <- argv[[1]]
args <- argv[-1]

# ------------------------------------------------------------
# Dispatch
# ------------------------------------------------------------
if (cmd == "scavenge") {
  if (length(args)==0 || args[[1]] %in% c("-h","--help")) { usage(); quit(status=0) }

  if (length(args) >= 2 && !startsWith(args[[1]], "--")) {
    # Legacy positional
    params <- list(InputQueryList = args[[1]], OutputDirectory = args[[2]])
    if (length(args) >= 3)  params$ExpandRORegion <- b(args[[3]])
    if (length(args) >= 4)  params$BPfactor       <- b(args[[4]])
    if (length(args) >= 5)  params$DefaultSizes   <- b(args[[5]])
    if (length(args) >= 6)  params$xs             <- n(args[[6]])
    if (length(args) >= 7)  params$xl             <- n(args[[7]])
    if (length(args) >= 8)  params$y1s            <- n(args[[8]])
    if (length(args) >= 9)  params$y1l            <- n(args[[9]])
    if (length(args) >= 10) params$y2s            <- n(args[[10]])
    if (length(args) >= 11) params$y2l            <- n(args[[11]])
  } else {
    p <- parse_flags(args)
    f <- p$flags
    if (is.null(f[['in']]) || is.null(f[['out']])) stop("scavenge requires --in and --out")
    params <- list(
      InputQueryList = f[['in']],
      OutputDirectory = f[['out']]
    )
    if (!is.null(f[['expand']]))       params$ExpandRORegion <- b(f[['expand']])
    if (!is.null(f[['bpfactor']]))     params$BPfactor       <- b(f[['bpfactor']])
    if (!is.null(f[['defaultsizes']])) params$DefaultSizes   <- b(f[['defaultsizes']])
    if (!is.null(f[['xs']]))           params$xs             <- n(f[['xs']])
    if (!is.null(f[['xl']]))           params$xl             <- n(f[['xl']])
    if (!is.null(f[['y1s']]))          params$y1s            <- n(f[['y1s']])
    if (!is.null(f[['y1l']]))          params$y1l            <- n(f[['y1l']])
    if (!is.null(f[['y2s']]))          params$y2s            <- n(f[['y2s']])
    if (!is.null(f[['y2l']]))          params$y2l            <- n(f[['y2l']])
  } # <-- IMPORTANT: closes the positional-vs-flag if/else
  # (and also closes the outer else properly now)

  odir <- params$OutputDirectory
  dir.create(odir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(odir)) stop("Cannot create OutputDirectory: ", odir)

  cat("[svcrows] Calling Scavenge...\n", file=stderr())
  Scavenge(
    InputQueryList = params$InputQueryList,
    OutputDirectory = params$OutputDirectory,
    ExpandRORegion = if (!is.null(params$ExpandRORegion)) params$ExpandRORegion else FALSE,
    BPfactor       = if (!is.null(params$BPfactor)) params$BPfactor else TRUE,
    DefaultSizes   = if (!is.null(params$DefaultSizes)) params$DefaultSizes else FALSE,
    xs  = if (!is.null(params$xs))  params$xs  else 5000,
    xl  = if (!is.null(params$xl))  params$xl  else 25000,
    y1s = if (!is.null(params$y1s)) params$y1s else 500,
    y1l = if (!is.null(params$y1l)) params$y1l else 2500,
    y2s = if (!is.null(params$y2s)) params$y2s else 50,
    y2l = if (!is.null(params$y2l)) params$y2l else 80
  )
  cat("[svcrows] Done\n", file=stderr())

} else if (cmd == "hunt") {
  if (length(args)==0 || args[[1]] %in% c("-h","--help")) { usage(); quit(status=0) }

  if (length(args) >= 3 && !startsWith(args[[1]], "--")) {
    # Legacy positional
    params <- list(InputQueryList=args[[1]], FeatureList=args[[2]], OutputDirectory=args[[3]])
    if (length(args) >= 4)  params$BPfactor     <- b(args[[4]])
    if (length(args) >= 5)  params$DefaultSizes <- b(args[[5]])
    if (length(args) >= 6)  params$xs           <- n(args[[6]])
    if (length(args) >= 7)  params$xl           <- n(args[[7]])
    if (length(args) >= 8)  params$y1s          <- n(args[[8]])
    if (length(args) >= 9)  params$y1l          <- n(args[[9]])
    if (length(args) >= 10) params$y2s          <- n(args[[10]])
    if (length(args) >= 11) params$y2l          <- n(args[[11]])
  } else {
    p <- parse_flags(args)
    f <- p$flags
    if (is.null(f[['in']]) || is.null(f[['features']]) || is.null(f[['out']])) {
      stop("hunt requires --in, --features, and --out")
    }
    params <- list(InputQueryList=f[['in']], FeatureList=f[['features']], OutputDirectory=f[['out']])
    if (!is.null(f[['bpfactor']]))     params$BPfactor     <- b(f[['bpfactor']])
    if (!is.null(f[['defaultsizes']])) params$DefaultSizes <- b(f[['defaultsizes']])
    if (!is.null(f[['xs']]))           params$xs           <- n(f[['xs']])
    if (!is.null(f[['xl']]))           params$xl           <- n(f[['xl']])
    if (!is.null(f[['y1s']]))          params$y1s          <- n(f[['y1s']])
    if (!is.null(f[['y1l']]))          params$y1l          <- n(f[['y1l']])
    if (!is.null(f[['y2s']]))          params$y2s          <- n(f[['y2s']])
    if (!is.null(f[['y2l']]))          params$y2l          <- n(f[['y2l']])
  }

  odir <- params$OutputDirectory
  dir.create(odir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(odir)) stop("Cannot create OutputDirectory: ", odir)

  cat("[svcrows] Calling Hunt...\n", file=stderr())
  Hunt(
    InputQueryList = params$InputQueryList,
    FeatureList = params$FeatureList,
    OutputDirectory = params$OutputDirectory,
    BPfactor     = if (!is.null(params$BPfactor)) params$BPfactor else TRUE,
    DefaultSizes = if (!is.null(params$DefaultSizes)) params$DefaultSizes else FALSE,
    xs  = if (!is.null(params$xs))  params$xs  else 3000,
    xl  = if (!is.null(params$xl))  params$xl  else 6000,
    y1s = if (!is.null(params$y1s)) params$y1s else 300,
    y1l = if (!is.null(params$y1l)) params$y1l else 600,
    y2s = if (!is.null(params$y2s)) params$y2s else 30,
    y2l = if (!is.null(params$y2l)) params$y2l else 60
  )
  cat("[svcrows] Done\n", file=stderr())

} else if (cmd == "summarize") {
  if (length(args) < 2 || args[[1]] %in% c("-h","--help")) {
    cat("Usage: svcrows summarize OUTPUT_DIR OUTPUT_FILE\n", file=stderr())
    quit(status=0)
  }
  SVCSummarize(OutputDirectory = args[[1]], OutputFile = args[[2]])

} else if (cmd == "consensus-to-query") {
  if (length(args) < 2 || args[[1]] %in% c("-h","--help")) {
    cat("Usage: svcrows consensus-to-query CONSENSUS_TSV QUERY_OUT_TSV\n", file=stderr())
    quit(status=0)
  }
  ConsensusToQuery(args[[1]], args[[2]])

} else if (cmd == "concat-samples") {
  if (length(args) < 2 || args[[1]] %in% c("-h","--help")) {
    cat("Usage: svcrows concat-samples IN_DIR OUT_DIR\n", file=stderr())
    quit(status=0)
  }
  ConcatSamples(args[[1]], args[[2]])

} else {
  cat("Unknown command: ", cmd, "\n\n", sep="", file=stderr())
  usage()
  quit(status=2)
}
