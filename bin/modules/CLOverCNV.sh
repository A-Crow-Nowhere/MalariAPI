#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# CLOverCNV.sh
#   Wrapper pipeline:
#     Model_Train -> Model_Apply -> Model_Finalize (bw QC)
#     CNV_Probe -> CNV_Segment -> CNV_Boundary -> CNV_Confidence -> CNV_Finalize
#
#   Also supports:
#     mapi modules CLOverCNV step <AliasOrInternal> [options...]
#     mapi modules CLOverCNV run  [pipeline options...]
#
# R handling notes:
#   - Conda R is supported (YAML includes r-changepoint, etc.)
#   - For HPC module R, this wrapper:
#       * auto-selects a writable R_LIBS_USER (HOME or /standard)
#       * installs missing CRAN pkgs (data.table, R.utils, changepoint) if needed
#       * does NOT touch .R scripts
# ============================================================

# -----------------------------
# Logging / errors
# -----------------------------
LOG_TS_FMT='+%Y-%m-%d %H:%M:%S'
msg() { echo "[CLOverCNV] $(date "$LOG_TS_FMT") $*" >&2; }
die() { msg "ERROR: $*"; exit 1; }

# -----------------------------
# Determine MAPI_ROOT (robust)
# -----------------------------
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAPI_ROOT="${MAPI_ROOT:-"$(cd "$THIS_DIR/../.." && pwd)"}"
export MAPI_ROOT

msg "DEBUG: wrapper path = ${BASH_SOURCE[0]}"
msg "DEBUG: MAPI_ROOT = $MAPI_ROOT"

# Hidden tool dir (modules live here, wrapper stays visible)
HIDDEN_DIR="$MAPI_ROOT/bin/modules/.CLOverCNV"
[[ -d "$HIDDEN_DIR" ]] || die "Missing hidden dir: $HIDDEN_DIR (create + move helpers there)"

# Conda (used for PY steps; NOT for R in CNV_Segment by default)
CONDA="$MAPI_ROOT/tools/miniconda3/bin/conda"
[[ -x "$CONDA" ]] || die "conda not found/executable at: $MAPI_ROOT/tools/miniconda3/bin/conda"

# CLOverCNV env path (python deps etc.)
CNV_ENV="$MAPI_ROOT/envs/CLOverCNV"
[[ -d "$CNV_ENV" ]] || msg "WARN: env dir not found (will still try PATH tools): $CNV_ENV"

# -----------------------------
# Defaults (run)
# -----------------------------
SAMPLE=""
GENOME=""
GFF=""
OUT_ROOT=""
BAM=""
BAM_DIR=""
MODEL_IN=""

# Model train/apply defaults (match your modules)
TRAIN_WINDOW=100
TRAIN_MAPQ=20
TRAIN_THREADS=8
TRAIN_DROP_DUP=1
TRAIN_EXCLUDE_CONTIGS=()

APPLY_BINSIZE=50
APPLY_THREADS=8
APPLY_EMIT_RAW=0
APPLY_PRIMARY_PREFER="primary"  # primary|aligned

# R runtime selection for HPC R steps
R_MODULES="off"   # off|none|auto|"<mods...>"
R_AUTOINSTALL="${R_AUTOINSTALL:-1}"
MAPI_R_LIBS_USER_ROOT="${MAPI_R_LIBS_USER_ROOT:-}"

# CNV defaults
TILE_BP=75
MIN_TILE_BP=50
EXCLUDE_BED=""
LAMBDA=25
MIN_PROBES_PER_SEG=10
FLANK=1250
AGG="sum"
CONF_METHOD="min"

FINAL_WEAK_RATIO=2.0
FINAL_STRONG_RATIO=2.3
FINAL_WEAK_Z=1.5
FINAL_STRONG_Z=2.0
FINAL_KEEP_MODE="both"

# NEW: CNV finalize fusion options (post-call consolidation)
FINAL_FUSE=0
FINAL_FUSE_MAX_GAP=0

# NEW: CNV finalize BAM-based counting options
FINAL_COUNT_FLANK=1250
FINAL_COUNT_ALL_ALIGN=0
FINAL_COUNT_DROP_DUP=0
FINAL_COUNT_INCLUDE_SUPP=0
FINAL_COUNT_INCLUDE_SECONDARY=0

FORCE=0
KEEP_TMP=0

START_AT="cov_model_train"
STOP_AFTER="cnv_finalize_segments"

# -----------------------------
# Step registry + aliases
# -----------------------------
declare -A STEP_ALIAS=(
  ["Model_Train"]="cov_model_train"
  ["Model_Apply"]="cov_model_apply"
  ["Model_Finalize"]="cov_bw_qc"
  ["CNV_Probe"]="make_probes"
  ["CNV_Segment"]="fused_lasso_segment"
  ["CNV_Boundary"]="boundary_support_from_bigwigs"
  ["CNV_Confidence"]="segment_confidence_from_boundaries"
  ["CNV_Finalize"]="cnv_finalize_segments"
)

steps=(
  cov_model_train
  cov_model_apply
  cov_bw_qc
  make_probes
  fused_lasso_segment
  boundary_support_from_bigwigs
  segment_confidence_from_boundaries
  cnv_finalize_segments
)

resolve_step() {
  local s="$1"
  if [[ -n "${STEP_ALIAS[$s]:-}" ]]; then
    echo "${STEP_ALIAS[$s]}"
    return 0
  fi
  for x in "${steps[@]}"; do
    [[ "$x" == "$s" ]] && { echo "$s"; return 0; }
  done
  echo ""
}

idx() {
  local needle="$1"
  for i in "${!steps[@]}"; do
    [[ "${steps[$i]}" == "$needle" ]] && { echo "$i"; return 0; }
  done
  echo "-1"
}

# -----------------------------
# Usage
# -----------------------------
usage() {
  cat <<'EOF'
CLOverCNV

Subcommands:
  run       Run the full pipeline (or a slice) for one sample
  step      Run a single step (alias or internal name) for one sample dir

Step aliases (recommended):
  Model_Train      Train coverage model from a BAM
  Model_Apply      Apply model to a BAM directory -> corr_bw/ raw_bw/ factors/
  Model_Finalize   QC + optional dedup bigWigs summary (cov_bw_qc)
  CNV_Probe        Make CDS-tiled probes from corrected primary bigWig
  CNV_Segment      Segment probes with PELT (R changepoint; penalty = --lambda)
  CNV_Boundary     Score boundaries with split/disco bigWigs
  CNV_Confidence   Convert boundary support into per-segment confidence
  CNV_Finalize     Call CNV states + emit bedGraphs

Glossary (what these inputs/outputs mean)

  BAM (--bam):
    The primary BAM for the sample (typically your main “primary” or “aligned” BAM).
    Used for model training. Must be indexed (.bai).

  BAM directory (--bam-dir):
    A directory containing *all* BAMs for the same sample that Model_Apply may use
    (e.g. primary/aligned + splitters + discordant, and optionally supplementary/unmapped).
    This should be one sample’s BAM set, not a mixed multi-sample folder. Must include .bai indices.

  Bin size (--apply-binsize):
    The window size (bp) used to summarize coverage into bigWig tracks during model application.
    Smaller bins = higher resolution but noisier; larger bins = smoother but less sensitive to small events.

  Probes (--tile-bp / CNV_Probe):
    Fixed-size tiles laid across CDS (exome-focused). Each probe gets a normalized coverage value.
    These probe values are the *only* signal used for CNV segmentation (not full-genome depth).

  Segments (CNV_Segment output):
    Contiguous runs of probes with a similar mean log2 ratio (y). Segments can span large genomic ranges
    because probes exist only where CDS exists; the segment “span” is a visualization envelope.

  Boundaries (CNV_Boundary / CNV_Confidence):
    The change-points between adjacent segments. Boundaries do not have copy-number themselves;
    they are events scored using split/disco evidence in flanking windows (e.g. ±--flank bp).

  Ratio / y (CNV_Finalize):
    y is log2(ratio). ratio > 1 indicates gains; ratio < 1 indicates losses.
    Effect tiers come from ratio thresholds; support tiers come from boundary z-scores.

  Optional fusion (--final-fuse):
    If enabled, adjacent *kept* CNVs of the same direction (gain/loss) are merged when they are close
    (gap ≤ --final-fuse-max-gap). Fused CN is probe-weighted; confidence uses only the outer boundary z-scores.
    Split/disco biology from internal fused-away boundaries is preserved by summarizing internal supports.

Run (pipeline) minimal inputs:
  mapi modules CLOverCNV run \
    --sample <SAMPLE> \
    --genome <GENOME_KEY|/path/to.fa[.gz]> \
    --gff <genes.gff> \
    --bam <primary.bam> \
    --bam-dir <dir_with_all_bams_for_sample> \
    --out-root <OUTDIR>

Common run options:
  --start-at <Alias|internal>     [Model_Train]
  --stop-after <Alias|internal>   [CNV_Finalize]
  --force                         overwrite step outputs if present
  --keep-tmp                      keep _work temp dirs

Normalization knobs:
  --train-window <int>            [100] model training window bp
  --train-mapq <int>              [20]
  --train-threads <int>           [8]
  --train-drop-dup                drop duplicates for training depth only (flag 1024)
  --train-exclude-contigs <list>  repeatable; e.g. "PfDd2_MT,PfDd2_API"
  --r-modules <mode|string>       R runtime for *HPC module* steps (NOT conda R)
                                  off|none     : do not module load, use PATH Rscript
                                  auto|R       : auto-select latest R/<ver> + prereqs
                                  "<mods...>"  : explicit modules, e.g. "gcc/... openmpi/... R/4.4.1"

  --apply-binsize <int>           [50]
  --apply-threads <int>           [8]
  --apply-emit-raw                emit raw bigWigs
  --apply-primary-prefer <p|a>    [primary] or aligned

CNV knobs:
  --tile-bp <int>                 [75]
  --min-tile-bp <int>             [50]
  --exclude-bed <bed>             optional exclude for probes
  --lambda <float>                [25] PELT manual penalty
  --min-probes-per-seg <int>      [10]
  --flank <int>                   [1250] boundary support flank bp
  --agg <sum|mean>                [sum]
  --confidence-method <min|mean|max> [min]
  --final-weak-ratio <float>      [2.0]
  --final-strong-ratio <float>    [2.3]
  --final-weak-z <float>          [1.5]
  --final-strong-z <float>        [2.0]
  --final-keep-mode <ratio_only|z_only|both> [both]

Finalize extras:
  --final-fuse                    enable fusion of adjacent kept CNVs (same direction)
  --final-fuse-max-gap <int>      [0] max bp gap allowed between fused segments

  --final-count-flank <int>       [1250] flank bp for BAM-based split/disco boundary read counting
  --final-count-all-alignments    count all alignments (default counts read1 only for PE)
  --final-count-drop-dup          exclude duplicate-marked reads from BAM counts
  --final-count-include-supp      include supplementary alignments in BAM counts
  --final-count-include-secondary include secondary alignments in BAM counts

Early stop convenience:
  --stop-after Model_Finalize  will stop after cov_bw_qc (no CNV calling)

Step (single-step) usage:
  mapi modules CLOverCNV step <Alias|internal> --sample <S> --out-root <OUTROOT> [other inputs...]

EOF
}

# -----------------------------
# R user-library helpers (HPC-safe)
# -----------------------------
pick_r_user_lib_root() {
  if [[ -n "${MAPI_R_LIBS_USER_ROOT:-}" ]]; then
    echo "$MAPI_R_LIBS_USER_ROOT"
    return 0
  fi

  if [[ -d "/standard" && -n "${HPC_ALLOCATION:-}" ]]; then
    local cand="/standard/${HPC_ALLOCATION}/${USER}/R/goolf"
    if mkdir -p "$cand" >/dev/null 2>&1; then
      echo "$cand"
      return 0
    fi
  fi

  echo "$HOME/R/goolf"
}

ensure_r_userlib_and_pkgs() {
  command -v Rscript >/dev/null 2>&1 || die "Rscript not found (cannot configure R user library)"

  local r_mm root
  r_mm="$(Rscript --vanilla -e 'cat(paste0(R.version$major,".",strsplit(R.version$minor,"\\.")[[1]][1]))')"
  root="$(pick_r_user_lib_root)"

  export R_LIBS_USER="${root}/${r_mm}"
  mkdir -p "$R_LIBS_USER" || die "Failed to create R_LIBS_USER: $R_LIBS_USER"
  msg "DEBUG: R_LIBS_USER=$R_LIBS_USER"

  if [[ "${R_AUTOINSTALL:-1}" -eq 1 ]]; then
    if ! Rscript --vanilla -e '
req <- c("data.table","R.utils","changepoint")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  message("[R preflight] installing: ", paste(miss, collapse=", "))
  install.packages(miss, repos="https://cloud.r-project.org")
} else {
  message("[R preflight] ok")
}
' >/dev/null 2>&1; then
      return 1
    fi
  fi

  Rscript --vanilla -e '
req <- c("data.table","R.utils","changepoint")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  cat("MISSING:", paste(miss, collapse=" "), "\n")
  print(.libPaths())
  quit(status=2)
}
' >/dev/null 2>&1 || return 1

  if [[ "${MAPI_DEBUG:-0}" == "1" ]]; then
    Rscript --vanilla -e 'cat("DEBUG .libPaths():\n"); print(.libPaths())' >&2 || true
  fi
}

# -----------------------------
# Module system helpers (for R)
# -----------------------------
maybe_enable_modules() {
  if command -v module >/dev/null 2>&1; then return 0; fi
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh && command -v module >/dev/null 2>&1 && return 0
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash && command -v module >/dev/null 2>&1 && return 0
  return 1
}

pick_latest_r_module() {
  module spider R 2>&1 \
    | awk '/Versions:/{flag=1;next} flag && $1 ~ /^R\/[0-9]/{print $1}' \
    | sort -V \
    | tail -n 1
}

extract_prereq_combos() {
  local rmod="$1"
  module spider "$rmod" 2>&1 | awk '
    function trim(s){ sub(/^[ \t]+/,"",s); sub(/[ \t]+$/,"",s); return s }
    /You will need to load all module\(s\) on any one of the lines below/ {flag=1; next}
    flag && NF==0 {exit}
    flag {
      line=$0
      if (match(line, /^[ \t]+/) && prev != "") {
        prev = prev " " trim(line)
      } else {
        if (prev != "") print trim(prev)
        prev = trim(line)
      }
      next
    }
    END { if (prev != "") print trim(prev) }
  ' | awk '{$1=$1;print}' | sed '/^$/d'
}

try_module_load_silent() {
  local spec="$1"
  # shellcheck disable=SC2086
  module load $spec >/dev/null 2>&1
}

load_r_with_prereqs() {
  local rmod="$1"
  if try_module_load_silent "$rmod"; then return 0; fi
  mapfile -t combos < <(extract_prereq_combos "$rmod" || true)
  [[ "${#combos[@]}" -gt 0 ]] || return 1
  for combo in "${combos[@]}"; do
    msg "Trying prerequisites for $rmod: $combo"
    module purge >/dev/null 2>&1 || true
    if ! try_module_load_silent "$combo"; then
      continue
    fi
    if try_module_load_silent "$rmod"; then
      return 0
    fi
  done
  return 1
}

run_r_file_hpc() {
  local script="$1"; shift
  [[ -f "$script" ]] || die "R script not found: $script"

  export TMPDIR="$WORK/tmp_R"
  mkdir -p "$TMPDIR"

  local rm
  rm="$(echo "${R_MODULES:-off}" | awk '{$1=$1;print}')"
  case "${rm,,}" in
    ""|"off"|"none") rm="" ;;
    "auto"|"r") rm="AUTO" ;;
    "conda") rm="" ;;  # treat as "no module loads" (still uses PATH Rscript)
    *) : ;;
  esac

  if [[ -n "$rm" ]]; then
    maybe_enable_modules || die "R modules requested but module system not available on this node"
    module purge >/dev/null 2>&1 || true

    if [[ "$rm" == "AUTO" ]]; then
      local rpick
      rpick="$(pick_latest_r_module || true)"
      [[ -n "$rpick" ]] || die "Could not find R/<version> via 'module spider R'"
      msg "Auto-selected R module: $rpick"
      load_r_with_prereqs "$rpick" || die "Failed loading $rpick (auto). Use explicit --r-modules \"gcc/... openmpi/... $rpick\""
    else
      msg "Loading HPC modules for R: $rm"
      # shellcheck disable=SC2086
      module load $rm || die "Failed to load modules: $rm"
    fi
  fi

  command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH (after module handling)"
  msg "DEBUG: Rscript=$(command -v Rscript)"
  msg "DEBUG: TMPDIR=$TMPDIR"

  if ! ensure_r_userlib_and_pkgs; then
    Rscript --vanilla -e 'cat("R pkg preflight failed.\n"); print(.libPaths())' >&2 || true
    die "R runtime missing required packages (data.table/R.utils/changepoint) and auto-install failed. If compute nodes lack CRAN access, pre-install once on a login node into: $R_LIBS_USER"
  fi

  Rscript --vanilla "$script" "$@"
}

# -----------------------------
# Execution helper (logging)
# -----------------------------
run_step() {
  local name="$1"; shift
  local logfile="$1"; shift
  mkdir -p "$(dirname "$logfile")"
  msg ">>> $name"
  { "$@"; } 2>&1 | tee "$logfile"
}

# -----------------------------
# Small utilities
# -----------------------------
need_hidden() {
  [[ -d "$HIDDEN_DIR" ]] || die "Missing hidden dir: $HIDDEN_DIR"
}

ensure_env_python() {
  if [[ -d "$CNV_ENV" ]]; then
    echo "$CONDA run -p $CNV_ENV --no-capture-output"
  else
    echo ""
  fi
}

pick_primary_corr_bw() {
  local corr_dir="$1"
  local sample="$2"
  local bins="$3"

  local cand=""
  cand="$(ls -1 "$corr_dir/${sample}.primary.corr.bins${bins}.bw" 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$corr_dir/${sample}.aligned.corr.bins${bins}.bw" 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }

  cand="$(ls -1 "$corr_dir"/*.bw 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

# NEW: pick BAMs for finalize counting
pick_primary_bam() {
  local bam_dir="$1"
  local sample="$2"

  local cand=""
  cand="$(ls -1 "$bam_dir/${sample}."*primary*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir/${sample}."*aligned*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }

  # fallback: any bam that looks "main-ish"
  cand="$(ls -1 "$bam_dir/${sample}."*.bam 2>/dev/null | grep -Ev 'split|disco|discord|supp|unmap|unmapped' | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }

  # last resort: first bam
  cand="$(ls -1 "$bam_dir"/*.bam 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

pick_split_bam() {
  local bam_dir="$1"
  local sample="$2"
  local cand=""
  cand="$(ls -1 "$bam_dir/${sample}."*split*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir"/*split*.bam 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

pick_disco_bam() {
  local bam_dir="$1"
  local sample="$2"
  local cand=""
  cand="$(ls -1 "$bam_dir/${sample}."*disco*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir/${sample}."*discord*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir"/*disco*.bam 2>/dev/null | head -n 1 || true)"
  [[ -n "$cand" ]] && { echo "$cand"; return 0; }
  cand="$(ls -1 "$bam_dir"/*discord*.bam 2>/dev/null | head -n 1 || true)"
  echo "$cand"
}

# -----------------------------
# Parse subcommand
# -----------------------------
if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

subcmd="$1"; shift

# -----------------------------
# Parse args (shared)
# -----------------------------
parse_common_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample) SAMPLE="${2:-}"; shift 2 ;;
      --genome) GENOME="${2:-}"; shift 2 ;;
      --gff) GFF="${2:-}"; shift 2 ;;
      --out-root) OUT_ROOT="${2:-}"; shift 2 ;;
      --bam) BAM="${2:-}"; shift 2 ;;
      --bam-dir) BAM_DIR="${2:-}"; shift 2 ;;
      --model) MODEL_IN="${2:-}"; shift 2 ;;

      --start-at) START_AT="${2:-}"; shift 2 ;;
      --stop-after) STOP_AFTER="${2:-}"; shift 2 ;;
      --force) FORCE=1; shift ;;
      --keep-tmp) KEEP_TMP=1; shift ;;

      --train-window) TRAIN_WINDOW="${2:-}"; shift 2 ;;
      --train-mapq) TRAIN_MAPQ="${2:-}"; shift 2 ;;
      --train-threads) TRAIN_THREADS="${2:-}"; shift 2 ;;
      --train-drop-dup) TRAIN_DROP_DUP=1; shift ;;
      --train-no-drop-dup) TRAIN_DROP_DUP=0; shift ;;
      --train-exclude-contigs) TRAIN_EXCLUDE_CONTIGS+=("${2:-}"); shift 2 ;;

      --apply-binsize) APPLY_BINSIZE="${2:-}"; shift 2 ;;
      --apply-threads) APPLY_THREADS="${2:-}"; shift 2 ;;
      --apply-emit-raw) APPLY_EMIT_RAW=1; shift ;;
      --apply-primary-prefer) APPLY_PRIMARY_PREFER="${2:-}"; shift 2 ;;

      --r-modules) R_MODULES="${2:-}"; shift 2 ;;

      --tile-bp) TILE_BP="${2:-}"; shift 2 ;;
      --min-tile-bp) MIN_TILE_BP="${2:-}"; shift 2 ;;
      --exclude-bed) EXCLUDE_BED="${2:-}"; shift 2 ;;
      --lambda) LAMBDA="${2:-}"; shift 2 ;;
      --min-probes-per-seg) MIN_PROBES_PER_SEG="${2:-}"; shift 2 ;;
      --flank) FLANK="${2:-}"; shift 2 ;;
      --agg) AGG="${2:-}"; shift 2 ;;
      --confidence-method) CONF_METHOD="${2:-}"; shift 2 ;;
      --final-weak-ratio) FINAL_WEAK_RATIO="${2:-}"; shift 2 ;;
      --final-strong-ratio) FINAL_STRONG_RATIO="${2:-}"; shift 2 ;;
      --final-weak-z) FINAL_WEAK_Z="${2:-}"; shift 2 ;;
      --final-strong-z) FINAL_STRONG_Z="${2:-}"; shift 2 ;;
      --final-keep-mode) FINAL_KEEP_MODE="${2:-}"; shift 2 ;;

      # finalize fusion
      --final-fuse) FINAL_FUSE=1; shift ;;
      --final-fuse-max-gap) FINAL_FUSE_MAX_GAP="${2:-}"; shift 2 ;;

      # finalize BAM counting
      --final-count-flank) FINAL_COUNT_FLANK="${2:-}"; shift 2 ;;
      --final-count-all-alignments) FINAL_COUNT_ALL_ALIGN=1; shift ;;
      --final-count-drop-dup) FINAL_COUNT_DROP_DUP=1; shift ;;
      --final-count-include-supp) FINAL_COUNT_INCLUDE_SUPP=1; shift ;;
      --final-count-include-secondary) FINAL_COUNT_INCLUDE_SECONDARY=1; shift ;;

      -h|--help) usage; exit 0 ;;
      *) die "Unknown argument: $1" ;;
    esac
  done
}

# -----------------------------
# Step implementations
# -----------------------------
step_cov_model_train() {
  [[ -n "$BAM" ]] || die "--bam is required for Model_Train"
  [[ -f "$BAM" ]] || die "BAM not found: $BAM"
  [[ -n "$GENOME" ]] || die "--genome is required"
  local out="$WORK/train"
  mkdir -p "$out"

  local cmd=( bash "$HIDDEN_DIR/cov_model_train.sh"
    --bam "$BAM"
    --genome "$GENOME"
    --out-dir "$out"
    --window "$TRAIN_WINDOW"
    --mapq "$TRAIN_MAPQ"
    --threads "$TRAIN_THREADS"
  )

  if [[ "$TRAIN_DROP_DUP" -eq 1 ]]; then
    cmd+=( --drop-dup-for-train )
  fi

  for x in "${TRAIN_EXCLUDE_CONTIGS[@]}"; do
    cmd+=( --exclude-contigs "$x" )
  done

  cmd+=( --r-modules "$R_MODULES" )
  run_step "Model_Train" "$LOG/cov_model_train.log" "${cmd[@]}"
}

step_cov_model_apply() {
  [[ -n "$BAM_DIR" ]] || die "--bam-dir is required for Model_Apply"
  [[ -d "$BAM_DIR" ]] || die "BAM dir not found: $BAM_DIR"
  [[ -n "$GENOME" ]] || die "--genome is required"

  local model="$MODEL_IN"
  if [[ -z "$model" ]]; then
    model="$WORK/train/model.rds"
  fi
  [[ -f "$model" ]] || die "Model .rds not found: $model (provide --model or run Model_Train)"

  local out="$WORK/apply"
  mkdir -p "$out"

  local cmd=( bash "$HIDDEN_DIR/cov_model_apply.sh"
    --dir "$BAM_DIR"
    --sample-key "$SAMPLE"
    --model "$model"
    --genome "$GENOME"
    --out-dir "$out"
    --binsize "$APPLY_BINSIZE"
    --threads "$APPLY_THREADS"
    --primary-prefer "$APPLY_PRIMARY_PREFER"
    --r-modules "$R_MODULES"
  )
  if [[ "$APPLY_EMIT_RAW" -eq 1 ]]; then
    cmd+=( --emit-raw )
  fi
  if [[ "$KEEP_TMP" -eq 1 ]]; then
    cmd+=( --keep-tmp )
  fi

  run_step "Model_Apply" "$LOG/cov_model_apply.log" "${cmd[@]}"

  mkdir -p "$SAMPLE_DIR/corr_bw" "$SAMPLE_DIR/factors"
  rsync -a "$out/$SAMPLE/corr_bw/" "$SAMPLE_DIR/corr_bw/" || true
  rsync -a "$out/$SAMPLE/factors/" "$SAMPLE_DIR/factors/" || true
  if [[ -d "$out/$SAMPLE/raw_bw" ]]; then
    mkdir -p "$SAMPLE_DIR/raw_bw"
    rsync -a "$out/$SAMPLE/raw_bw/" "$SAMPLE_DIR/raw_bw/" || true
  fi
}

step_cov_bw_qc() {
  local out_tsv="$WORK/cov_bw_qc.summary.tsv"
  run_step "Model_Finalize" "$LOG/cov_bw_qc.log" bash "$HIDDEN_DIR/cov_bw_qc.sh" \
    --outdir "$WORK/apply" \
    --out "$out_tsv" \
    --bam-root "$BAM_DIR"
  cp "$WORK/cov_bw_qc.summary.tsv" "$SAMPLE_DIR/$SAMPLE.cov_bw_qc.summary.tsv"
}

step_make_probes() {
  [[ -n "$GENOME" ]] || die "--genome required"
  [[ -n "$GFF" ]] || die "--gff required"

  local corr_dir="$SAMPLE_DIR/corr_bw"
  [[ -d "$corr_dir" ]] || die "Missing corr_bw dir: $corr_dir (run Model_Apply first)"

  local primary_bw
  primary_bw="$(pick_primary_corr_bw "$corr_dir" "$SAMPLE" "$APPLY_BINSIZE")"
  [[ -n "$primary_bw" && -f "$primary_bw" ]] || die "Could not select primary corrected BW under: $corr_dir"
  msg "Primary BW selected: $primary_bw"

  local out="$WORK/probes"
  mkdir -p "$out"
  local py="$HIDDEN_DIR/make_probes.py"
  [[ -f "$py" ]] || die "Missing make_probes.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local cmd=( python "$py"
    --genome "$GENOME"
    --gff "$GFF"
    --primary-bw "$primary_bw"
    --out-prefix "$out"
    --tile-bp "$TILE_BP"
    --min-tile-bp "$MIN_TILE_BP"
  )
  if [[ -n "$EXCLUDE_BED" ]]; then
    cmd+=( --exclude-bed "$EXCLUDE_BED" )
  fi

  if [[ -n "$runner" ]]; then
    run_step "CNV_Probe" "$LOG/make_probes.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Probe" "$LOG/make_probes.log" "${cmd[@]}"
  fi

}

step_fused_lasso_segment() {
  local probe_tsv="$WORK/probes/probe_table.tsv.gz"
  [[ -f "$probe_tsv" ]] || die "Missing probe table (run CNV_Probe first): $probe_tsv"

  local out="$WORK/segment"
  mkdir -p "$out"

  local r_script="$HIDDEN_DIR/fused_lasso_segment.R"
  [[ -f "$r_script" ]] || die "Missing R script on this node: $r_script"

  msg "Running fused_lasso_segment (PELT) penalty(lambda)=$LAMBDA min_probes=$MIN_PROBES_PER_SEG"
  msg "R modules mode: ${R_MODULES:-off}"
  msg "R TMPDIR: $WORK/tmp_R"

  run_step "CNV_Segment" "$LOG/fused_lasso_segment.log" \
    run_r_file_hpc "$r_script" \
      --probe-tsv "$probe_tsv" \
      --out-prefix "$out" \
      --lambda "$LAMBDA" \
      --min-probes-per-seg "$MIN_PROBES_PER_SEG" \
      --tmpdir "$WORK/tmp_R"

  cp "$out/segments.tsv.gz" "$SAMPLE_DIR/$SAMPLE.segments.tsv.gz"
  cp "$out/segments.bed.gz" "$SAMPLE_DIR/$SAMPLE.segments.bed.gz"

}

step_boundary_support() {
  local boundaries="$WORK/segment/boundaries.bed.gz"
  [[ -f "$boundaries" ]] || die "Missing boundaries (run CNV_Segment first): $boundaries"

  local corr_dir="$SAMPLE_DIR/corr_bw"
  [[ -d "$corr_dir" ]] || die "Missing corr_bw dir: $corr_dir"

  local split_bw disco_bw
  split_bw="$(ls -1 "$corr_dir/${SAMPLE}."*split*"corr.bins${APPLY_BINSIZE}.bw" 2>/dev/null | head -n 1 || true)"
  disco_bw="$(ls -1 "$corr_dir/${SAMPLE}."*disco*"corr.bins${APPLY_BINSIZE}.bw" 2>/dev/null | head -n 1 || true)"
  [[ -n "$split_bw" && -f "$split_bw" ]] || die "Could not find splitters corr BW in: $corr_dir"
  [[ -n "$disco_bw" && -f "$disco_bw" ]] || die "Could not find discordant corr BW in: $corr_dir"

  msg "Split BW selected: $split_bw"
  msg "Disco BW selected: $disco_bw"

  local out="$WORK/boundary_support"
  mkdir -p "$out"
  local py="$HIDDEN_DIR/boundary_support_from_bigwigs.py"
  [[ -f "$py" ]] || die "Missing boundary_support_from_bigwigs.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_tsv="$out/boundary_support.tsv.gz"
  local out_bg="$out/boundary_support.z.bedGraph.gz"

  local cmd=( python "$py"
    --boundaries-bed "$boundaries"
    --genome "$GENOME"
    --split-bw "$split_bw"
    --disco-bw "$disco_bw"
    --flank "$FLANK"
    --agg "$AGG"
    --out-tsv "$out_tsv"
    --out-bedgraph "$out_bg"
    --bedgraph-field combined_z_chr
  )

  if [[ -n "$runner" ]]; then
    run_step "CNV_Boundary" "$LOG/boundary_support_from_bigwigs.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Boundary" "$LOG/boundary_support_from_bigwigs.log" "${cmd[@]}"
  fi

}

step_segment_confidence() {
  local seg_bed="$WORK/segment/segments.bed.gz"
  local bsup="$WORK/boundary_support/boundary_support.tsv.gz"
  [[ -f "$seg_bed" ]] || die "Missing segments.bed.gz: $seg_bed"
  [[ -f "$bsup" ]] || die "Missing boundary_support.tsv.gz: $bsup"

  local out="$WORK/confidence"
  mkdir -p "$out"
  local py="$HIDDEN_DIR/segment_confidence_from_boundaries.py"
  [[ -f "$py" ]] || die "Missing segment_confidence_from_boundaries.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_tsv="$out/${SAMPLE}_segments.confidence.tsv.gz"
  local out_bed="$out/${SAMPLE}_segments.confidence.bed.gz"

  local cmd=( python "$py"
    --segments-bed "$seg_bed"
    --boundary-support-tsv "$bsup"
    --out-tsv "$out_tsv"
    --out-bed "$out_bed"
    --method "$CONF_METHOD"
    --score-from z
  )

  if [[ -n "$runner" ]]; then
    run_step "CNV_Confidence" "$LOG/segment_confidence_from_boundaries.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Confidence" "$LOG/segment_confidence_from_boundaries.log" "${cmd[@]}"
  fi
}

step_cnv_finalize() {
  local seg_bed="$WORK/segment/segments.bed.gz"
  local bsup="$WORK/boundary_support/boundary_support.tsv.gz"
  [[ -f "$seg_bed" ]] || die "Missing segments.bed.gz: $seg_bed"
  [[ -f "$bsup" ]] || die "Missing boundary_support.tsv.gz: $bsup"

  local out="$WORK/final"
  mkdir -p "$out"

  local py="$HIDDEN_DIR/cnv_finalize_segments.py"
  [[ -f "$py" ]] || die "Missing cnv_finalize_segments.py: $py"

  local runner
  runner="$(ensure_env_python)"

  local out_calls="$out/segments.calls.tsv.gz"
  local out_y="$out/${SAMPLE}_cn.y.bedGraph.gz"
  local out_ratio="$out/${SAMPLE}_cn.ratio.bedGraph.gz"

  # ---- Pick BAMs for BAM-based counting (finalize) ----
  local primary_bam split_bam disco_bam
  primary_bam="$BAM"
  if [[ -z "$primary_bam" && -n "$BAM_DIR" && -d "$BAM_DIR" ]]; then
    primary_bam="$(pick_primary_bam "$BAM_DIR" "$SAMPLE")"
  fi

  split_bam=""
  disco_bam=""
  if [[ -n "$BAM_DIR" && -d "$BAM_DIR" ]]; then
    split_bam="$(pick_split_bam "$BAM_DIR" "$SAMPLE")"
    disco_bam="$(pick_disco_bam "$BAM_DIR" "$SAMPLE")"
  fi

  # Only pass BAM args if we actually found files (pysam will require .bai)
  local cmd=( python "$py"
    --segments-bed "$seg_bed"
    --boundary-support-tsv "$bsup"
    --out-tsv "$out_calls"
    --weak-ratio "$FINAL_WEAK_RATIO"
    --strong-ratio "$FINAL_STRONG_RATIO"
    --weak-z "$FINAL_WEAK_Z"
    --strong-z "$FINAL_STRONG_Z"
    --keep-mode "$FINAL_KEEP_MODE"
    --cn-base 1.0
    --out-bedgraph-y "$out_y"
    --out-bedgraph-ratio "$out_ratio"
  )

  # Fusion controls
  if [[ "$FINAL_FUSE" -eq 1 ]]; then
    cmd+=( --fuse --fuse-max-gap "$FINAL_FUSE_MAX_GAP" )
  fi

  # BAM-based read counts (segment + boundary windows)
  if [[ -n "$primary_bam" && -f "$primary_bam" ]]; then
    msg "Finalize counts: primary_bam=$primary_bam"
    cmd+=( --primary-bam "$primary_bam" )
  else
    msg "WARN: finalize counts: primary_bam not found (skipping seg_primary_reads)"
  fi

  if [[ -n "${split_bam:-}" && -f "${split_bam:-}" ]]; then
    msg "Finalize counts: split_bam=$split_bam"
    cmd+=( --split-bam "$split_bam" )
  else
    msg "WARN: finalize counts: split_bam not found (skipping split counts)"
  fi

  if [[ -n "${disco_bam:-}" && -f "${disco_bam:-}" ]]; then
    msg "Finalize counts: disco_bam=$disco_bam"
    cmd+=( --disco-bam "$disco_bam" )
  else
    msg "WARN: finalize counts: disco_bam not found (skipping disco counts)"
  fi

  cmd+=( --count-flank "$FINAL_COUNT_FLANK" )

  if [[ "$FINAL_COUNT_ALL_ALIGN" -eq 1 ]]; then
    cmd+=( --count-all-alignments )
  fi
  if [[ "$FINAL_COUNT_DROP_DUP" -eq 1 ]]; then
    cmd+=( --count-drop-dup )
  fi
  if [[ "$FINAL_COUNT_INCLUDE_SUPP" -eq 1 ]]; then
    cmd+=( --count-include-supp )
  fi
  if [[ "$FINAL_COUNT_INCLUDE_SECONDARY" -eq 1 ]]; then
    cmd+=( --count-include-secondary )
  fi

  if [[ -n "$runner" ]]; then
    run_step "CNV_Finalize" "$LOG/cnv_finalize_segments.log" $runner "${cmd[@]}"
  else
    run_step "CNV_Finalize" "$LOG/cnv_finalize_segments.log" "${cmd[@]}"
  fi

for f in "$WORK/final/"/*; do
  [[ -f "$f" ]] || continue
  bn="$(basename "$f")"
  cp "$f" "$SAMPLE_DIR/${SAMPLE}.${bn}"
done

  
  
}

run_internal_step() {
  local internal="$1"
  case "$internal" in
    cov_model_train) step_cov_model_train ;;
    cov_model_apply) step_cov_model_apply ;;
    cov_bw_qc) step_cov_bw_qc ;;
    make_probes) step_make_probes ;;
    fused_lasso_segment) step_fused_lasso_segment ;;
    boundary_support_from_bigwigs) step_boundary_support ;;
    segment_confidence_from_boundaries) step_segment_confidence ;;
    cnv_finalize_segments) step_cnv_finalize ;;
    *) die "Unknown internal step: $internal" ;;
  esac
}

# -----------------------------
# Subcommand: step
# -----------------------------
if [[ "$subcmd" == "step" ]]; then
  [[ $# -ge 1 ]] || die "step requires a step name (alias or internal)"
  step_name="$1"; shift

  parse_common_args "$@"

  [[ -n "$SAMPLE" ]] || die "--sample required"
  [[ -n "$OUT_ROOT" ]] || die "--out-root required"

  internal="$(resolve_step "$step_name")"
  [[ -n "$internal" ]] || die "Unknown step: $step_name"

  SAMPLE_DIR="$OUT_ROOT/$SAMPLE/$SAMPLE"
  WORK="$SAMPLE_DIR/_work"
  LOG="$SAMPLE_DIR/_log"
  mkdir -p "$SAMPLE_DIR" "$WORK" "$LOG"

  msg "Sample        : $SAMPLE"
  msg "Output root   : $OUT_ROOT"
  msg "Sample dir    : $SAMPLE_DIR"
  msg "Step          : $step_name -> $internal"

  run_internal_step "$internal"
  exit 0
fi

# -----------------------------
# Subcommand: run
# -----------------------------
if [[ "$subcmd" != "run" ]]; then
  die "Unknown subcommand: $subcmd (expected run|step)"
fi

parse_common_args "$@"

[[ -n "$SAMPLE" ]] || die "--sample required"
[[ -n "$GENOME" ]] || die "--genome required"
[[ -n "$OUT_ROOT" ]] || die "--out-root required"

START_AT="$(resolve_step "$START_AT")"
STOP_AFTER="$(resolve_step "$STOP_AFTER")"
[[ -n "$START_AT" ]] || die "Unknown --start-at step"
[[ -n "$STOP_AFTER" ]] || die "Unknown --stop-after step"

SAMPLE_DIR="$OUT_ROOT/$SAMPLE/$SAMPLE"
WORK="$SAMPLE_DIR/_work"
LOG="$SAMPLE_DIR/_log"
mkdir -p "$SAMPLE_DIR" "$WORK" "$LOG"

msg "Sample        : $SAMPLE"
msg "Output root   : $OUT_ROOT"
msg "Sample dir    : $SAMPLE_DIR"
msg "Start at      : $START_AT"
msg "Stop after   : $STOP_AFTER"

sidx="$(idx "$START_AT")"
eidx="$(idx "$STOP_AFTER")"
[[ "$sidx" -ge 0 ]] || die "Bad start step: $START_AT"
[[ "$eidx" -ge 0 ]] || die "Bad stop step: $STOP_AFTER"
[[ "$eidx" -ge "$sidx" ]] || die "--stop-after must be >= --start-at in pipeline order"

for i in $(seq "$sidx" "$eidx"); do
  internal="${steps[$i]}"

  if [[ "$FORCE" -ne 1 ]]; then
    case "$internal" in
      cov_model_train)
        [[ -f "$WORK/train/model.rds" ]] && { msg "Skip cov_model_train (exists)"; continue; }
        ;;
      cov_model_apply)
        [[ -d "$SAMPLE_DIR/corr_bw" && -n "$(ls -1 "$SAMPLE_DIR/corr_bw"/*.bw 2>/dev/null | head -n 1 || true)" ]] && { msg "Skip cov_model_apply (corr_bw exists)"; continue; }
        ;;
      cov_bw_qc)
        [[ -f "$SAMPLE_DIR/cov_bw_qc.summary.tsv" ]] && { msg "Skip cov_bw_qc (summary exists)"; continue; }
        ;;
      make_probes)
        [[ -f "$WORK/probes/probe_table.tsv.gz" ]] && { msg "Skip make_probes (probe table exists)"; continue; }
        ;;
      fused_lasso_segment)
        [[ -f "$WORK/segment/segments.bed.gz" ]] && { msg "Skip fused_lasso_segment (segments exist)"; continue; }
        ;;
      boundary_support_from_bigwigs)
        [[ -f "$WORK/boundary_support/boundary_support.tsv.gz" ]] && { msg "Skip boundary_support (exists)"; continue; }
        ;;
      segment_confidence_from_boundaries)
        [[ -f "$WORK/confidence/${SAMPLE}_segments.confidence.tsv.gz" ]] && { msg "Skip segment_confidence (exists)"; continue; }
        ;;
      cnv_finalize_segments)
        [[ -f "$WORK/final/segments.calls.tsv.gz" ]] && { msg "Skip cnv_finalize (calls exist)"; continue; }
        ;;
    esac
  fi

  run_internal_step "$internal"
done

msg "CLOverCNV run complete: $SAMPLE_DIR"
