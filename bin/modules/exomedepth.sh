#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# exomedepth_call.sh
#   - HPC module-R runner (R via Lmod)
#   - avoids rtracklayer entirely
#   - GFF->BED conversion in bash
#   - R preflight/install via heredoc (no quoting bugs)
#   - AUTO selects best R/* and loads prereqs as spider requires
# ============================================================

msg() { echo "[exomedepth_call] $*"; }
die() { msg "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Usage:
  mapi modules exomedepth_call \
    --test-bam test.bam \
    --ref-bams ref1.bam,ref2.bam[,ref3.bam...] \
    --out-dir /path/to/out \
    (--exons-bed exons.bed | --gff anno.gff[.gz]) \
    [--features exon,CDS] \
    [--sample NAME] \
    [--r-modules off|auto|R/4.3.1|"gcc/... openmpi/... R/4.3.1"] \
    [--r-autoinstall 1|0]

Notes:
  - Avoids rtracklayer (prevents XML/RCurl/libxml build chain).
  - Uses Lmod module R; AUTO selects the newest R/4.3.* if available, else newest R/*.
  - Loads prerequisite toolchain modules as required by `module spider R/<ver>`.

Outputs:
  exomedepth.calls.tsv
  exomedepth.calls.bed
  exomedepth.log
  exons.from_gff.bed (if --gff used)
USAGE
}

# -----------------------------
# Args
# -----------------------------
TEST_BAM=""
REF_BAMS_CSV=""
EXONS_BED=""
GFF=""
FEATURES="exon,CDS"
OUT_DIR=""
SAMPLE=""

R_MODULES="auto"
R_AUTOINSTALL=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --test-bam) TEST_BAM="$2"; shift 2 ;;
    --ref-bams) REF_BAMS_CSV="$2"; shift 2 ;;
    --exons-bed) EXONS_BED="$2"; shift 2 ;;
    --gff) GFF="$2"; shift 2 ;;
    --features) FEATURES="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --r-modules) R_MODULES="$2"; shift 2 ;;
    --r-autoinstall) R_AUTOINSTALL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$TEST_BAM" ]] || die "--test-bam is required"
[[ -n "$REF_BAMS_CSV" ]] || die "--ref-bams is required"
[[ -n "$OUT_DIR" ]] || die "--out-dir is required"
[[ -f "$TEST_BAM" ]] || die "Test BAM not found: $TEST_BAM"
[[ -f "${TEST_BAM}.bai" || -f "${TEST_BAM%.bam}.bai" ]] || die "Index not found for test BAM: $TEST_BAM"

mkdir -p "$OUT_DIR"

if [[ -z "${SAMPLE:-}" ]]; then
  SAMPLE="$(basename "$TEST_BAM")"
  SAMPLE="${SAMPLE%.bam}"
fi

IFS="," read -r -a REF_BAMS <<< "${REF_BAMS_CSV:-}"
[[ ${#REF_BAMS[@]} -ge 1 ]] || die "No ref BAMs parsed from --ref-bams"
for b in "${REF_BAMS[@]}"; do
  [[ -f "$b" ]] || die "Ref BAM not found: $b"
  [[ -f "${b}.bai" || -f "${b%.bam}.bai" ]] || die "Index not found for ref BAM: $b"
done

WORK="$OUT_DIR"
mkdir -p "$WORK/tmp_R"

LOG="$OUT_DIR/exomedepth.log"
TSV="$OUT_DIR/exomedepth.calls.tsv"
BED="$OUT_DIR/exomedepth.calls.bed"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_HELPER="$SCRIPT_DIR/.exomedepth_call/run_exomedepth.R"
[[ -f "$R_HELPER" ]] || die "Missing helper R script: $R_HELPER"

# ------------------------------------------------------------
# If --gff provided, build exon-like BED into OUT_DIR
# ------------------------------------------------------------
if [[ -n "${GFF:-}" ]]; then
  [[ -f "$GFF" ]] || die "GFF not found: $GFF"
  EXONS_BED="$OUT_DIR/exons.from_gff.bed"

  msg "building exons BED from GFF: $GFF"
  msg "features: $FEATURES"

  features_re="$(
    printf '%s\n' "$FEATURES" \
      | awk -F',' '{
          for (i=1; i<=NF; i++) {
            gsub(/^[ \t]+|[ \t]+$/, "", $i)
            printf "%s%s", (i==1 ? "" : "|"), $i
          }
        }'
  )"
  fre="^(${features_re})$"

  gff_cat() {
    case "$1" in
      *.gz) gzip -cd -- "$1" ;;
      *)    cat -- "$1" ;;
    esac
  }

  gff_cat "$GFF" \
    | awk -v OFS="\t" -v fre="$fre" '
        $0 ~ /^#/ { next }
        NF < 9 { next }
        $3 ~ fre {
          chrom=$1
          start=$4-1; if(start<0) start=0
          end=$5
          attr=$9

          id="."
          if (match(attr, /(^|;)ID=([^;]+)/, m)) id=m[2]
          else if (match(attr, /(^|;)Parent=([^;]+)/, m)) id=m[2]
          else if (match(attr, /(^|;)Name=([^;]+)/, m)) id=m[2]

          print chrom, start, end, id
        }
      ' \
    | sort -k1,1 -k2,2n -k3,3n \
    > "$EXONS_BED"

  [[ -s "$EXONS_BED" ]] || die "GFF->BED produced empty file: $EXONS_BED"
fi

[[ -n "${EXONS_BED:-}" ]] || die "Provide --exons-bed or --gff"
[[ -f "$EXONS_BED" ]] || die "Exons BED not found: $EXONS_BED"

# ============================================================
# Lmod helpers
# ============================================================
maybe_enable_modules() {
  if command -v module >/dev/null 2>&1; then return 0; fi
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh && command -v module >/dev/null 2>&1 && return 0
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash && command -v module >/dev/null 2>&1 && return 0
  return 1
}

list_r_modules() {
  (module spider R 2>&1 || true) \
    | awk '
        /Versions:/{flag=1; next}
        flag && $1 ~ /^R\/[0-9]/ {print $1}
      ' \
    | sed 's/[[:space:]]*$//' \
    | sort -V | uniq
}

pick_best_r_module() {
  mapfile -t mods < <(list_r_modules)
  [[ "${#mods[@]}" -gt 0 ]] || return 1

  local best43
  best43="$(printf "%s\n" "${mods[@]}" | awk -F/ '$2 ~ /^4\.3(\.|$)/ {print $0}' | sort -V | tail -n 1)"
  if [[ -n "$best43" ]]; then
    echo "$best43"
    return 0
  fi
  printf "%s\n" "${mods[@]}" | sort -V | tail -n 1
}

# Parse the prereq line(s) from `module spider R/<ver>`
# Example spider output:
#   You will need to load all module(s) on any one of the lines below ...
#     gcc/11.4.0  openmpi/4.1.4
extract_prereq_lines() {
  local rmod="$1"
  (module spider "$rmod" 2>&1 || true) | awk '
    function trim(s){ sub(/^[ \t]+/,"",s); sub(/[ \t]+$/,"",s); return s }
    /You will need to load all module\(s\) on any one of the lines below/ {flag=1; next}
    flag && /^[ \t]*$/ {exit}
    flag && $0 ~ /^[ \t]+/ {
      line = trim($0)
      if (line != "") print line
    }
  ' | sed '/^R\/[0-9]/d' | sed '/^This module provides/d' | sed '/^Description/d' | sed '/^Extensions/d' | sed '/^$/d'
}

load_r_with_prereqs() {
  local rmod="$1"

  # Try direct
  if module load "$rmod" >/dev/null 2>&1; then
    return 0
  fi

  mapfile -t prereq < <(extract_prereq_lines "$rmod" || true)

  # If spider gave us nothing, fail loudly (better than silently doing wrong thing)
  if [[ "${#prereq[@]}" -eq 0 ]]; then
    return 1
  fi

  # Try each prereq line (usually just one)
  for line in "${prereq[@]}"; do
    msg "Loading prereqs for $rmod: $line"
    module purge >/dev/null 2>&1 || true
    # shellcheck disable=SC2086
    if module load $line >/dev/null 2>&1; then
      if module load "$rmod" >/dev/null 2>&1; then
        return 0
      fi
    fi
  done
  return 1
}

# ============================================================
# R library + preflight
# ============================================================
pick_r_user_lib_root() {
  if [[ -n "${MAPI_R_LIBS_USER_ROOT:-}" ]]; then
    echo "$MAPI_R_LIBS_USER_ROOT"; return 0
  fi
  if [[ -d "/standard" && -n "${HPC_ALLOCATION:-}" ]]; then
    local cand="/standard/${HPC_ALLOCATION}/${USER}/R/goolf"
    if mkdir -p "$cand" >/dev/null 2>&1; then
      echo "$cand"; return 0
    fi
  fi
  echo "$HOME/R/goolf"
}

r_preflight_install() {
  Rscript --vanilla - <<'RS' >> "$LOG" 2>&1
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))

if (!requireNamespace("data.table", quietly=TRUE)) {
  install.packages("data.table", lib=lib)
}
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", lib=lib)
}

# DO NOT install rtracklayer. Keep the dependency surface minimal.
pkgs <- c("ExomeDepth","GenomicRanges","Rsamtools","IRanges","S4Vectors")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(need)) {
  message("[R preflight] installing Bioc: ", paste(need, collapse=", "))
  BiocManager::install(need, lib=lib, ask=FALSE, update=FALSE)
} else {
  message("[R preflight] ok")
}
RS
}

r_preflight_verify() {
  Rscript --vanilla - <<'RS' >> "$LOG" 2>&1
lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(lib, .libPaths()))
req <- c("ExomeDepth","GenomicRanges","Rsamtools")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  cat("MISSING:", paste(miss, collapse=" "), "\n")
  print(.libPaths())
  quit(status=2)
} else {
  cat("OK\n")
}
RS
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
    r_preflight_install || return 1
  fi
  r_preflight_verify || return 1
}

run_r_file_hpc() {
  local script="$1"; shift
  [[ -f "$script" ]] || die "R script not found: $script"

  export TMPDIR="$WORK/tmp_R"
  mkdir -p "$TMPDIR"

  local rm
  rm="$(echo "${R_MODULES:-auto}" | awk '{$1=$1;print}')"
  rm="${rm%\"}"; rm="${rm#\"}"
  rm="${rm%\'}"; rm="${rm#\'}"

  maybe_enable_modules || die "module system not available on this node"
  module purge >/dev/null 2>&1 || true

  case "${rm,,}" in
    ""|"off"|"none")
      ;;
    "auto")
      local rpick
      rpick="$(pick_best_r_module || true)"
      [[ -n "$rpick" ]] || die "Could not find any R/<version> via 'module spider R'"
      msg "Auto-selected R module: $rpick"
      load_r_with_prereqs "$rpick" || die "Failed to load $rpick with prerequisites. See: module spider $rpick"
      ;;
    R/*)
      msg "Loading R module: $rm"
      load_r_with_prereqs "$rm" || die "Failed to load $rm with prerequisites. See: module spider $rm"
      ;;
    *)
      # Treat as explicit full module line(s), e.g. "gcc/... openmpi/... R/4.3.1"
      msg "Loading explicit module line(s): $rm"
      # shellcheck disable=SC2086
      module load $rm >/dev/null 2>&1 || die "Failed to load explicit modules: $rm"
      ;;
  esac

  command -v Rscript >/dev/null 2>&1 || die "Rscript not found in PATH (after module handling)"
  msg "DEBUG: Rscript=$(command -v Rscript)"
  msg "DEBUG: TMPDIR=$TMPDIR"

  if ! ensure_r_userlib_and_pkgs; then
    Rscript --vanilla -e 'cat("R pkg preflight failed.\n"); print(.libPaths())' >> "$LOG" 2>&1 || true
    die "R runtime missing required packages (ExomeDepth/GenomicRanges/Rsamtools) and auto-install failed. See $LOG"
  fi

  Rscript --vanilla "$script" "$@"
}

# -----------------------------
# Run
# -----------------------------
msg "sample    : $SAMPLE"
msg "test bam  : $TEST_BAM"
msg "ref bams  : ${#REF_BAMS[@]}"
msg "exons bed : $EXONS_BED"
msg "out       : $OUT_DIR"
msg "r-modules : $R_MODULES"
msg "r-autoinstall : $R_AUTOINSTALL"

set -x
run_r_file_hpc "$R_HELPER" \
  --sample "$SAMPLE" \
  --test-bam "$TEST_BAM" \
  --ref-bams "$REF_BAMS_CSV" \
  --exons-bed "$EXONS_BED" \
  --out-tsv "$TSV" \
  --out-bed "$BED" \
  >> "$LOG" 2>&1
set +x

msg "done"


