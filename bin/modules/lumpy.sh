#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# lumpy_call.sh
# ============================================================

# -----------------------------
# Logging helpers
# -----------------------------
say() { echo "[lumpy_call] $*"; }
die() { say "ERROR: $*"; exit 1; }

# -----------------------------
# Help
# -----------------------------
usage() {
  cat <<'EOF'
Usage:
  mapi modules lumpy_call \
    --bam sample.bam \
    --out-dir /path/to/out \
    [--sample NAME] \
    [--threads 8]

Notes:
  - Expects a coordinate-sorted, indexed BAM (.bam + .bai).
  - Uses lumpyexpress to produce a VCF of breakpoint calls.
  - This is not a pure read-depth CNV caller, but good as a breakpoint-based comparator.

Outputs (in --out-dir):
  lumpy.vcf
  lumpy.log
EOF
}

# -----------------------------
# Args
# -----------------------------
BAM=""
OUT_DIR=""
SAMPLE=""
THREADS=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$BAM" ]] || die "--bam is required"
[[ -n "$OUT_DIR" ]] || die "--out-dir is required"
[[ -f "$BAM" ]] || die "BAM not found: $BAM"
[[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" ]] || die "BAM index (.bai) not found for: $BAM"

mkdir -p "$OUT_DIR"
TMP="$OUT_DIR/_tmp"
mkdir -p "$TMP"

if [[ -z "$SAMPLE" ]]; then
  SAMPLE="$(basename "$BAM")"
  SAMPLE="${SAMPLE%.bam}"
fi

OUT_VCF="$OUT_DIR/lumpy.vcf"
LOG="$OUT_DIR/lumpy.log"

say "bam      : $BAM"
say "sample   : $SAMPLE"
say "threads  : $THREADS"
say "out      : $OUT_VCF"

# lumpyexpress uses its own extraction helpers internally when given -B.
# Keep it simple: feed BAM, get VCF.
# If your lumpy build requires explicit split/discordant BAMs, we can adapt later.
set -x
lumpyexpress \
  -B "$BAM" \
  -o "$OUT_VCF" \
  -T "$TMP" \
  2> "$LOG"
set +x

say "done"
