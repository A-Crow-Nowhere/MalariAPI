#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# cnvnator_call.sh
# ============================================================

say() { echo "[cnvnator_call] $*"; }
die() { say "ERROR: $*"; exit 1; }

usage() {
  cat <<'EOF'
Usage:
  mapi modules cnvnator_call \
    --bam sample.bam \
    --ref reference.fa \
    --out-dir /path/to/out \
    [--bin 100] \
    [--sample NAME] \
    [--refdir DIR]

Notes:
  - CNVnator is read-depth based.
  - Requires coordinate-sorted, indexed BAM (.bam + .bai).
  - CNVnator expects -d to point to a directory containing one FASTA per contig,
    named exactly like the contigs in the BAM header (e.g. PfDd2_01.fa).
    This module will auto-build that directory from --ref unless --refdir is provided.

Outputs (in --out-dir):
  cnvnator.<sample>.root
  cnvnator.<sample>.calls.tsv
  cnvnator.<sample>.calls.bed
  cnvnator.<sample>.log
  cnvnator_refdir/   (if auto-built)
EOF
}

# -----------------------------
# Args
# -----------------------------
BAM=""
REF=""
OUT_DIR=""
BIN=100
SAMPLE=""
REFDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --bin) BIN="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --refdir) REFDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$BAM" ]] || die "--bam is required"
[[ -n "$REF" ]] || die "--ref is required"
[[ -n "$OUT_DIR" ]] || die "--out-dir is required"
[[ -f "$BAM" ]] || die "BAM not found: $BAM"
[[ -f "$REF" ]] || die "REF not found: $REF"

# BAM index can be either sample.bam.bai or sample.bai
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  die "BAM index (.bai) not found for: $BAM"
fi

mkdir -p "$OUT_DIR"

if [[ -z "$SAMPLE" ]]; then
  SAMPLE="$(basename "$BAM")"
  SAMPLE="${SAMPLE%.bam}"
fi

ROOT="$OUT_DIR/cnvnator.${SAMPLE}.root"
LOG="$OUT_DIR/cnvnator.${SAMPLE}.log"
CALLS="$OUT_DIR/cnvnator.${SAMPLE}.calls.tsv"
BED="$OUT_DIR/cnvnator.${SAMPLE}.calls.bed"

TMP="$OUT_DIR/_tmp"
mkdir -p "$TMP"
export TMPDIR="$TMP"

say "bam    : $BAM"
say "ref    : $REF"
say "bin    : $BIN"
say "sample : $SAMPLE"
say "out    : $OUT_DIR"

# Infer chrom list from BAM header
CHROMS_FILE="$OUT_DIR/_chroms.${SAMPLE}.txt"
samtools view -H "$BAM" \
  | awk -F'\t' '$1=="@SQ"{for(i=1;i<=NF;i++) if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}' \
  > "$CHROMS_FILE"
[[ -s "$CHROMS_FILE" ]] || die "Failed to infer chrom list from BAM header"

CHROMS="$(paste -sd' ' "$CHROMS_FILE")"

# Build CNVnator refdir if not provided
if [[ -z "$REFDIR" ]]; then
  REFDIR="$OUT_DIR/cnvnator_refdir"
fi
mkdir -p "$REFDIR"

# Ensure reference indexed
if [[ ! -f "${REF}.fai" ]]; then
  say "indexing reference: $REF"
  samtools faidx "$REF"
fi

# Populate per-contig FASTAs
say "building per-contig FASTAs in: $REFDIR"
while read -r c; do
  [[ -n "$c" ]] || continue
  fa="$REFDIR/$c.fa"
  if [[ ! -s "$fa" ]]; then
    samtools faidx "$REF" "$c" > "$fa" || die "failed to extract contig '$c' from $REF"
  fi
done < "$CHROMS_FILE"

# CNVnator pipeline
set -x
cnvnator -root "$ROOT" -tree "$BAM" -chrom $CHROMS 2> "$LOG"
cnvnator -root "$ROOT" -his "$BIN" -d "$REFDIR" >> "$LOG" 2>&1
cnvnator -root "$ROOT" -stat "$BIN" >> "$LOG" 2>&1
cnvnator -root "$ROOT" -partition "$BIN" >> "$LOG" 2>&1
cnvnator -root "$ROOT" -call "$BIN" > "$CALLS"
set +x

# Make a pragmatic BED from cnvnator calls
# CNVnator call lines typically include "chr:start-end" token; we extract first such token.
awk 'BEGIN{OFS="\t"}
{
  for(i=1;i<=NF;i++){
    if($i ~ /:/ && $i ~ /-/){
      split($i,a,":"); chrom=a[1];
      split(a[2],b,"-");
      start=b[1]-1; end=b[2];
      if(start<0) start=0;
      # name=type plus a little metadata; keep full line in last column for debugging
      print chrom, start, end, $1, $0;
      break
    }
  }
}' "$CALLS" > "$BED" || true

say "done"
