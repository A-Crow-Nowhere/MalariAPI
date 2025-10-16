#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
mapi run --sample-dir <dir> --ref <ref.fa> [--threads N] [--prefix NAME]
EOF
  exit 1
}

SAMPLEDIR=""; REF=""; THREADS=8; PREFIX=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-dir) SAMPLEDIR="$2"; shift 2 ;;
    --ref)        REF="$2";        shift 2 ;;
    --threads)    THREADS="$2";    shift 2 ;;
    --prefix)     PREFIX="$2";     shift 2 ;;
    -h|--help)    usage ;;
    *) echo "Unknown arg: $1"; usage ;;
  esac
done

[[ -n "$SAMPLEDIR" && -d "$SAMPLEDIR" ]] || { echo "[mapi] sample dir missing"; usage; }
[[ -n "$REF" && -f "$REF" ]] || { echo "[mapi] ref fasta missing"; usage; }

shopt -s nullglob
readarray -t R1s < <(ls "$SAMPLEDIR"/*_R1*.fastq* 2>/dev/null || true)
readarray -t R2s < <(ls "$SAMPLEDIR"/*_R2*.fastq* 2>/dev/null || true)
readarray -t SEs < <(ls "$SAMPLEDIR"/*.fastq* 2>/dev/null || true)

paired=0; single=0; R1=""; R2=""; SE=""
if (( ${#R1s[@]} >= 1 && ${#R2s[@]} >= 1 )); then
  paired=1; R1="${R1s[0]}"; R2="${R2s[0]}"
elif (( ${#SEs[@]} >= 1 )); then
  single=1; SE="${SEs[0]}"
else
  echo "[mapi] No FASTQ files found in $SAMPLEDIR"; exit 1
fi

if [[ -z "$PREFIX" ]]; then
  base="$(basename "$SAMPLEDIR")"
  PREFIX="${base%/}"
fi

OUTROOT="$(pwd)/${PREFIX}_MAPI"
mkdir -p "$OUTROOT"/{logs,qc_pre,filter_fastp,qc_post,align_bwa,lumpy}
export THREADS

log(){ echo "[$(date '+%F %T')] $*"; }

# 1) pre-QC
log "pre-QC (fastqc)"
if (( paired )); then
  mapi fastqc -o "$OUTROOT/qc_pre" "$R1" "$R2" > "$OUTROOT/logs/01_fastqc_pre.log" 2>&1
else
  mapi fastqc -o "$OUTROOT/qc_pre" "$SE"        > "$OUTROOT/logs/01_fastqc_pre.log" 2>&1
fi

# 2) filtering
log "filtering (fastp)"
if (( paired )); then
  mapi fastp -o "$OUTROOT/filter_fastp" -p "$PREFIX" -1 "$R1" -2 "$R2" > "$OUTROOT/logs/02_fastp.log" 2>&1
  F1="$OUTROOT/filter_fastp/${PREFIX}_R1.filtered.fastq.gz"
  F2="$OUTROOT/filter_fastp/${PREFIX}_R2.filtered.fastq.gz"
else
  mapi fastp -o "$OUTROOT/filter_fastp" -p "$PREFIX" -s "$SE" > "$OUTROOT/logs/02_fastp.log" 2>&1
  FSE="$OUTROOT/filter_fastp/${PREFIX}.filtered.fastq.gz"
fi

# 3) post-QC
log "post-QC (fastqc)"
if (( paired )); then
  mapi fastqc -o "$OUTROOT/qc_post" "$F1" "$F2" > "$OUTROOT/logs/03_fastqc_post.log" 2>&1
else
  mapi fastqc -o "$OUTROOT/qc_post" "$FSE"      > "$OUTROOT/logs/03_fastqc_post.log" 2>&1
fi

# 4) alignment
log "alignment (bwa-mem2 + samtools)"
if (( paired )); then
  mapi bwa -r "$REF" -o "$OUTROOT/align_bwa" -p "$PREFIX" -1 "$F1" -2 "$F2" > "$OUTROOT/logs/04_bwa.log" 2>&1
else
  mapi bwa -r "$REF" -o "$OUTROOT/align_bwa" -p "$PREFIX" -s "$FSE"         > "$OUTROOT/logs/04_bwa.log" 2>&1
fi
BAM="$OUTROOT/align_bwa/${PREFIX}.sorted.bam"

# 5) LUMPY (stub)
log "lumpy (placeholder)"
mapi lumpy -i "$BAM" -o "$OUTROOT/lumpy" -p "$PREFIX" > "$OUTROOT/logs/05_lumpy.log" 2>&1 || true

log "DONE. Outputs in: $OUTROOT"
