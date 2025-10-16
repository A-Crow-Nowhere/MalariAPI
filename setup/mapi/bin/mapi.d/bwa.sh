#!/usr/bin/env bash
set -euo pipefail
. "$HOME/tools/mapi/lib.sh"

ENV_NAME="bwa-mem2-env"
YAML_BASE="bwa-mem2"
THREADS="${THREADS:-8}"

usage() {
  cat <<EOF
Usage:
  Paired: mapi bwa -r <ref.fa> -o <outdir> -p <prefix> -1 R1.fq.gz -2 R2.fq.gz
  Single: mapi bwa -r <ref.fa> -o <outdir> -p <prefix> -s Reads.fq.gz
EOF
  exit 1
}

REF=""; OUTDIR=""; PREFIX=""; R1=""; R2=""; SE=""
while getopts ":r:o:p:1:2:s:" opt; do
  case "$opt" in
    r) REF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    s) SE="$OPTARG" ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

[[ -z "$REF" || -z "$OUTDIR" || -z "$PREFIX" ]] && usage
[[ -n "$SE" && (-n "$R1" || -n "$R2") ]] && usage
[[ -z "$SE" && ( -z "$R1" || -z "$R2" ) ]] && usage
[[ -f "$REF" ]] || { echo "[bwa] Reference not found: $REF"; exit 1; }

mkdir -p "$OUTDIR"
ensure_env "$ENV_NAME" "$YAML_BASE"

# auto-index if missing
need_idx=0
for ext in .0123 .amb .ann .bwt.2bit.64 .pac .sa; do
  [[ -f "${REF}${ext}" ]] || { need_idx=1; break; }
done
(( need_idx )) && run_in_env "$ENV_NAME" bwa-mem2 index "$REF"

BAM="$OUTDIR/${PREFIX}.sorted.bam"
if [[ -n "$SE" ]]; then
  run_in_env "$ENV_NAME" bash -lc "bwa-mem2 mem -t $THREADS '$REF' '$SE' | samtools sort -@ $THREADS -o '$BAM' -"
else
  run_in_env "$ENV_NAME" bash -lc "bwa-mem2 mem -t $THREADS '$REF' '$R1' '$R2' | samtools sort -@ $THREADS -o '$BAM' -"
fi
run_in_env "$ENV_NAME" samtools index -@ "$THREADS" "$BAM"
echo "[bwa] Output → $BAM"
