#!/usr/bin/env bash
set -euo pipefail
. "$HOME/tools/mapi/lib.sh"

ENV_NAME="fastp-env"
YAML_BASE="fastp"
THREADS="${THREADS:-4}"

usage() {
  cat <<EOF
Usage:
  Paired: mapi fastp -o <outdir> -p <prefix> -1 R1.fq.gz -2 R2.fq.gz
  Single: mapi fastp -o <outdir> -p <prefix> -s Reads.fq.gz
EOF
  exit 1
}

OUTDIR=""; PREFIX=""; R1=""; R2=""; SE=""
while getopts ":o:p:1:2:s:" opt; do
  case "$opt" in
    o) OUTDIR="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    s) SE="$OPTARG" ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

[[ -z "$OUTDIR" || -z "$PREFIX" ]] && usage
[[ -n "$SE" && (-n "$R1" || -n "$R2") ]] && usage
[[ -z "$SE" && ( -z "$R1" || -z "$R2" ) ]] && usage

mkdir -p "$OUTDIR"
ensure_env "$ENV_NAME" "$YAML_BASE"

JSON="$OUTDIR/${PREFIX}.fastp.json"
HTML="$OUTDIR/${PREFIX}.fastp.html"

if [[ -n "$SE" ]]; then
  OUT="$OUTDIR/${PREFIX}.filtered.fastq.gz"
  run_in_env "$ENV_NAME" fastp -w "$THREADS" -i "$SE" -o "$OUT" -j "$JSON" -h "$HTML"
  echo "[fastp] Single-end → $OUT"
else
  OUT1="$OUTDIR/${PREFIX}_R1.filtered.fastq.gz"
  OUT2="$OUTDIR/${PREFIX}_R2.filtered.fastq.gz"
  run_in_env "$ENV_NAME" fastp -w "$THREADS" -i "$R1" -I "$R2" -o "$OUT1" -O "$OUT2" -j "$JSON" -h "$HTML"
  echo "[fastp] Paired-end → $OUT1 / $OUT2"
fi
