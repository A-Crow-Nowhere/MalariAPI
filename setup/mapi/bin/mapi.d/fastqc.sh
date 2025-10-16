#!/usr/bin/env bash
set -euo pipefail
. "$HOME/tools/mapi/lib.sh"

ENV_NAME="fastqc-env"
YAML_BASE="fastqc"
THREADS="${THREADS:-4}"

usage(){ echo "Usage: mapi fastqc -o <outdir> <reads1.fastq[.gz]> [reads2.fastq[.gz] ...]"; exit 1; }

OUTDIR=""
while getopts ":o:" opt; do
  case "$opt" in
    o) OUTDIR="$OPTARG" ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))
[[ -z "${OUTDIR}" || "$#" -lt 1 ]] && usage

mkdir -p "$OUTDIR"
ensure_env "$ENV_NAME" "$YAML_BASE"
run_in_env "$ENV_NAME" fastqc -t "$THREADS" -o "$OUTDIR" "$@"
echo "[fastqc] Done → $OUTDIR"
