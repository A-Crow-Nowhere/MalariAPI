#!/usr/bin/env bash
set -euo pipefail
. "$HOME/tools/mapi/lib.sh"

ENV_NAME="lumpy-env"
YAML_BASE="lumpy"

usage(){ echo "Usage: mapi lumpy -i <sorted.bam> -o <outdir> -p <prefix>   [PLACEHOLDER]"; exit 1; }

BAM=""; OUTDIR=""; PREFIX=""
while getopts ":i:o:p:" opt; do
  case "$opt" in
    i) BAM="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

[[ -z "$BAM" || -z "$OUTDIR" || -z "$PREFIX" ]] && usage
mkdir -p "$OUTDIR"
ensure_env "$ENV_NAME" "$YAML_BASE"

echo "[lumpy] TODO: implement splitters/discordants extraction and lumpy command line."
touch "$OUTDIR/${PREFIX}.lumpy.STUB.txt"
