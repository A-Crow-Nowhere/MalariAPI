#!/usr/bin/env bash
# Place copy at: bin/pipelines/<name>.sh
# Uses metadata:  bin/pipelines/yaml/<name>.yml
set -euo pipefail

source "$(dirname "$0")/../_mapi_conda_lock.sh" 2>/dev/null || true
REPO_ROOT="${MAPI_REPO_ROOT:-"$(cd "$(dirname "$0")/../.." && pwd)"}"
YAML_DIR="$REPO_ROOT/bin/pipelines/yaml"

self="$(basename "$0")"; pipe="${self%.sh}"
meta="${MAPI_PIPELINE_YAML:-$YAML_DIR/$pipe.yml}"
[[ -f "$meta" ]] || { echo "[pipeline:$pipe] Missing metadata: $meta" >&2; exit 3; }

# Optional: activate a wrapper env if specified in the pipeline YAML
env_yaml="$(awk -F: '/^env:/ {f=1} f&&/conda:/ {gsub(/[\"'\'' ]/,""); print $2; exit}' "$meta")"
if [[ -n "${env_yaml:-}" && -f "$REPO_ROOT/$env_yaml" ]]; then
  env_name="$(basename "${env_yaml%.yml}")"
  conda env list >/dev/null 2>&1 && conda activate "$env_name" || true
fi

SAMPLE="" OUTDIR="$(pwd)/run" THREADS="${THREADS:-8}"
usage(){ echo "Usage: $self --sample sample.yaml [--out run_dir] [--threads N]"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --out)    OUTDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown: $1" >&2; exit 2;;
  esac
done
[[ -f "$SAMPLE" ]] || { echo "Missing --sample"; exit 2; }
mkdir -p "$OUTDIR/logs"

# ---- CALL MODULES VIA mapi ----
# Example skeleton (replace inputs/refs as needed):
"$REPO_ROOT/bin/mapi" run modules fastp -- --in R1.fq.gz --in2 R2.fq.gz --out "$OUTDIR/fastp" --threads "$THREADS"
"$REPO_ROOT/bin/mapi" run modules bwa   -- --r1 R1.fq.gz --r2 R2.fq.gz --ref "$REPO_ROOT/genomes/pf3d7/pf3d7_GCF000002765v6.fa" --out "$OUTDIR/aln" --threads "$THREADS"

echo "[$pipe] done -> $OUTDIR"
