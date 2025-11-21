#!/usr/bin/env bash
# MAPI pipeline template (validator-friendly; mirrors module style)
# Recommended final location: bin/pipelines/<name>.sh
# Expected metadata: pipeline/yaml/<name>.yml
set -euo pipefail

# Identify self early so Usage can print even if lock isn’t reachable yet
self="$(basename "$0")"
pipe="${self%.sh}"

# ------------------- Usage (designed to satisfy validator) -------------------
usage(){ cat <<EOF
Usage: $self --prefix <run_prefix> [--threads N] [--in <file1>] [--in2 <file2>] [-- ...tool flags...]
Notes:
  • Envs: 'env_name' MUST equal pipeline name; env prefix = envs/<name>/
  • Outputs go under: {run_prefix}_mapi-out/  (per-pipeline run prefix)
  • Filenames SHOULD follow: {run_prefix}.<TAG>.<ext>  (default TAG shown below)
  • Default TAG for this pipeline: '\${DEFAULT_TAG}'
  • YAML: pipeline/yaml/<name>.yml with embedded env: block

Examples:
  $self --prefix cohortA --threads 16 --in samples.tsv
  $self --prefix test --in R1.fq.gz --in2 R2.fq.gz -- --extra-flag
EOF
}

# Hardcoded pipeline tag (documented in Usage above)
DEFAULT_TAG="${DEFAULT_TAG:-mapiPipe}"

# Let --help/-h succeed before sourcing anything
if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then usage; exit 0; fi

# ---------------- Enforce MAPI miniconda (compatible no-op if mapi already handled) ----------------
# shellcheck disable=SC1090
source "$(dirname "$0")/../_mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../bin/_mapi_conda_lock.sh" \
  || true

# ----------------------- Repo root & metadata discovery ----------------------
REPO_ROOT="${MAPI_REPO_ROOT:-"$(cd "$(dirname "$0")/../.." 2>/dev/null || pwd)"}"
meta="${MAPI_PIPELINE_YAML:-"$REPO_ROOT/pipeline/yaml/$pipe.yml"}"

if [[ ! -f "$meta" ]]; then
  echo "[warn] Missing pipeline YAML (expected: $meta). Template will still run." >&2
fi

# ----------------------------- Defaults & CLI --------------------------------
THREADS="${THREADS:-8}"
PREFIX=""
IN1="" IN2=""
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)                PREFIX="$2"; shift 2;;
    --threads|-t)            THREADS="$2"; shift 2;;
    --in|--in1|-1)           IN1="$2"; shift 2;;
    --in2|-2)                IN2="$2"; shift 2;;
    --)                      shift; EXTRA+=("$@"); break;;
    -h|--help)               usage; exit 0;;
    *) echo "Unknown: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$PREFIX" ]] || { echo "[error] missing --prefix" >&2; usage; exit 2; }

# Enforce outdir & filename conventions
OUTDIR="${PREFIX}_mapi-out"
mkdir -p "$OUTDIR"

# Helper to produce tagged filenames: ${OUTDIR}/${PREFIX}.${DEFAULT_TAG}.<ext>
out_path(){ printf "%s/%s.%s.%s\n" "$OUTDIR" "$PREFIX" "$DEFAULT_TAG" "$1"; }

# ------------------------- IMPLEMENT YOUR PIPELINE HERE ----------------------
# Orchestrate modules, etc. Example scaffold:
# mapi modules fastp      -- --in "$IN1" ${IN2:+--in2 "$IN2"} --threads "$THREADS"
# mapi modules bwa_mem2   -- --in "$IN1" ${IN2:+--in2 "$IN2"} --threads "$THREADS"
# ... then collect/emit a pipeline product:
echo "[$pipe] meta=${meta:-none} threads=$THREADS prefix=$PREFIX in1=${IN1:-none} in2=${IN2:-none} out=$OUTDIR tag=$DEFAULT_TAG" \
  | tee "$(out_path log)"
echo "TEMPLATE OUTPUT for pipeline $pipe" > "$(out_path out)"
