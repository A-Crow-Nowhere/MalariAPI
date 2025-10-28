#!/usr/bin/env bash
# MAPI module template (validator-friendly)
# Recommended final location: bin/modules/<name>.sh
# Expected metadata: bin/modules/yaml/<name>.yml
set -euo pipefail

# Identify self early so Usage can print even if lock isn’t reachable yet
self="$(basename "$0")"
mod="${self%.sh}"

# ------------------- Usage (designed to satisfy validator) -------------------
usage(){ cat <<EOF
Usage: $self --in <reads.fastq> [--in2 <reads2.fastq>] --out <outdir> [--threads N] [-- ...tool flags...]
Aliases: --in1 ≡ --in ≡ -1 ; --in2 ≡ -2 ; --out ≡ --outdir ≡ -o ; --threads ≡ -t
Notes:
  • This is a template; real modules should live under bin/modules/<name>.sh
  • Metadata YAML is expected at bin/modules/yaml/<name>.yml unless overridden by \$MAPI_MODULE_YAML

EOF
}

# Let --help/-h succeed before sourcing anything (so validator can read Usage)
if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then
  usage; exit 0
fi

# ---------------- Enforce MAPI miniconda (validator-friendly) ----------------
# The FIRST line below MUST remain exactly as shown to satisfy E105.
# shellcheck disable=SC1090
source "$(dirname "$0")/../_mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../bin/_mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../._mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../bin/._mapi_conda_lock.sh" \
  || { echo "[ERROR] Could not source _mapi_conda_lock(.sh) from parent or bin/." >&2; exit 3; }

# ----------------------- Repo root & metadata discovery ----------------------
REPO_ROOT="${MAPI_REPO_ROOT:-"$(cd "$(dirname "$0")/../.." 2>/dev/null || pwd)"}"
meta="${MAPI_MODULE_YAML:-"$REPO_ROOT/bin/modules/yaml/$mod.yml"}"

if [[ ! -f "$meta" ]]; then
  echo "[warn] Missing metadata YAML (expected: $meta). Template will still run." >&2
fi

# Optionally activate env from YAML, if present
if [[ -f "$meta" ]]; then
  env_yaml="$(awk -F: '/^env:/ {f=1} f&&/conda:/ {gsub(/[\"'\'' ]/,""); print $2; exit}' "$meta" || true)"
  if [[ -n "${env_yaml:-}" && -f "$REPO_ROOT/$env_yaml" ]]; then
    env_name="$(basename "${env_yaml%.yml}")"
    conda env list >/dev/null 2>&1 && conda activate "$env_name" || true
  fi
fi

# ----------------------------- Defaults & CLI --------------------------------
THREADS="${THREADS:-8}"
OUTDIR="$(pwd)/out"
IN1="" IN2=""

EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in1|--in|-1)            IN1="$2"; shift 2;;
    --in2|-2)                 IN2="$2"; shift 2;;
    --out|--outdir|-o)        OUTDIR="$2"; shift 2;;
    --threads|-t)             THREADS="$2"; shift 2;;
    --)                       shift; EXTRA+=("$@"); break;;
    -h|--help)                usage; exit 0;;
    *) echo "Unknown: $1" >&2; usage; exit 2;;
  esac
done

# Validator-visible guards (documented in Usage above)
[[ -n "$IN1" ]]    || { echo "[error] missing --in/--in1" >&2; usage; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "[error] missing --out/--outdir" >&2; usage; exit 2; }
mkdir -p "$OUTDIR"

# ------------------------- IMPLEMENT YOUR TOOL HERE --------------------------
# Example (replace with your real command):
# tool_binary --threads "$THREADS" -i "$IN1" ${IN2:+-I "$IN2"} -o "$OUTDIR/result.txt" "${EXTRA[@]:-}"

# Safe no-op so the template runs:
echo "[$mod] meta=${meta:-none} env=${env_name:-none} threads=$THREADS in1=$IN1 in2=${IN2:-none} out=$OUTDIR" \
  | tee "$OUTDIR/${mod}.log"
echo "TEMPLATE OUTPUT for $mod" > "$OUTDIR/${mod}.out"
