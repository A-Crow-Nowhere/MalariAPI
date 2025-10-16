#!/usr/bin/env bash
# Template module for MAPI
# Copy to:   ~/bin/mapi.d/yourtool.sh
# Requires:  ~/envs/yaml/yourtool.yaml  (with name: yourtool-env)
# Contract:  https://github.com/<YOUR_GH_USER>/MalariAPI/blob/main/docs/module-contract.md

set -euo pipefail
. "$HOME/tools/mapi/lib.sh"

# ---- module -> env mapping (must match your YAML) ----
ENV_NAME="yourtool-env"      # equals the 'name:' field in your YAML
YAML_BASE="yourtool"         # equals the YAML filename (basename) in ~/envs/yaml/
THREADS="${THREADS:-4}"      # honor global THREADS; set a sane default

usage(){
  cat <<'EOF'
Usage:
  mapi yourtool -i <input> -o <outdir> -p <prefix> [-- ...extra tool flags...]

Description:
  This is a starter template for a MAPI module. Fill in the actual tool binary
  and CLI flags where indicated below. Follow the standard I/O conventions:
    - Write primary outputs under <outdir>/, named with <prefix>.
    - Keep stdout minimal; the pipeline captures logs.
    - Do not cd; use absolute/quoted paths.
    - Respect THREADS if the tool supports threading.

Outputs (example; rename as appropriate):
  <outdir>/<prefix>.yourtool.out    # primary artifact
  <outdir>/_mapi_meta.tsv           # optional sidecar metadata
  <outdir>/_mapi_done               # optional success sentinel

Examples:
  mapi yourtool -i data.bam -o out -p SampleX -- --flag1 val --flag2
EOF
  exit 1
}

# ---- parse args (edit as needed) ----
IN=""
OUTDIR=""
PREFIX=""
EXTRA=()    # passthrough flags to underlying tool after '--'

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   IN="$2"; shift 2 ;;
    -o|--outdir)  OUTDIR="$2"; shift 2 ;;
    -p|--prefix)  PREFIX="$2"; shift 2 ;;
    --) shift; EXTRA+=("$@"); break ;;  # everything after -- is passed through
    -h|--help)    usage ;;
    *)            EXTRA+=("$1"); shift ;;
  esac
done

[[ -n "$IN" && -n "$OUTDIR" && -n "$PREFIX" ]] || usage
mkdir -p "$OUTDIR"

# ---- ensure environment exists (created from ~/envs/yaml/yourtool.yaml) ----
ensure_env "$ENV_NAME" "$YAML_BASE"

# ---- define outputs (rename extension to fit your tool) ----
OUT_PRIMARY="$OUTDIR/${PREFIX}.yourtool.out"
META="$OUTDIR/_mapi_meta.tsv"

# ---- main run (REPLACE 'yourtool' and flags with the real command) ----
# Example pattern:
# run_in_env "$ENV_NAME" yourtool --threads "$THREADS" -i "$IN" -o "$OUT_PRIMARY" "${EXTRA[@]:-}"

# Placeholder execution so the template is runnable without errors:
echo "[yourtool] (TEMPLATE) Would run with:"
echo "  env:    $ENV_NAME"
echo "  input:  $IN"
echo "  output: $OUT_PRIMARY"
echo "  flags:  ${EXTRA[*]:-(none)}" 

# create a stub file so pipeline steps don't break during dry tests
echo "TEMPLATE OUTPUT for yourtool on \$IN=$IN" > "$OUT_PRIMARY"

# ---- optional sidecars ----
{
  echo -e "tool\tyourtool (TEMPLATE)"
  echo -e "threads\t$THREADS"
  echo -e "input\t$IN"
  echo -e "output\t$OUT_PRIMARY"
} > "$META"
touch "$OUTDIR/_mapi_done"

echo "[yourtool] Output → $OUT_PRIMARY"
