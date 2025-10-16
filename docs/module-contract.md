# MAPI Module (Subcommand) Contract

This document defines the **required structure and behavior** for MAPI modules (the bash subcommands in `~/bin/mapi.d/`). The goal is to keep tools swappable and pipelines predictable.

## Location, Naming, Invocation

- One file per tool at `~/bin/mapi.d/<tool>.sh`
- Must be executable (`chmod +x`).
- Invoked as: `mapi <tool> [options]` (via the dispatcher at `~/bin/mapi`)

## Mandatory Structure (in order)

```bash
#!/usr/bin/env bash
set -euo pipefail
. "$HOME/tools/mapi/lib.sh"        # shared helpers

ENV_NAME="yourtool-env"            # MUST match YAML 'name:'
YAML_BASE="yourtool"               # MUST match ~/envs/yaml/yourtool.yaml
THREADS="${THREADS:-4}"            # honor global THREADS, default sensible

usage(){ cat <<'EOF'
Usage:
  mapi yourtool -i <input> -o <outdir> -p <prefix> [--flags...]

Description:
  Briefly describe inputs/outputs and what the module does.
EOF
exit 1; }

# ---- Parse args (required) ----
IN=""; OUTDIR=""; PREFIX=""
EXTRA=()                           # optional passthrough flags

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   IN="$2"; shift 2 ;;
    -o|--outdir)  OUTDIR="$2"; shift 2 ;;
    -p|--prefix)  PREFIX="$2"; shift 2 ;;
    --) shift; EXTRA+=("$@"); break ;;
    -h|--help)    usage ;;
    *)            EXTRA+=("$1"); shift ;;
  esac
done

[[ -n "$IN" && -n "$OUTDIR" && -n "$PREFIX" ]] || usage
mkdir -p "$OUTDIR"

ensure_env "$ENV_NAME" "$YAML_BASE"

# ---- Define outputs with standard names ----
OUT_PRIMARY="$OUTDIR/${PREFIX}.yourtool.out"   # change extension appropriately
META="$OUTDIR/_mapi_meta.tsv"

# ---- Run the tool ----
run_in_env "$ENV_NAME" yourtool --threads "$THREADS" -i "$IN" -o "$OUT_PRIMARY" "${EXTRA[@]:-}"

# ---- Optional sidecars ----
{
  echo -e "tool\tyourtool"
  echo -e "threads\t$THREADS"
  echo -e "input\t$IN"
  echo -e "output\t$OUT_PRIMARY"
} > "$META"

touch "$OUTDIR/_mapi_done"
echo "[yourtool] Output → $OUT_PRIMARY"
