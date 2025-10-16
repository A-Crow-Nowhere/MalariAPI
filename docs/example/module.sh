
---

## `MalariAPI/docs/templates.md`

```markdown
# MAPI Templates

Use these templates to add new tools quickly and consistently.

---

## Module Template (`templates/yourtool.sh`)

```bash
#!/usr/bin/env bash
set -euo pipefail
. "$HOME/tools/mapi/lib.sh"

ENV_NAME="yourtool-env"        # must equal YAML 'name:'
YAML_BASE="yourtool"           # equals ~/envs/yaml/yourtool.yaml (basename)
THREADS="${THREADS:-4}"

usage(){
  cat <<'EOF'
Usage:
  mapi yourtool -i <input> -o <outdir> -p <prefix> [--flags]

Description:
  What this module does, what it expects, and what it writes.
Outputs:
  <outdir>/<prefix>.yourtool.out          # primary artifact (rename appropriately)
  <outdir>/_mapi_meta.tsv (optional)
  <outdir>/_mapi_done     (optional)
EOF
  exit 1
}

IN=""; OUTDIR=""; PREFIX=""
EXTRA=()

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

OUT_PRIMARY="$OUTDIR/${PREFIX}.yourtool.out"
META="$OUTDIR/_mapi_meta.tsv"

run_in_env "$ENV_NAME" yourtool --threads "$THREADS" -i "$IN" -o "$OUT_PRIMARY" "${EXTRA[@]:-}"

{
  echo -e "tool\tyourtool"
  echo -e "threads\t$THREADS"
  echo -e "input\t$IN"
  echo -e "output\t$OUT_PRIMARY"
} > "$META"

touch "$OUTDIR/_mapi_done"
echo "[yourtool] Output → $OUT_PRIMARY"
