#!/usr/bin/env bash
set -euo pipefail

# ===========================
# MAPI installer (bash-first)
# ===========================
# Creates local dirs, fetches files from your repo, sets permissions, and
# enables the `mapi` dispatcher + subcommands.
#
# Usage:
#   bash install.sh                          # uses default BASEURL below
#   bash install.sh --baseurl URL            # override base (e.g., your fork/branch)
#   bash install.sh --dry-run                # show what would happen
#
# Re-run safe: idempotent.

# --- EDIT ME: point to your raw GitHub base once you publish ---
BASEURL_DEFAULT="https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/setup/mapi"
# Example if using a branch:
# BASEURL_DEFAULT="https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/my-branch/setup/mapi"

DRYRUN=0
BASEURL="$BASEURL_DEFAULT"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseurl) BASEURL="$2"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  endesac
done

# --- target locations (user may already have ~/bin) ---
BIN="$HOME/bin"
CMDDIR="$BIN/mapi.d"
ENVYAML="$HOME/envs/yaml"
TOOLS="$HOME/tools/mapi"

# --- helpers ---
say(){ printf '%s\n' "$*" ; }
doit(){ if (( DRYRUN )); then echo "DRYRUN: $*"; else eval "$@"; fi }
need_cmd(){ command -v "$1" >/dev/null 2>&1 || { say "ERROR: '$1' not found in PATH"; exit 1; }; }

# --- pick fetcher ---
if command -v curl >/dev/null 2>&1; then
  FETCH_CMD="curl -fsSL"
elif command -v wget >/dev/null 2>&1; then
  FETCH_CMD="wget -qO-"
else
  say "ERROR: need 'curl' or 'wget' to download files."
  exit 1
fi

fetch_to() {
  local url="$1" dest="$2"
  if (( DRYRUN )); then
    echo "DRYRUN: fetch $url -> $dest"
    return
  fi
  # shellcheck disable=SC2086
  $FETCH_CMD "$url" > "$dest"
}

# --- sanity: conda/mamba recommended ---
if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1 ; then
  say "WARN: conda/mamba not found on PATH. You can still install, but envs won't auto-create until conda/mamba is set up."
fi

say "==> Creating directories"
doit "mkdir -p '$BIN' '$CMDDIR' '$ENVYAML' '$TOOLS'"

# --- files to install: map of repo-path -> local-dest ---
# Adjust this list only if you change repo layout.
declare -A FILEMAP=(
  # dispatcher
  ["bin/mapi"]="$BIN/mapi"

  # subcommands
  ["bin/mapi.d/fastqc.sh"]="$CMDDIR/fastqc.sh"
  ["bin/mapi.d/fastp.sh"]="$CMDDIR/fastp.sh"
  ["bin/mapi.d/bwa.sh"]="$CMDDIR/bwa.sh"
  ["bin/mapi.d/lumpy.sh"]="$CMDDIR/lumpy.sh"
  ["bin/mapi.d/run.sh"]="$CMDDIR/run.sh"

  # shared lib
  ["tools/mapi/lib.sh"]="$TOOLS/lib.sh"

  # env yamls
  ["envs/yaml/fastqc.yaml"]="$ENVYAML/fastqc.yaml"
  ["envs/yaml/fastp.yaml"]="$ENVYAML/fastp.yaml"
  ["envs/yaml/bwa-mem2.yaml"]="$ENVYAML/bwa-mem2.yaml"
  ["envs/yaml/lumpy.yaml"]="$ENVYAML/lumpy.yaml"

  # docs (optional)
  ["README_MAPI.md"]="$TOOLS/README_MAPI.md"

  # templates (optional but handy)
  ["templates/yourtool.yaml"]="$ENVYAML/yourtool.yaml"
  ["templates/yourtool.sh"]="$CMDDIR/yourtool.sh"
)

say "==> Downloading files from: $BASEURL"
for rel in "${!FILEMAP[@]}"; do
  src="$BASEURL/$rel"
  dst="${FILEMAP[$rel]}"
  doit "mkdir -p \"$(dirname "$dst")\""
  fetch_to "$src" "$dst"
done

say "==> Setting execute bits"
for f in "$BIN/mapi" "$CMDDIR"/*.sh "$TOOLS/lib.sh"; do
  [[ -f "$f" ]] && doit "chmod +x \"$f\""
done

say "==> Creating convenience alias 'pipeline' → 'run'"
if [[ -f "$CMDDIR/run.sh" ]]; then
  doit "ln -sf \"$CMDDIR/run.sh\" \"$CMDDIR/pipeline.sh\""
fi

say "==> Install complete."

# quick usage hint
cat <<'TIP'

Try:
  mapi --help
  mapi fastqc -o ./qc your_R1.fastq.gz your_R2.fastq.gz
  mapi run --sample-dir ./Sample123 --ref /path/to/Pf3D7_v3.fa --threads 8

If 'mapi' isn't found:
  export PATH="$HOME/bin:$PATH"    # add to your shell rc to persist
TIP
