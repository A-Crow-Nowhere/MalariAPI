#!/usr/bin/env bash
set -euo pipefail

# ===========================
# MAPI installer (bash-first)
# ===========================
# Creates local dirs, fetches files from your repo, sets permissions, and
# enables the `mapi` dispatcher + subcommands.

# --- EDIT ME: point to your raw GitHub base once you publish ---
BASEURL_DEFAULT="https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/setup/tools/mapi"

DRYRUN=0
BASEURL="$BASEURL_DEFAULT"

# Ensure we have bash ≥ 4 if we use associative arrays later
if ! command -v bash >/dev/null 2>&1; then
  echo "ERROR: bash not found"; exit 1
fi
# shellcheck disable=SC2001
BASH_MAJ="${BASH_VERSION%%.*}"
if [[ "${BASH_MAJ:-0}" -lt 4 ]]; then
  echo "WARN: Bash < 4 detected; proceeding with a portable path list (no assoc arrays)."
fi

# --- arg parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseurl) BASEURL="$2"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    -h|--help)
      cat <<EOF
Usage: bash install.sh [--baseurl URL] [--dry-run]
Default BASEURL:
  $BASEURL_DEFAULT
EOF
      exit 0
      ;;
    *) echo "Unknown arg: $1"; exit 2 ;;
  esac
done

# --- target locations ---
BIN="$HOME/bin"
CMDDIR="$BIN/mapi.d"
ENVYAML="$HOME/envs/yaml"
TOOLS="$HOME/tools/mapi"

say(){ printf '%s\n' "$*"; }
doit(){ if (( DRYRUN )); then echo "DRYRUN: $*"; else eval "$@"; fi }
need_cmd(){ command -v "$1" >/dev/null 2>&1 || { say "ERROR: '$1' not found in PATH"; exit 1; }; }

# --- pick fetcher ---
if command -v curl >/dev/null 2>&1; then
  FETCH() { curl -fsSL "$1"; }
elif command -v wget >/dev/null 2>&1; then
  FETCH() { wget -qO- "$1"; }
else
  say "ERROR: need 'curl' or 'wget' to download files."; exit 1
fi

fetch_to() {
  local url="$1" dest="$2"
  if (( DRYRUN )); then
    echo "DRYRUN: fetch $url -> $dest"
    return
  fi
  mkdir -p "$(dirname "$dest")"
  FETCH "$url" > "$dest"
}

# --- sanity: conda/mamba recommended ---
if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1 ; then
  say "WARN: conda/mamba not found on PATH. Envs will not auto-create until conda/mamba is set."
fi

say "==> Creating directories"
doit "mkdir -p '$BIN' '$CMDDIR' '$ENVYAML' '$TOOLS'"

# --- file list: repo relative path ::: local destination ---
# (string list avoids assoc arrays for Bash 3.x compatibility)
FILES=$(
  cat <<'EOF'
bin/mapi:::${BIN}/mapi
bin/mapi.d/fastqc.sh:::${CMDDIR}/fastqc.sh
bin/mapi.d/fastp.sh:::${CMDDIR}/fastp.sh
bin/mapi.d/bwa.sh:::${CMDDIR}/bwa.sh
bin/mapi.d/lumpy.sh:::${CMDDIR}/lumpy.sh
bin/mapi.d/run.sh:::${CMDDIR}/run.sh
tools/mapi/lib.sh:::${TOOLS}/lib.sh
envs/yaml/fastqc.yaml:::${ENVYAML}/fastqc.yaml
envs/yaml/fastp.yaml:::${ENVYAML}/fastp.yaml
envs/yaml/bwa-mem2.yaml:::${ENVYAML}/bwa-mem2.yaml
envs/yaml/lumpy.yaml:::${ENVYAML}/lumpy.yaml
templates/yourtool.yaml:::${ENVYAML}/yourtool.yaml
templates/yourtool.sh:::${CMDDIR}/yourtool.sh
EOF
)

say "==> Downloading files from: $BASEURL"
# shellcheck disable=SC2162
while read line; do
  [[ -z "$line" ]] && continue
  src_rel="${line%%:::*}"
  dst_tmpl="${line#*:::}"
  # expand ${VAR} in the destination path template
  eval "dst=\"$dst_tmpl\""
  src="$BASEURL/$src_rel"
  fetch_to "$src" "$dst"
done <<< "$FILES"

say "==> Setting execute bits"
for f in "$BIN/mapi" "$CMDDIR"/*.sh "$TOOLS/lib.sh"; do
  [[ -f "$f" ]] && doit "chmod +x \"$f\""
done

say "==> Creating convenience alias 'pipeline' → 'run'"
if [[ -f "$CMDDIR/run.sh" ]]; then
  doit "ln -sf \"$CMDDIR/run.sh\" \"$CMDDIR/pipeline.sh\""
fi

cat <<'TIP'

Install complete.

Try:
  mapi --help
  mapi doctor
  mapi fastqc -o ./qc your_R1.fastq.gz your_R2.fastq.gz
  mapi run --sample-dir ./Sample123 --ref /path/to/Pf3D7_v3.fa --threads 8

If 'mapi' isn't found:
  export PATH="$HOME/bin:$PATH"    # add to your shell rc to persist
TIP
