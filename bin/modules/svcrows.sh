#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# svcrows.sh (MAPI module)
# ============================================================

# -----------------------------
# Logging helpers
# -----------------------------
say() { echo "[svcrows] $*"; }
die() { say "ERROR: $*"; exit 1; }

# -----------------------------
# Usage
# -----------------------------
usage() {
  cat <<'USAGE'
Usage:
  mapi modules svcrows <command> [args...]

Commands (R-backed):
  scavenge                Run SVCROWS Scavenge mode
  hunt                    Run SVCROWS Hunt mode
  summarize               Summarize all FCL output in an output directory
  consensus-to-query      Convert consensus output to clean query list
  concat-samples          Concatenate sample files into one file

Commands (format conversion helpers):
  append-crows            Append / merge SVCROWS tables (helper)

Notes:
  • R-backed commands run inside the Conda env: $MAPI_ROOT/envs/svcrows
  • Most conversion commands are pure bash/awk and do not require R.

Examples:
  mapi modules svcrows scavenge --in ./SVCROWSin --out ./SVCROWSout --expand TRUE --bpfactor TRUE --defaultsizes FALSE --xs 5000 --xl 25000 --y1s 500 --y1l 2500 --y2s 50 --y2l 80
  mapi modules svcrows hunt --in ./SVCROWSin --features ./SVCROWSin/Featurelist.tsv --out ./SVCROWSout
  mapi modules svcrows summarize ./SVCROWSout ./SVCROWSout/summary.tsv
USAGE
}

# -----------------------------
# Resolve MAPI root
# -----------------------------
# Try to locate MAPI_ROOT from environment; fall back to script location.
MAPI_ROOT="${MAPI_ROOT:-}" 
if [[ -z "$MAPI_ROOT" ]]; then
  # wrapper lives under $MAPI_ROOT/bin/modules
  this="$(readlink -f "${BASH_SOURCE[0]}")"
  MAPI_ROOT="$(dirname "$(dirname "$(dirname "$this")")")"
fi

HELPER_DIR="$MAPI_ROOT/bin/modules/.svcrows"
[[ -d "$HELPER_DIR" ]] || die "Missing helper dir: $HELPER_DIR"

ENV_NAME="svcrows"
ENV_PATH="$MAPI_ROOT/envs/$ENV_NAME"

run_r() {
  [[ -d "$ENV_PATH" ]] || die "Missing conda env at: $ENV_PATH (build it from tools/yaml/svcrows.yml)"
  # Ensure a predictable tempdir (important on some HPCs / shared filesystems)
  local tmp="${TMPDIR:-}" 
  if [[ -z "$tmp" ]]; then
    tmp="${MAPI_SCRATCH:-/tmp}/tmp.svcrows.$$"
    mkdir -p "$tmp" || true
  fi
  TMPDIR="$tmp" conda run --no-capture-output -p "$ENV_PATH" Rscript "$HELPER_DIR/svcrows_cli.R" "$@"
}

run_sh() {
  local script="$1"; shift
  [[ -x "$HELPER_DIR/$script" ]] || die "Missing helper script: $HELPER_DIR/$script"
  "$HELPER_DIR/$script" "$@"
}

# -----------------------------
# Main
# -----------------------------
if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

cmd="$1"; shift || true

case "$cmd" in
  scavenge|hunt|summarize|consensus-to-query|concat-samples)
    run_r "$cmd" "$@"
    ;;
  append-crows)
    run_sh "append_crows.sh" "$@"
    ;;

  *)
    say "Unknown command: $cmd"
    usage
    exit 2
    ;;
esac
