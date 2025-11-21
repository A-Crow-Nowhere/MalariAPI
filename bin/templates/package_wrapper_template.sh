#!/usr/bin/env bash
# MAPI package wrapper
# Location: bin/packages/<name>
# Package files: bin/packages/<name>/...
set -euo pipefail

self="$(basename "$0")"
pkg="$self"  # 1:1 name = env = wrapper

usage(){ cat <<EOF
Usage: $self [args...]
Notes:
  • Package dir: bin/packages/$pkg/
  • Env prefix : envs/$pkg  (1:1 mapping)
  • YAML       : bin/packages/$pkg/yaml/$pkg.yml  (generated via tools/gen_env.sh)
  • Normally launched via: mapi packages $pkg -- [args...]
EOF
}

if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then usage; exit 0; fi

# Locate repo and package dir
SCRIPT_PATH="$(readlink -f "$0" 2>/dev/null || python3 - <<'PY'
import os,sys; print(os.path.realpath(sys.argv[1]))
PY
"$0")"
REPO_ROOT="$(cd "$(dirname "$SCRIPT_PATH")/../.." && pwd)"
PKG_DIR="$REPO_ROOT/bin/packages/$pkg"

[[ -d "$PKG_DIR" ]] || { echo "[error] package dir missing: $PKG_DIR" >&2; exit 4; }

# If invoked outside mapi, try to re-exec inside the env prefix
if [[ -z "${MAPI_REPO_ROOT:-}" ]]; then
  # Try to find conda
  resolve_conda_root() {
    for cand in \
      "${MAPI_CONDA_ROOT:-}" \
      "$REPO_ROOT/tools/miniconda3" \
      "$HOME/MalariAPI/tools/miniconda3" \
      "$HOME/tools/miniconda3"
    do
      [[ -n "$cand" && -x "$cand/bin/conda" ]] && { echo "$cand"; return 0; }
    done
    command -v conda >/dev/null 2>&1 && python3 - <<'PY'
import os, shutil
p = shutil.which("conda")
print(os.path.realpath(os.path.join(os.path.dirname(p), "..")) if p else "", end="")
PY
  }
  CONDA_ROOT="$(resolve_conda_root || true)"
  ENV_PREFIX="$REPO_ROOT/envs/$pkg"
  if [[ -n "${CONDA_ROOT:-}" && -x "$CONDA_ROOT/bin/conda" && -d "$ENV_PREFIX" ]]; then
    exec "$CONDA_ROOT/bin/conda" run -p "$ENV_PREFIX" -- "$SCRIPT_PATH" "$@"
  fi
fi

# -------------------- REAL PACKAGE ENTRY --------------------
# Choose ONE of these patterns and remove the others:

# A) Python package entrypoint (module)
if [[ -f "$PKG_DIR/__main__.py" ]]; then
  exec python3 -m "bin.packages.$pkg" "$@"   # requires repo-root on sys.path (adjust if needed)
fi

# B) Python script file
if [[ -f "$PKG_DIR/main.py" ]]; then
  exec python3 "$PKG_DIR/main.py" "$@"
fi

# C) Bash entrypoint
if [[ -x "$PKG_DIR/main.sh" ]]; then
  exec bash "$PKG_DIR/main.sh" "$@"
fi

echo "[error] No package entry found. Expected one of: __main__.py, main.py, main.sh in $PKG_DIR" >&2
exit 5
