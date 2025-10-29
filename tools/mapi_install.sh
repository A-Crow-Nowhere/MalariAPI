THIS="$(readlink -f "$0" 2>/dev/null || python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$0")"
REPO_ROOT="$(cd "$(dirname "$THIS")/.." && pwd)"

resolve_conda_root() {
  for cand in \
    "${MAPI_CONDA_ROOT:-}" \
    "$REPO_ROOT/tools/miniconda3" \
    "$HOME/MalariAPI/tools/miniconda3" \
    "$HOME/tools/miniconda3"
  do
    [[ -n "$cand" && -x "$cand/bin/conda" ]] && { echo "$cand"; return 0; }
  done
  if command -v conda >/dev/null 2>&1; then
    python3 - <<'PY'
import os, shutil
p = shutil.which("conda")
print(os.path.realpath(os.path.join(os.path.dirname(p), "..")) if p else "", end="")
PY
    return 0
  fi
  return 1
}

CONDA_ROOT="$(resolve_conda_root)" || { echo "[error] conda not found (set MAPI_CONDA_ROOT)"; exit 1; }
CONDA_BIN="$CONDA_ROOT/bin/conda"
export CONDA_ENVS_DIRS="$REPO_ROOT/envs"
[[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]] && source "$CONDA_ROOT/etc/profile.d/conda.sh"
mkdir -p "$REPO_ROOT/envs"
bash "$REPO_ROOT/tools/install_envs.sh" "$@"
