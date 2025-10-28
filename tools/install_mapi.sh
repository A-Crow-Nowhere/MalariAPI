#!/usr/bin/env bash
# MAPI one-shot installer (layout-aware)
# - Pulls repo (git if available; else GitHub auto-tarball)
# - Ensures bin/mapi + all bin/**/*.sh are executable
# - Adds MalariAPI bin paths to PATH (bash/zsh)
# - Builds conda/mamba envs from tools/yaml/*.yml and pipeline/yaml/*.yml
#   NOTE: envs/ holds actual runtime env dirs; we DO NOT copy YAMLs there.
set -euo pipefail

# ======== CONFIG (override via env vars before running) ========
REPO_OWNER="${REPO_OWNER:-A-Crow-Nowhere}"
REPO_NAME="${REPO_NAME:-MalariAPI}"
BRANCH="${BRANCH:-main}"
ROOT="${ROOT:-$HOME/MalariAPI}"

# Lock a specific Miniconda if you want (e.g., /home/you/tools/miniconda3)
MINICONDA_HOME="${MINICONDA_HOME:-}"
# ==============================================================

say(){ printf '%s\n' "$*" >&2; }
have(){ command -v "$1" >/dev/null 2>&1; }

# Safer dir sync helper (rsync if present)
sync_dir() {
  local src="$1" dst="$2"
  mkdir -p "$dst"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --delete "$src"/ "$dst"/
  else
    rm -rf "$dst"/*
    cp -a "$src"/. "$dst"/
  fi
}

add_path_line(){
  local rc="$1"
  local line='export PATH="$HOME/MalariAPI/bin:$HOME/MalariAPI/bin/scripts:$PATH"'
  [[ -f "$rc" ]] || return 0
  grep -Fq 'MalariAPI/bin' "$rc" || printf '\n# MalariAPI executables\n%s\n' "$line" >> "$rc"
}

# --- 0) Prep tree --------------------------------------------------------------
say "==> Target install directory: $ROOT"
mkdir -p "$ROOT"

# --- 1) Acquire/refresh repo ---------------------------------------------------
if [[ -d "$ROOT/.git" ]]; then
  ( cd "$ROOT"
    git fetch origin "$BRANCH" || true
    git checkout "$BRANCH" 2>/dev/null || git checkout -b "$BRANCH"
    # ensure sparse mode & patterns are applied on existing clones
    git sparse-checkout init --no-cone || true
    git sparse-checkout set '/*' '!:**/docs/**' '!:**/[Rr][Ee][Aa][Dd][Mm][Ee]*'
    git pull --rebase --autostash origin "$BRANCH" || true
    git sparse-checkout reapply || true
  )
else
  rm -rf "$ROOT"
  # sparse, partial clone (smaller & faster)
  git clone --filter=blob:none --sparse --branch "$BRANCH" \
    "git@github.com:${REPO_OWNER}/${REPO_NAME}.git" "$ROOT" 2>/dev/null \
  || git clone --filter=blob:none --sparse --branch "$BRANCH" \
    "https://github.com/${REPO_OWNER}/${REPO_NAME}.git" "$ROOT"

  ( cd "$ROOT"
    git sparse-checkout init --no-cone
    git sparse-checkout set '/*' '!:**/docs/**' '!:**/[Rr][Ee][Aa][Dd][Mm][Ee]*'
  )
fi

# --- 2) Executables ------------------------------------------------------------
if [[ -d "$ROOT/bin" ]]; then
  say "==> Marking executables in bin/ as +x"
  find "$ROOT/bin" -type f \( -name "mapi" -o -name "*.sh" \) -exec chmod +x {} \;
fi

# --- 3) PATH -------------------------------------------------------------------
say "==> Ensuring MalariAPI bin paths are on your PATH"
add_path_line "$HOME/.bashrc"
add_path_line "$HOME/.zshrc"
export PATH="$HOME/MalariAPI/bin:$HOME/MalariAPI/bin/scripts:$PATH"

# --- 4) Find YAML specs (NEW layout) -------------------------------------------
# YAMLs now live here:
YAML_DIRS=()
[[ -d "$ROOT/tools/yaml" ]]    && YAML_DIRS+=("$ROOT/tools/yaml")
[[ -d "$ROOT/pipeline/yaml" ]] && YAML_DIRS+=("$ROOT/pipeline/yaml")

if [[ "${#YAML_DIRS[@]}" -eq 0 ]]; then
  say "!! No YAML directories found (expected tools/yaml and/or pipeline/yaml). Skipping env build."
  exit 0
fi

# --- 5) Conda/Mamba boot -------------------------------------------------------
# Use locked Miniconda if provided; else whatever conda in PATH.
if [[ -n "$MINICONDA_HOME" ]]; then
  # shellcheck disable=SC1090
  source "$MINICONDA_HOME/etc/profile.d/conda.sh" 2>/dev/null || true
fi

# If conda is available, source its profile to ensure 'conda env' works in non-login shells
if have conda; then
  # shellcheck disable=SC1090
  source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
fi

if ! have conda && ! have mamba; then
  cat <<'NOTE'
!! No conda/mamba detected in PATH.
   After installing Miniconda/Mambaforge, run this block to build envs:

for y in "$HOME/MalariAPI/tools/yaml/"*.yml "$HOME/MalariAPI/tools/yaml/"*.yaml \
         "$HOME/MalariAPI/pipeline/yaml/"*.yml "$HOME/MalariAPI/pipeline/yaml/"*.yaml; do
  [[ -e "$y" ]] || continue
  n="$(basename "${y%.*}")"
  if conda env list | awk '{print $1}' | grep -qx "$n"; then
    conda env update -n "$n" -f "$y" --prune
  else
    conda env create -n "$n" -f "$y"
  fi
done
NOTE
  exit 0
fi

CMD="conda"; have mamba && CMD="mamba"

# --- 6) Create/Update envs directly from YAMLs ---------------------------------
say "==> Creating/updating Conda envs from tools/yaml and pipeline/yaml"
shopt -s nullglob
for dir in "${YAML_DIRS[@]}"; do
  for y in "$dir"/*.yml "$dir"/*.yaml; do
    [[ -e "$y" ]] || continue
    name="$(basename "${y%.*}")"
    if conda env list | awk '{print $1}' | grep -qx "$name"; then
      say "    updating env: $name"
      "$CMD" env update -n "$name" -f "$y" --prune
    else
      say "    creating env: $name"
      "$CMD" env create -n "$name" -f "$y"
    fi
  done
done
shopt -u nullglob

# --- 7) Verify -----------------------------------------------------------------
say "==> Verifying commands on PATH"
command -v mapi >/dev/null && say "    OK: mapi at $(command -v mapi)" || say "    MISSING: mapi"
command -v bed_to_vcf.sh >/dev/null && say "    OK: bed_to_vcf.sh on PATH" || say "    MISSING: bed_to_vcf.sh"

say "==> Done. Open a new shell or 'source ~/.bashrc' / 'source ~/.zshrc' for PATH persistence."
