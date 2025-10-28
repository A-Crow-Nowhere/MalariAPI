#!/usr/bin/env bash
# MAPI one-shot installer
# - Pulls repo (git if available; else GitHub auto-tarball)
# - Ensures bin/mapi + all bin/**/*.sh are executable
# - Adds MalariAPI bin paths to PATH (bash/zsh)
# - Collects env YAMLs from modules/yaml + pipelines/yaml into ~/MalariAPI/envs/yaml
# - Creates/updates conda/mamba envs for each YAML
set -euo pipefail

# ======== CONFIG (override by exporting before running, e.g. BRANCH=v1.0) ========
REPO_OWNER="${REPO_OWNER:-A-Crow-Nowhere}"
REPO_NAME="${REPO_NAME:-MalariAPI}"
BRANCH="${BRANCH:-main}"
ROOT="${ROOT:-$HOME/MalariAPI}"
SYNC_DIRS_DEFAULT="bin modules packages pipelines scripts templates tools"
SYNC_DIRS=(${SYNC_DIRS:-$SYNC_DIRS_DEFAULT})
# ================================================================================

say() { printf '%s\n' "$*" >&2; }
have() { command -v "$1" >/dev/null 2>&1; }

# Safer copy/sync helper (rsync if present; otherwise hard replace)
sync_dir() {
  local src="$1" dst="$2"
  mkdir -p "$dst"
  if have rsync; then
    rsync -a --delete "$src"/ "$dst"/
  else
    rm -rf "$dst"/*
    cp -a "$src"/. "$dst"/
  fi
}

# Add PATH line idempotently
add_path_line() {
  local rc="$1"
  local line='export PATH="$HOME/MalariAPI/bin:$HOME/MalariAPI/bin/scripts:$PATH"'
  [[ -f "$rc" ]] || return 0
  grep -Fq 'MalariAPI/bin' "$rc" || printf '\n# MalariAPI executables\n%s\n' "$line" >> "$rc"
}

# --- 0) Prep target tree --------------------------------------------------------
say "==> Target install directory: $ROOT"
mkdir -p "$ROOT"

# --- 1) Acquire/refresh repo into $ROOT ----------------------------------------
if have git; then
  say "==> Using git to clone/pull ${REPO_OWNER}/${REPO_NAME}@${BRANCH}"
  if [[ -d "$ROOT/.git" ]]; then
    ( cd "$ROOT"
      git fetch origin "$BRANCH" || true
      git checkout "$BRANCH" 2>/dev/null || git checkout -b "$BRANCH"
      git pull --rebase --autostash origin "$BRANCH" || true
    )
  else
    rm -rf "$ROOT"
    # Try SSH first; fallback to HTTPS
    git clone --branch "$BRANCH" "git@github.com:${REPO_OWNER}/${REPO_NAME}.git" "$ROOT" 2>/dev/null \
    || git clone --branch "$BRANCH" "https://github.com/${REPO_OWNER}/${REPO_NAME}.git" "$ROOT"
  fi
else
  say "==> git not found; using GitHub auto-generated tarball"
  TARBALL_URL="https://codeload.github.com/${REPO_OWNER}/${REPO_NAME}/tar.gz/refs/heads/${BRANCH}"
  TMPDIR="$(mktemp -d)"; trap 'rm -rf "$TMPDIR"' EXIT
  curl -fsSL "$TARBALL_URL" -o "$TMPDIR/repo.tgz"
  tar -xzf "$TMPDIR/repo.tgz" -C "$TMPDIR"
  EXTRACTED_DIR="$(find "$TMPDIR" -maxdepth 1 -type d -name "${REPO_NAME}-*" | head -n1)"
  # Sync selected top-level folders
  for d in "${SYNC_DIRS[@]}"; do
    [[ -d "$EXTRACTED_DIR/$d" ]] || { say "    (skip) $d not present"; continue; }
    sync_dir "$EXTRACTED_DIR/$d" "$ROOT/$d"
  done
fi

# --- 2) Make executables usable everywhere -------------------------------------
if [[ -d "$ROOT/bin" ]]; then
  say "==> Marking executables in bin/ as +x"
  find "$ROOT/bin" -type f \( -name "mapi" -o -name "*.sh" \) -exec chmod +x {} \;
fi

# --- 3) Put MAPI on PATH (persist & current shell) ------------------------------
say "==> Ensuring MalariAPI bin paths are on your PATH"
add_path_line "$HOME/.bashrc"
add_path_line "$HOME/.zshrc"
export PATH="$HOME/MalariAPI/bin:$HOME/MalariAPI/bin/scripts:$PATH"

# --- 4) Collect env YAMLs into ~/MalariAPI/envs/yaml ---------------------------
say "==> Collecting environment YAMLs from modules/yaml and pipelines/yaml"
mkdir -p "$ROOT/envs/yaml"

copy_yaml_dir() {
  local src="$1"
  [[ -d "$src" ]] || return 0
  shopt -s nullglob
  for y in "$src"/*.yml "$src"/*.yaml; do
    cp -f "$y" "$ROOT/envs/yaml/"
    say "    copied env: $(basename "$y")"
  done
  shopt -u nullglob
}

copy_yaml_dir "$ROOT/modules/yaml"
copy_yaml_dir "$ROOT/pipelines/yaml"

# --- 5) (Optional) generate minimal env YAMLs via tools/gen_env.sh -------------
if [[ -x "$ROOT/tools/gen_env.sh" ]]; then
  say "==> Running tools/gen_env.sh to generate minimal env YAMLs (if any)"
  "$ROOT/tools/gen_env.sh" || true
  # If it writes to env/yaml, mirror to envs/yaml
  [[ -d "$ROOT/env/yaml" ]] && sync_dir "$ROOT/env/yaml" "$ROOT/envs/yaml"
fi

# --- 6) Create/Update Conda envs from envs/yaml/*.yml --------------------------
say "==> Creating/updating Conda envs from envs/yaml"
if have conda || have mamba; then
  CMD="conda"; have mamba && CMD="mamba"
  # ensure base conda is initialized (no-op if already)
  if have conda; then
    # shellcheck disable=SC1090
    source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
  fi
  shopt -s nullglob
  for y in "$ROOT/envs/yaml/"*.yml "$ROOT/envs/yaml/"*.yaml; do
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
  shopt -u nullglob
else
  cat <<'NOTE'
!! No conda/mamba detected in PATH.
   After installing Miniconda/Mambaforge, run this block to build envs:

for y in "$HOME/MalariAPI/envs/yaml/"*.yml "$HOME/MalariAPI/envs/yaml/"*.yaml; do
  [[ -e "$y" ]] || continue
  n="$(basename "${y%.*}")"
  if conda env list | awk '{print $1}' | grep -qx "$n"; then
    conda env update -n "$n" -f "$y" --prune
  else
    conda env create -n "$n" -f "$y"
  fi
done
NOTE
fi

# --- 7) Verify -----------------------------------------------------------------
say "==> Verifying commands on PATH"
command -v mapi >/dev/null && say "    OK: mapi at $(command -v mapi)" || say "    MISSING: mapi"
command -v bed_to_vcf.sh >/dev/null && say "    OK: bed_to_vcf.sh on PATH" || say "    MISSING: bed_to_vcf.sh"

say "==> Done. Open a new shell or 'source ~/.bashrc' / 'source ~/.zshrc' for PATH persistence."
