#!/usr/bin/env bash
# MAPI one-shot installer (layout-aware + git-aware)
# - Clones/updates MalariAPI
# - Configures git user.name / user.email (optional)
# - Sets up origin (fork) + upstream (canonical) remotes
# - Ensures bin/mapi + bin/**/*.sh are executable
# - Adds MalariAPI bin paths to PATH (bash/zsh)
# - Ensures Miniconda exists and updates base env from envs/base.yml
#
# NOTE:
#   - Does NOT build all module/pipeline envs; that’s done later via MAPI tools.

set -euo pipefail

say(){ printf '%s\n' "$*" >&2; }
have(){ command -v "$1" >/dev/null 2>&1; }

# ======== DEFAULT CONFIG (overridable via env or flags) ========
REPO_OWNER="${REPO_OWNER:-A-Crow-Nowhere}"   # where we actually clone from
UPSTREAM_OWNER="${UPSTREAM_OWNER:-A-Crow-Nowhere}"  # canonical / “yours”
REPO_NAME="${REPO_NAME:-MalariAPI}"
BRANCH="${BRANCH:-main}"
ROOT="${ROOT:-$HOME/MalariAPI}"

# If empty, we install Miniconda into "$ROOT/tools/miniconda3"
MINICONDA_HOME="${MINICONDA_HOME:-}"

# Git identity (local to this repo)
GIT_USER_NAME="${GIT_USER_NAME:-}"
GIT_USER_EMAIL="${GIT_USER_EMAIL:-}"
# ==============================================================

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --branch, -b <name>          Git branch to use (default: $BRANCH)
  --root <path>                Install root (default: $ROOT)
  --repo-owner <name>          GitHub owner to clone from (default: $REPO_OWNER)
  --upstream-owner <name>      Canonical GitHub owner (default: $UPSTREAM_OWNER)
  --repo-name <name>           GitHub repo name (default: $REPO_NAME)
  --miniconda-home <path>      Existing Miniconda root to reuse
  --git-name <name>            git config user.name for this repo
  --git-email <email>          git config user.email for this repo
  -h, --help                   Show this help and exit

Typical contributor workflow:

  # If they forked your repo to <their-username>/MalariAPI:
  ./install_mapi.sh \\
    --repo-owner <their-username> \\
    --upstream-owner A-Crow-Nowhere \\
    --branch main \\
    --git-name "Their Name" \\
    --git-email "their_email@example.com"

EOF
}

# ----- Parse CLI args ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --branch|-b)
      BRANCH="$2"; shift 2 ;;
    --root)
      ROOT="$2"; shift 2 ;;
    --repo-owner)
      REPO_OWNER="$2"; shift 2 ;;
    --upstream-owner)
      UPSTREAM_OWNER="$2"; shift 2 ;;
    --repo-name)
      REPO_NAME="$2"; shift 2 ;;
    --miniconda-home)
      MINICONDA_HOME="$2"; shift 2 ;;
    --git-name)
      GIT_USER_NAME="$2"; shift 2 ;;
    --git-email)
      GIT_USER_EMAIL="$2"; shift 2 ;;
    -h|--help)
      usage
      exit 0 ;;
    *)
      say "Unknown option: $1"
      usage
      exit 1 ;;
  esac
done

if ! have git; then
  say "!! git is not installed or not on PATH."
  say "   Please install git and rerun this installer."
  exit 1
fi

# Normalize ROOT (portable, no weird heredoc side effects)
if ROOT="$(cd "$ROOT" 2>/dev/null && pwd)"; then
  :
else
  say "!! Could not resolve ROOT directory from '$ROOT'"
  exit 1
fi

# ----- Shell rc detection + PATH helpers --------------------------------------

detect_shell_rc() {
  local shell_rc=""

  # 1) Trust \$SHELL if present
  case "${SHELL:-}" in
    *bash) shell_rc="$HOME/.bashrc" ;;
    *zsh)  shell_rc="$HOME/.zshrc" ;;
  esac

  # 2) Fallback: inspect the running shell
  if [[ -z "$shell_rc" ]]; then
    local comm
    comm="$(ps -p $$ -o comm= 2>/dev/null || echo "")"
    case "$comm" in
      bash) shell_rc="$HOME/.bashrc" ;;
      zsh)  shell_rc="$HOME/.zshrc" ;;
    esac
  fi

  # 3) Last resort: default to bashrc (never a directory)
  [[ -n "$shell_rc" ]] || shell_rc="$HOME/.bashrc"

  printf '%s\n' "$shell_rc"
}

add_path_line(){
  local rc="$1"
  local line="export PATH=\"$ROOT/bin:$ROOT/bin/scripts:\$PATH\""

  # Only proceed if rc is a regular file (not dir, not empty)
  if [[ -z "${rc:-}" ]] || [[ ! -f "$rc" ]] || [[ -d "$rc" ]]; then
    return 0
  fi

  if ! grep -Fq "$ROOT/bin" "$rc"; then
    printf '\n# MalariAPI executables\n%s\n' "$line" >> "$rc"
  fi
}

add_conda_init_line() {
  local rc="$1"

  # rc must be a regular file
  if [[ -z "${rc:-}" ]] || [[ ! -f "$rc" ]] || [[ -d "$rc" ]]; then
    return 0
  fi

  [[ -n "$MINICONDA_HOME" ]] || return 0

  local init="source \"$MINICONDA_HOME/etc/profile.d/conda.sh\""
  if ! grep -Fq "$MINICONDA_HOME/etc/profile.d/conda.sh" "$rc"; then
    printf '\n# MalariAPI Miniconda init\n%s\n' "$init" >> "$rc"
  fi
}

configure_git_repo() {
  (
    cd "$ROOT"

    # Local identity (does not touch global config)
    if [[ -n "$GIT_USER_NAME" ]]; then
      say "==> Setting git user.name to '$GIT_USER_NAME' (local to this repo)"
      git config user.name "$GIT_USER_NAME"
    fi
    if [[ -n "$GIT_USER_EMAIL" ]]; then
      say "==> Setting git user.email to '$GIT_USER_EMAIL' (local to this repo)"
      git config user.email "$GIT_USER_EMAIL"
    fi

    # Set up upstream remote if this looks like a fork
    # origin → REPO_OWNER, upstream → UPSTREAM_OWNER
    if [[ "$REPO_OWNER" != "$UPSTREAM_OWNER" ]]; then
      local upstream_url
      upstream_url="git@github.com:${UPSTREAM_OWNER}/${REPO_NAME}.git"
      if ! git remote get-url upstream >/dev/null 2>&1; then
        say "==> Adding 'upstream' remote → $upstream_url"
        git remote add upstream "$upstream_url"
      fi
    fi

    say "==> Git remotes:"
    git remote -v || true
  )
}

ensure_miniconda() {
  # Decide where Miniconda lives if not provided
  if [[ -z "$MINICONDA_HOME" ]]; then
    MINICONDA_HOME="$ROOT/tools/miniconda3"
  fi

  if [[ -x "$MINICONDA_HOME/bin/conda" ]]; then
    say "==> Reusing existing Miniconda at $MINICONDA_HOME"
  else
    say "==> Installing Miniconda into $MINICONDA_HOME"
    mkdir -p "$(dirname "$MINICONDA_HOME")"
    local install_dir
    install_dir="$(dirname "$MINICONDA_HOME")"
    (
      cd "$install_dir"
      local INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
      local URL="https://repo.anaconda.com/miniconda/$INSTALLER"

      if have curl; then
        curl -fsSLo "$INSTALLER" "$URL"
      elif have wget; then
        wget -O "$INSTALLER" "$URL"
      else
        say "!! Neither curl nor wget is available; cannot download Miniconda."
        exit 1
      fi

      bash "$INSTALLER" -b -p "$MINICONDA_HOME"
      rm -f "$INSTALLER"
    )
  fi

  # Load conda into this shell
  if [[ -f "$MINICONDA_HOME/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "$MINICONDA_HOME/etc/profile.d/conda.sh"
  else
    say "!! Miniconda installed but etc/profile.d/conda.sh not found at $MINICONDA_HOME"
  fi

  # Persist conda init if Miniconda is under the MAPI tree
  case "$MINICONDA_HOME" in
    "$ROOT"/*)
      add_conda_init_line "$HOME/.bashrc"
      add_conda_init_line "$HOME/.zshrc"
      ;;
  esac
}

# --- 0) Prep tree --------------------------------------------------------------
say "==> Target install directory: $ROOT"
mkdir -p "$ROOT"

# --- 1) Acquire/refresh repo ---------------------------------------------------
if [[ -d "$ROOT/.git" ]]; then
  say "==> Existing git repo found; updating $BRANCH"
  (
    cd "$ROOT"
    git fetch origin "$BRANCH" || true
    git checkout "$BRANCH" 2>/dev/null || git checkout -b "$BRANCH"
    git sparse-checkout init --no-cone || true
    git sparse-checkout set '/*' '!:**/docs/**' '!:**/[Rr][Ee][Aa][Dd][Mm][Ee]*'
    git pull --rebase --autostash origin "$BRANCH" || true
    git sparse-checkout reapply || true
  )
else
  say "==> Cloning repo $REPO_OWNER/$REPO_NAME (branch: $BRANCH)"
  rm -rf "$ROOT"
  git clone --filter=blob:none --sparse --branch "$BRANCH" \
    "git@github.com:${REPO_OWNER}/${REPO_NAME}.git" "$ROOT" 2>/dev/null \
  || git clone --filter=blob:none --sparse --branch "$BRANCH" \
    "https://github.com/${REPO_OWNER}/${REPO_NAME}.git" "$ROOT"

  (
    cd "$ROOT"
    git sparse-checkout init --no-cone
    git sparse-checkout set '/*' '!:**/docs/**' '!:**/[Rr][Ee][Aa][Dd][Mm][Ee]*'
  )
fi

# --- 2) Git config / remotes ---------------------------------------------------
configure_git_repo

# --- 3) Executables ------------------------------------------------------------
if [[ -d "$ROOT/bin" ]]; then
  say "==> Marking executables in bin/ as +x"
  find "$ROOT/bin" -type f \( -name "mapi" -o -name "*.sh" \) -exec chmod +x {} \;
fi

# --- 4) PATH -------------------------------------------------------------------
say "==> Ensuring MalariAPI bin paths are on your PATH"

SHELL_RC="$(detect_shell_rc)"
add_path_line "$SHELL_RC"
# Also attempt to update both common rc files if they exist (harmless if not)
add_path_line "$HOME/.bashrc"
add_path_line "$HOME/.zshrc"

export PATH="$ROOT/bin:$ROOT/bin/scripts:$PATH"

# --- 5) Miniconda / base env ---------------------------------------------------
say "==> Ensuring Miniconda is available"
ensure_miniconda

if ! have conda; then
  say "!! conda command is still not available after Miniconda setup."
  say "   You may need to 'source ~/.bashrc' or restart your shell, then run:"
  say "     conda env update -n base -f \"$ROOT/envs/base.yml\"  (if present)"
else
  BASE_YAML="$ROOT/envs/base.yml"
  if [[ -f "$BASE_YAML" ]]; then
    say "==> Updating base environment from $BASE_YAML"
    # If you have CONDA_ROOT available, this is the most explicit:
    #   "$CONDA_ROOT/bin/conda" env update -p "$CONDA_ROOT" -f "$BASE_YAML" --prune
    # Using -n base is also fine:
    conda env update -n base -f "$BASE_YAML" --prune
  else
    say "==> No base.yml found at $BASE_YAML; skipping base env update."
  fi
fi
mv ~/MalariAPI/envs/base.yml ~/MalariAPI/envs/.base.yml


# Add MAPI envs to Conda search path
if ! grep -q 'CONDA_ENVS_DIRS' "$HOME/.bashrc"; then
  echo 'export CONDA_ENVS_DIRS="$HOME/MalariAPI/envs"' >> "$HOME/.bashrc"
  echo "[installer] Added CONDA_ENVS_DIRS to ~/.bashrc"
else
  echo "[installer] CONDA_ENVS_DIRS already configured in ~/.bashrc"
fi

bash ~/MalariAPI/tools/install_envs.sh ~/MalariAPI/envs/default.yml

conda activate default


# --- 6) Post-install guidance --------------------------------------------------
cat <<'NOTE'

chmod +x ~/MalariAPI/tools/git/update
chmod +x ~/MalariAPI/tools/git/upload
chmod +x ~/MalariAPI/tools/git/switch

==> Initial MAPI install complete.

Git:
  - Local user.name/user.email have been configured for this repo if provided.
  - If you cloned your own fork (REPO_OWNER != UPSTREAM_OWNER),
    an 'upstream' remote pointing at the canonical MalariAPI repo was added.

Conda:
  - Miniconda has been installed/reused.
  - The base environment was updated from envs/base.yml if present.

Next steps:
  - Open a new shell or run:
        source ~/.bashrc    # or ~/.zshrc
  - Create a feature branch and start hacking, e.g.:
        cd "$HOME/MalariAPI"
        git checkout -b feature/my_first_module

  - When ready, commit and push:
        git add -A
        git commit -m "Add my first MAPI module"
        git push -u origin feature/my_first_module

  - Then open a Pull Request on GitHub targeting the canonical main branch.

NOTE

say "==> Done. Remember to source your shell rc for PATH + conda persistence."
