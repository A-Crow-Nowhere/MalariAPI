#!/usr/bin/env bash
# MAPI one-shot installer (clean + idempotent)
# - Clones/updates MalariAPI
# - Ensures Miniconda in $ROOT/tools/miniconda3 (unless overridden)
# - Adds MAPI bin paths to PATH (bash/zsh)
# - Optionally adds conda init to rc (only if --global-conda)
# - Optionally builds default env via tools/install_envs.sh + tools/yaml/default.yml
#
# NON-GOALS (moved to install_hpc):
# - SSH/HPC host config
# - Remote rsync/miniconda/bootstrap
#
set -euo pipefail

say(){ printf '%s\n' "$*" >&2; }
die(){ say "ERROR: $*"; exit 1; }
have(){ command -v "$1" >/dev/null 2>&1; }

# ======== DEFAULT CONFIG (overridable via env or flags) ========
REPO_OWNER="${REPO_OWNER:-A-Crow-Nowhere}"          # where we clone from (fork/user)
UPSTREAM_OWNER="${UPSTREAM_OWNER:-A-Crow-Nowhere}"  # canonical upstream
REPO_NAME="${REPO_NAME:-MalariAPI}"
BRANCH="${BRANCH:-main}"
ROOT="${ROOT:-$HOME/MalariAPI}"

MINICONDA_HOME="${MINICONDA_HOME:-}"               # default: $ROOT/tools/miniconda3

# repo-local git identity (optional)
GIT_USER_NAME="${GIT_USER_NAME:-}"
GIT_USER_EMAIL="${GIT_USER_EMAIL:-}"

# env install behavior
INSTALL_DEFAULT_ENV="${INSTALL_DEFAULT_ENV:-1}"     # 1=yes, 0=no
DEFAULT_ENV_YAML="${DEFAULT_ENV_YAML:-}"            # default: auto-detect
ACTIVATE_DEFAULT_ENV="${ACTIVATE_DEFAULT_ENV:-0}"   # 1=yes, 0=no

# conda init behavior
GLOBAL_CONDA="${GLOBAL_CONDA:-0}"                   # 1=yes (persist conda init), 0=no

# git setup behavior (only if explicitly requested / args provided)
SETUP_GIT="${SETUP_GIT:-0}"                         # 1=yes, 0=no
# ===============================================================

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Core install options:
  --branch, -b <name>            Git branch to use (default: $BRANCH)
  --root <path>                  Install root (default: $ROOT)
  --repo-owner <name>            GitHub owner to clone from (default: $REPO_OWNER)
  --repo-name <name>             GitHub repo name (default: $REPO_NAME)
  --miniconda-home <path>        Existing Miniconda root to reuse
                                (default: \$ROOT/tools/miniconda3)

Git setup (ONLY runs if you pass these, or --setup-git):
  --git-name <name>              git config user.name for this repo (local only)
  --git-email <email>            git config user.email for this repo (local only)
  --upstream-owner <name>        Add upstream remote (git@github.com:<owner>/<repo>.git)
                                (default: $UPSTREAM_OWNER)
  --setup-git                    Force git setup even if no --git-name/--git-email provided

Environment build:
  --no-default-env               Skip building default env (tools/install_envs.sh)
  --default-env-yaml <path>      YAML spec for default env (default: auto-detect)
  --activate-default             After install, activate env 'default' in THIS shell

Conda init persistence:
  --global-conda                 Persist: source conda.sh + set CONDA_ENVS_DIRS in rc files
                                (default: off)

  -h, --help                     Show this help and exit

Examples:
  ./install_mapi.sh
  ./install_mapi.sh --repo-owner <yourfork> --upstream-owner A-Crow-Nowhere --setup-git
  ./install_mapi.sh --git-name "Your Name" --git-email you@domain.edu --setup-git
  ./install_mapi.sh --no-default-env
EOF
}

need_arg() {
  local opt="$1"
  [[ $# -ge 2 && -n "${2:-}" && "${2:-}" != --* ]] || die "$opt requires an argument"
}

# ----- Parse CLI args ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --branch|-b) need_arg "$1" "${2:-}"; BRANCH="$2"; shift 2;;
    --root) need_arg "$1" "${2:-}"; ROOT="$2"; shift 2;;
    --repo-owner) need_arg "$1" "${2:-}"; REPO_OWNER="$2"; shift 2;;
    --upstream-owner) need_arg "$1" "${2:-}"; UPSTREAM_OWNER="$2"; shift 2;;
    --repo-name) need_arg "$1" "${2:-}"; REPO_NAME="$2"; shift 2;;
    --miniconda-home) need_arg "$1" "${2:-}"; MINICONDA_HOME="$2"; shift 2;;

    --git-name) need_arg "$1" "${2:-}"; GIT_USER_NAME="$2"; SETUP_GIT=1; shift 2;;
    --git-email) need_arg "$1" "${2:-}"; GIT_USER_EMAIL="$2"; SETUP_GIT=1; shift 2;;
    --setup-git) SETUP_GIT=1; shift 1;;

    --no-default-env) INSTALL_DEFAULT_ENV=0; shift 1;;
    --default-env-yaml) need_arg "$1" "${2:-}"; DEFAULT_ENV_YAML="$2"; shift 2;;
    --activate-default) ACTIVATE_DEFAULT_ENV=1; shift 1;;

    --global-conda) GLOBAL_CONDA=1; shift 1;;

    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# ----- Basic prereqs ----------------------------------------------------------
have git || die "git is not installed. (sudo apt install -y git)"
have curl || have wget || die "Need curl or wget to download Miniconda."

# ----- Resolve ROOT safely ----------------------------------------------------
mkdir -p "$ROOT"
ROOT="$(cd "$ROOT" && pwd)"

# ----- Shell rc detection -----------------------------------------------------
detect_shell_rc() {
  local shell_rc=""
  case "${SHELL:-}" in
    *bash) shell_rc="$HOME/.bashrc" ;;
    *zsh)  shell_rc="$HOME/.zshrc" ;;
  esac
  if [[ -z "$shell_rc" ]]; then
    local comm=""
    comm="$(ps -p $$ -o comm= 2>/dev/null || true)"
    case "$comm" in
      bash) shell_rc="$HOME/.bashrc" ;;
      zsh)  shell_rc="$HOME/.zshrc" ;;
    esac
  fi
  [[ -n "$shell_rc" ]] || shell_rc="$HOME/.bashrc"
  printf '%s\n' "$shell_rc"
}

ensure_rc_file() {
  local rc="$1"
  [[ -n "$rc" ]] || return 0
  [[ -d "$rc" ]] && return 0
  [[ -f "$rc" ]] || : > "$rc"
}

append_once() {
  local rc="$1"
  local marker="$2"
  local line="$3"
  ensure_rc_file "$rc"
  grep -Fq "$marker" "$rc" 2>/dev/null || {
    printf '\n%s\n%s\n' "$marker" "$line" >> "$rc"
  }
}

add_mapi_path() {
  local rc="$1"
  local marker="# MalariAPI executables"
  local line="export PATH=\"$ROOT/bin:$ROOT/bin/scripts:\$PATH\""
  append_once "$rc" "$marker" "$line"
}

add_conda_init() {
  local rc="$1"
  local marker="# MalariAPI Miniconda init"
  local line="source \"$MINICONDA_HOME/etc/profile.d/conda.sh\""
  append_once "$rc" "$marker" "$line"
}

add_conda_envs_dirs() {
  local rc="$1"
  local marker="# MalariAPI conda envs dir"
  local line="export CONDA_ENVS_DIRS=\"$ROOT/envs\""
  append_once "$rc" "$marker" "$line"
}

# ----- Git repo config (ONLY if requested) -----------------------------------
configure_git_repo() {
  [[ "$SETUP_GIT" -eq 1 ]] || { say "==> Git setup skipped (no --git-* / --setup-git)."; return 0; }

  (
    cd "$ROOT"

    if [[ -n "$GIT_USER_NAME" ]]; then
      say "==> Setting git user.name (local): $GIT_USER_NAME"
      git config user.name "$GIT_USER_NAME"
    fi
    if [[ -n "$GIT_USER_EMAIL" ]]; then
      say "==> Setting git user.email (local): $GIT_USER_EMAIL"
      git config user.email "$GIT_USER_EMAIL"
    fi

    # Add upstream remote if owner differs OR user explicitly asked for git setup
    if [[ "$REPO_OWNER" != "$UPSTREAM_OWNER" || -n "$UPSTREAM_OWNER" ]]; then
      local upstream_url="git@github.com:${UPSTREAM_OWNER}/${REPO_NAME}.git"
      if ! git remote get-url upstream >/dev/null 2>&1; then
        say "==> Adding 'upstream' remote -> $upstream_url"
        git remote add upstream "$upstream_url" || true
      fi
    fi

    say "==> Git remotes:"
    git remote -v || true
  )
}

# ----- Miniconda --------------------------------------------------------------
ensure_miniconda() {
  if [[ -z "$MINICONDA_HOME" ]]; then
    MINICONDA_HOME="$ROOT/tools/miniconda3"
  fi

  if [[ -x "$MINICONDA_HOME/bin/conda" ]]; then
    say "==> Reusing Miniconda at $MINICONDA_HOME"
  else
    say "==> Installing Miniconda into $MINICONDA_HOME"
    mkdir -p "$(dirname "$MINICONDA_HOME")"
    (
      cd "$(dirname "$MINICONDA_HOME")"
      local INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
      local URL="https://repo.anaconda.com/miniconda/$INSTALLER"
      if have curl; then
        curl -fsSLo "$INSTALLER" "$URL"
      else
        wget -O "$INSTALLER" "$URL"
      fi
      bash "$INSTALLER" -b -p "$MINICONDA_HOME"
      rm -f "$INSTALLER"
    )
  fi

  [[ -f "$MINICONDA_HOME/etc/profile.d/conda.sh" ]] || die "conda.sh not found under $MINICONDA_HOME"
  # shellcheck disable=SC1091
  source "$MINICONDA_HOME/etc/profile.d/conda.sh"

  "$MINICONDA_HOME/bin/conda" config --set auto_activate_base false >/dev/null 2>&1 || true
}

# ----- Env YAML auto-detect ---------------------------------------------------
detect_default_env_yaml() {
  if [[ -f "$ROOT/tools/yaml/default.yml" ]]; then
    printf '%s\n' "$ROOT/tools/yaml/default.yml"; return 0
  fi
  # legacy fallbacks (harmless)
  if [[ -f "$ROOT/envs/default.yml" ]]; then
    printf '%s\n' "$ROOT/envs/default.yml"; return 0
  fi
  if [[ -f "$ROOT/envs/yaml/default.yml" ]]; then
    printf '%s\n' "$ROOT/envs/yaml/default.yml"; return 0
  fi
  return 1
}

install_default_env() {
  [[ "$INSTALL_DEFAULT_ENV" -eq 1 ]] || { say "==> Skipping default env install (--no-default-env)"; return 0; }

  local yml="$DEFAULT_ENV_YAML"
  if [[ -z "$yml" ]]; then
    yml="$(detect_default_env_yaml || true)"
  fi
  [[ -n "$yml" ]] || { say "==> No default env YAML found; skipping env build."; return 0; }

  [[ -x "$ROOT/tools/install_envs.sh" ]] || {
    say "==> tools/install_envs.sh missing or not executable; skipping env build."
    return 0
  }

  say "==> Building env(s) from: $yml"
  bash "$ROOT/tools/install_envs.sh" "$yml"

  if [[ "$ACTIVATE_DEFAULT_ENV" -eq 1 ]]; then
    say "==> Activating env 'default' in current shell"
    conda activate default || {
      say "!! Could not activate 'default'. It may have a different name than 'default'."
    }
  fi
}

# ----- 0) Prep ----------------------------------------------------------------
say "==> Target install directory: $ROOT"

# ----- 1) Acquire / refresh repo ---------------------------------------------
if [[ -d "$ROOT/.git" ]]; then
  say "==> Existing repo found; updating branch '$BRANCH'"
  (
    cd "$ROOT"
    git fetch origin "$BRANCH" || true
    git checkout "$BRANCH" 2>/dev/null || git checkout -b "$BRANCH"
    git sparse-checkout init --no-cone >/dev/null 2>&1 || true
    git sparse-checkout set '/*' '!:**/docs/**' >/dev/null 2>&1 || true
    git pull --rebase --autostash origin "$BRANCH" || true
    git sparse-checkout reapply >/dev/null 2>&1 || true
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
    git sparse-checkout set '/*' '!:**/docs/**'
  )
fi

# ----- 2) Git remotes / identity (optional) ----------------------------------
configure_git_repo

# ----- 3) Executables ---------------------------------------------------------
if [[ -d "$ROOT/bin" ]]; then
  say "==> Marking executables in bin/ as +x"
  find "$ROOT/bin" -type f \( -name "mapi" -o -name "*.sh" \) -exec chmod +x {} \;
fi
for f in "$ROOT/tools/git/update" "$ROOT/tools/git/upload" "$ROOT/tools/git/switch"; do
  [[ -f "$f" ]] && chmod +x "$f" || true
done

# ----- 4) PATH + conda init persistence --------------------------------------
say "==> Ensuring MalariAPI PATH is in your shell rc"
SHELL_RC="$(detect_shell_rc)"

for rc in "$HOME/.bashrc" "$HOME/.zshrc" "$SHELL_RC"; do
  add_mapi_path "$rc"
done

# ----- 5) Miniconda -----------------------------------------------------------
say "==> Ensuring Miniconda is available"
ensure_miniconda

# Make available immediately in THIS shell too
export PATH="$ROOT/bin:$ROOT/bin/scripts:$PATH"
export CONDA_ENVS_DIRS="$ROOT/envs"

# Persist conda init only if user opted in
if [[ "$GLOBAL_CONDA" -eq 1 ]]; then
  for rc in "$HOME/.bashrc" "$HOME/.zshrc" "$SHELL_RC"; do
    add_conda_init "$rc"
    add_conda_envs_dirs "$rc"
  done
  say "==> Global conda enabled: MAPI Miniconda will be available in new shells."
else
  say "==> Global conda NOT enabled: MAPI Miniconda will be used when running 'mapi' (and in this current shell)."
fi

# ----- 6) Build default env (optional) ---------------------------------------
install_default_env

# ----- 7) Post-install message ------------------------------------------------
cat <<NOTE

==> Initial MAPI install complete.

What was done:
  - Repo at:            $ROOT
  - PATH updated for:   $ROOT/bin and $ROOT/bin/scripts
  - Miniconda at:       ${MINICONDA_HOME:-$ROOT/tools/miniconda3}
  - CONDA_ENVS_DIRS set to: $ROOT/envs
  - Git setup:          $( [[ "$SETUP_GIT" -eq 1 ]] && echo "enabled" || echo "skipped" )
  - Conda init in rc:   $( [[ "$GLOBAL_CONDA" -eq 1 ]] && echo "enabled" || echo "skipped" )

Next:
  1) Restart your shell OR run:
       source ~/.bashrc    # (or ~/.zshrc)

  2) Verify:
       mapi --help
       conda --version

  3) If you skipped env creation and want it later:
       bash "$ROOT/tools/install_envs.sh" "$ROOT/tools/yaml/default.yml"

NOTE

say "==> Done."
