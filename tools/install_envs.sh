#!/usr/bin/env bash
# MAPI environment (re)installer - idempotent
# Builds/updates env prefixes under: <repo>/envs/<env_name>/
# Scans specs by default in:
#   - bin/modules/yaml/*.yml
#   - pipeline/yaml/*.yml
# OR, if one or more YAML paths are provided as positional arguments,
# processes *only* those paths, regardless of location.
set -Eeuo pipefail
[[ "${MAPI_DEBUG:-0}" == "1" ]] && set -x

# -------------------- Locate repo root --------------------
THIS="$(readlink -f "$0" 2>/dev/null || python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$0")"
REPO_ROOT="$(cd "$(dirname "$THIS")/.." && pwd)"

# -------------------- Conda auto-detect -------------------
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
[[ -x "$CONDA_BIN" ]] || { echo "[error] conda not executable at $CONDA_BIN" >&2; exit 1; }

# Detect whether we need / can use --environment-spec (newer Conda)
ENV_SPEC_ARGS=()
if "$CONDA_BIN" env create --help 2>&1 | grep -q -- '--environment-spec'; then
  ENV_SPEC_ARGS=(--environment-spec environment.yml)
fi

# Where we place built prefixes
ENV_PREFIX_ROOT="$REPO_ROOT/envs"
HASH_DIR="$ENV_PREFIX_ROOT/.mapi_envhash"
mkdir -p "$ENV_PREFIX_ROOT" "$HASH_DIR"

# Make sure conda looks here first
export CONDA_ENVS_DIRS="$ENV_PREFIX_ROOT"
# NOTE: we do NOT source conda.sh; conda env create/update is enough and faster.

# -------------------- Spec locations ----------------------
MODULE_SPECS_DIR="$REPO_ROOT/bin/modules/yaml"
PIPE_SPECS_DIR="$REPO_ROOT/pipeline/yaml"

# -------------------- CLI flags ---------------------------
FORCE=0
PRUNE=0
ONLY_LIST=""
SPEC_PATHS=()

usage() {
  cat <<EOF
Usage:
  tools/install_envs.sh [--force] [--prune] [--only env1,env2,...]
  tools/install_envs.sh [--force] [--prune] [--only env1,env2,...] SPEC.yml [SPEC2.yml ...]

Behavior:
  - With *no* SPEC paths:
      Build/update Conda envs from env specs in:
        - $MODULE_SPECS_DIR/*.yml
        - $PIPE_SPECS_DIR/*.yml
  - With one or more SPEC paths:
      Build/update envs *only* from those YAML files or directories,
      regardless of where they live on disk.

Supported spec formats:
  1) Plain Conda env YAML:
       name: crow
       channels: [...]
       dependencies: [...]
  2) Embedded MAPI-style spec:
       env_name: crow
       env:
         channels: [...]
         dependencies: [...]

Behavior details:
  - Idempotent via hash of each spec (+ env_name/name + conda major).
  - Creates/updates env prefix under: $ENV_PREFIX_ROOT/<env_name>/
  - Stores hashes under: $HASH_DIR/<env_name>.sha256

Options:
  --force           Rebuild all envs regardless of hash.
  --prune           Pass --prune to 'conda env update' (remove extras).
  --only a,b,c      Only process these env_name values (whether scanning or SPEC paths).
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE=1; shift;;
    --prune) PRUNE=1; shift;;
    --only)  ONLY_LIST="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        SPEC_PATHS+=("$1")
        shift
      done
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
    *)
      # Positional  treat as a spec path (file or directory)
      SPEC_PATHS+=("$1")
      shift
      ;;
  esac
done

# -------------------- YAML helpers ------------------------

# Detect kind of spec:
#   - "embedded" if it has env_name: and env:
#   - "plain"    if it has a top-level name:
detect_spec_kind() { # file
  local y="$1"
  if grep -qE '^env_name:' "$y" && grep -qE '^env:[[:space:]]*$' "$y"; then
    echo "embedded"
  elif grep -qE '^name:' "$y"; then
    echo "plain"
  else
    echo "unknown"
  fi
}

get_env_name_from_yaml() { # file
  local y="$1"
  # 1) embedded style: env_name:
  if grep -qE '^env_name:' "$y"; then
    grep -E '^env_name:' "$y" | head -n1 | sed 's/env_name:[[:space:]]*//'
  # 2) plain conda env: name:
  elif grep -qE '^name:' "$y"; then
    grep -E '^name:' "$y" | head -n1 | sed 's/name:[[:space:]]*//'
  else
    # fallback: basename without extension
    basename "$y" | sed 's/\.[Yy][Aa][Mm][Ll]$//; s/\.yml$//'
  fi
}

# Emit only the 'env:' block as a minimal env YAML (embedded style)
emit_env_block_embedded() { # spec out
  local spec="$1" out="$2"
  awk '
    BEGIN{inenv=0}
    /^env:[[:space:]]*$/ {inenv=1; print "name: placeholder"; next}
    inenv==1 && /^[^[:space:]]/ {inenv=0}
    inenv==1 {print}
  ' "$spec" | sed '1s/.*/name: placeholder/' > "$out"
}

# For plain env YAML, we just copy it as-is (but we don't rely on its name
# inside; we control the prefix path ourselves).
emit_env_block_plain() { # spec out
  local spec="$1" out="$2"
  cp "$spec" "$out"
}

# Compute a stable hash:
#   - Based on env block content + env_name + conda major + spec kind
spec_hash() { # spec env_name kind
  local spec="$1" env_name="$2" kind="$3"
  local tmp; tmp="$(mktemp --suffix=.yml)"
  if [[ "$kind" == "embedded" ]]; then
    emit_env_block_embedded "$spec" "$tmp"
  else
    emit_env_block_plain "$spec" "$tmp"
  fi
  local conda_major
  conda_major="$("$CONDA_BIN" --version | awk '{print $3}' | cut -d. -f1)"
  {
    echo "ENV_NAME=$env_name"
    echo "CONDA_MAJOR=$conda_major"
    echo "KIND=$kind"
    cat "$tmp"
  } | sha256sum | awk '{print $1}'
  rm -f "$tmp"
}

# -------------------- Core builder ------------------------
process_spec() { # spec file
  local spec="$1"
  [[ -f "$spec" ]] || return 0

  echo "[spec] $spec"

  local kind; kind="$(detect_spec_kind "$spec")"
  if [[ "$kind" == "unknown" ]]; then
    echo "[skip] Unknown spec format (no env_name/name + env): $spec"
    return 0
  fi

  local env_name; env_name="$(get_env_name_from_yaml "$spec")"
  if [[ -z "$env_name" ]]; then
    echo "[skip] Could not infer env_name for: $spec"
    return 0
  fi

  # --only filter
  if [[ -n "$ONLY_LIST" ]]; then
    IFS=',' read -r -a arr <<< "$ONLY_LIST"
    local keep=0
    for e in "${arr[@]}"; do [[ "$e" == "$env_name" ]] && keep=1; done
    [[ $keep -eq 1 ]] || { echo "[skip] $env_name (filtered by --only)"; return 0; }
  fi

  local prefix="$ENV_PREFIX_ROOT/$env_name"
  local new_hash; new_hash="$(spec_hash "$spec" "$env_name" "$kind")"
  local hash_file="$HASH_DIR/$env_name.sha256"
  local old_hash=""
  [[ -f "$hash_file" ]] && old_hash="$(cat "$hash_file" 2>/dev/null || true)"

  if [[ $FORCE -eq 0 && -n "$old_hash" && "$old_hash" == "$new_hash" && -d "$prefix" ]]; then
    echo "[skip] $env_name unchanged @ $prefix"
    return 0
  fi

  local tmp; tmp="$(mktemp --suffix=.yml)"
  if [[ "$kind" == "embedded" ]]; then
    emit_env_block_embedded "$spec" "$tmp"
  else
    emit_env_block_plain "$spec" "$tmp"
  fi

  if [[ ! -s "$tmp" ]]; then
    echo "[error] env block appeared empty for: $spec" >&2
    rm -f "$tmp"
    return 1
  fi

  if [[ -d "$prefix" ]]; then
    echo "[update] $env_name @ $prefix"
    if [[ $PRUNE -eq 1 ]]; then
      "$CONDA_BIN" env update "${ENV_SPEC_ARGS[@]}" -p "$prefix" -f "$tmp" --prune
    else
      "$CONDA_BIN" env update "${ENV_SPEC_ARGS[@]}" -p "$prefix" -f "$tmp"
    fi
  else
    echo "[create] $env_name @ $prefix"
    "$CONDA_BIN" env create "${ENV_SPEC_ARGS[@]}" -p "$prefix" -f "$tmp"
  fi

  echo "$new_hash" > "$hash_file"
  rm -f "$tmp"
}

# -------------------- Scan & build ------------------------
if [[ ${#SPEC_PATHS[@]} -gt 0 ]]; then
  # Explicit paths mode (files or dirs, anywhere)
  for path in "${SPEC_PATHS[@]}"; do
    if [[ -d "$path" ]]; then
      find "$path" -type f -name '*.yml' ! -name '.*' \
        | while read -r spec; do process_spec "$spec"; done
    else
      process_spec "$path"
    fi
  done
else
  # Default scanning mode: modules + pipelines in the repo
  if [[ -d "$MODULE_SPECS_DIR" ]]; then
    find "$MODULE_SPECS_DIR" -maxdepth 1 -type f -name '*.yml' ! -name '.*' \
      | while read -r spec; do process_spec "$spec"; done
  else
    echo "[info] No module env specs dir: $MODULE_SPECS_DIR"
  fi

  if [[ -d "$PIPE_SPECS_DIR" ]]; then
    find "$PIPE_SPECS_DIR" -maxdepth 1 -type f -name '*.yml' ! -name '.*' \
      | while read -r spec; do process_spec "$spec"; done
  else
    echo "[info] No pipeline env specs dir: $PIPE_SPECS_DIR"
  fi
fi

echo "Done. Envs live under: $ENV_PREFIX_ROOT"
