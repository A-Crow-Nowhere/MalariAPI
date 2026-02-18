#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# MAPI HPC installer (updated for current layout + "mapi summons conda" model)
#
# Key updates:
#  - NO remote ~/.bashrc conda init (conda should not be global; mapi summons it)
#  - Remote base env update reads YAML from tools/yaml/base.(yml|yaml) only
#  - Remote Miniconda download supports curl OR wget
#  - Excludes template uses __HPC_NAME__ placeholder (safer than <HPC_NAME>)
#  - Writes MAPI_CONDA_HOME into bin/.mapi.env (canonical conda location)
#  - Removes destructive rm -rf of tools/hpc_templates on remote
#  - Removes redundant tar exclude of tools/$HPC_NAME/submit
###############################################################################

###############################################################################
# 0) BASIC SETUP
###############################################################################

MAPI_ROOT="$HOME/MalariAPI"
CFG="$MAPI_ROOT/tools/hpc_config.json"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 3; }; }
need ssh; need rsync; need jq

mkdir -p "$MAPI_ROOT/tools"

###############################################################################
# 1) FIGURE OUT HPC_NAME, HOST, PARTITION, ALLOCATION
###############################################################################

HPC_NAME="${1:-}"
HOST="${2:-}"
PARTITION="${3:-}"
ALLOCATION="${4:-}"

if [[ -z "$HPC_NAME" || -z "$HOST" ]]; then
  if [[ -f "$CFG" ]]; then
    DEFAULT_HPC="$(jq -r '.default_hpc // empty' "$CFG")"
    [[ -n "$DEFAULT_HPC" ]] || {
      echo "No default_hpc in $CFG. Re-run as: $0 <hpc> <user@host> [partition] [allocation]" >&2
      exit 2
    }

    HOST_FROM_CFG="$(jq -r --arg h "$DEFAULT_HPC" '.clusters[$h].host // empty' "$CFG")"
    [[ -n "$HOST_FROM_CFG" ]] || {
      echo "Config missing host for '$DEFAULT_HPC'." >&2
      exit 2
    }

    HPC_NAME="$DEFAULT_HPC"
    HOST="$HOST_FROM_CFG"

    [[ -z "$PARTITION" ]] && PARTITION="$(jq -r --arg h "$HPC_NAME" '.clusters[$h].partition // empty' "$CFG")"
    [[ "$PARTITION" == "null" ]] && PARTITION=""

    [[ -z "$ALLOCATION" ]] && ALLOCATION="$(jq -r --arg h "$HPC_NAME" '.clusters[$h].allocation // empty' "$CFG")"
    [[ "$ALLOCATION" == "null" ]] && ALLOCATION=""

    echo "Using existing config ($CFG): HPC='$HPC_NAME', host='$HOST'"
  else
    echo "No config and no args. Usage: $0 <hpc> <user@host> [partition] [allocation]" >&2
    exit 2
  fi
fi

: "${PARTITION:=standard}"
: "${ALLOCATION:=}"

###############################################################################
# 2) RESOLVE REMOTE USER + REMOTE HOME (NO BatchMode YET)
###############################################################################

remote_user() {
  local host="$1"
  local u="${host%@*}"
  if [[ "$u" == "$host" || -z "$u" ]]; then
    u="$(ssh -o ConnectTimeout=8 "$host" 'printf "%s" "$USER"' 2>/dev/null || true)"
  fi
  printf "%s" "$u"
}

REMOTE_USER="$(remote_user "$HOST")"
REMOTE_HOME="$(ssh -o ConnectTimeout=8 "$HOST" 'printf "%s" "$HOME"' 2>/dev/null || true)"

[[ -n "$REMOTE_HOME" ]] || { echo "Could not determine remote HOME on $HOST" >&2; exit 4; }

REMOTE_MAPI_HOME="$REMOTE_HOME/MalariAPI"
REMOTE_SCRATCH="/scratch/$REMOTE_USER/MalariAPI"

SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=8)

###############################################################################
# 3) UPDATE CONFIG JSON + .mapi.env
###############################################################################

mkdir -p "$MAPI_ROOT/tools"

if [[ ! -f "$CFG" ]]; then
  cat >"$CFG" <<EOF
{"default_hpc":"$HPC_NAME","clusters":{"$HPC_NAME":{
  "host":"$HOST",
  "home_root":"$REMOTE_MAPI_HOME",
  "scratch_root":"$REMOTE_SCRATCH",
  "partition":"$PARTITION",
  "allocation":"$ALLOCATION"
}}}
EOF
else
  TMP="$(mktemp)"
  jq --arg h "$HPC_NAME" \
     --arg host "$HOST" \
     --arg rh "$REMOTE_MAPI_HOME" \
     --arg rs "$REMOTE_SCRATCH" \
     --arg part "$PARTITION" \
     --arg alloc "$ALLOCATION" '
    .default_hpc = ($h) |
    .clusters[$h] = {
      "host": $host,
      "home_root": $rh,
      "scratch_root": $rs,
      "partition": $part,
      "allocation": $alloc
    }
  ' "$CFG" > "$TMP"
  mv "$TMP" "$CFG"
fi

echo "Configured HPC '$HPC_NAME' ($HOST)"
echo "Remote HOME:    $REMOTE_MAPI_HOME"
echo "Remote SCRATCH: $REMOTE_SCRATCH"

###############################################################################
# 3b) MATERIALIZE .mapi.env
###############################################################################

mkdir -p "$MAPI_ROOT/bin"
MAPI_ENV_FILE="$MAPI_ROOT/bin/.mapi.env"
touch "$MAPI_ENV_FILE"

set_or_update() {
  local key="$1" val="$2"
  if grep -q "^export $key=" "$MAPI_ENV_FILE"; then
    sed -i "s|^export $key=.*$|export $key=\"$val\"|" "$MAPI_ENV_FILE"
  else
    echo "export $key=\"$val\"" >>"$MAPI_ENV_FILE"
  fi
}

add_if_missing() {
  local key="$1" val="$2"
  grep -q "^export $key=" "$MAPI_ENV_FILE" || echo "export $key=\"$val\"" >>"$MAPI_ENV_FILE"
}

add_if_missing "MAPI_ROOT"            "\$HOME/MalariAPI"
add_if_missing "MAPI_REMOTE_HOST"     "$HOST"
add_if_missing "REMOTE_MAPI_HOME"     "$REMOTE_MAPI_HOME"
add_if_missing "REMOTE_MAPI_SCRATCH"  "$REMOTE_SCRATCH"
add_if_missing "MAPI_SCRATCH"         "\$MAPI_ROOT/scratch"
add_if_missing "MAPI_SYNC_EXCLUDES_FILE" "$MAPI_ROOT/tools/.sync_excludes.txt"

# Canonical local conda location for MAPI (mapi summons this; no global shell changes)
add_if_missing "MAPI_CONDA_HOME"      "\$MAPI_ROOT/tools/miniconda3"

set_or_update "HPC_PARTITION" "$PARTITION"
set_or_update "HPC_ALLOCATION" "$ALLOCATION"

###############################################################################
# 3c) SAFE ~/.ssh/config ENTRY
###############################################################################

mkdir -p "$HOME/.ssh"
CONFIG_FILE="$HOME/.ssh/config"
[[ -f "$CONFIG_FILE" ]] || { touch "$CONFIG_FILE"; chmod 600 "$CONFIG_FILE"; }

BARE_HOST="${HOST#*@}"
if ! grep -qE "^Host[[:space:]]+$BARE_HOST\$" "$CONFIG_FILE"; then
  {
    echo ""
    echo "Host $BARE_HOST"
    echo "    HostName $BARE_HOST"
    [[ "$HOST" == *"@"* ]] && echo "    User ${HOST%@*}"
  } >> "$CONFIG_FILE"
  chmod 600 "$CONFIG_FILE"
fi

###############################################################################
# 4) ENSURE KEY-BASED LOGIN
###############################################################################

if ! ssh "${SSH_OPTS[@]}" "$HOST" true 2>/dev/null; then
  echo "No key-based login; installing key..."
  [[ -f "$HOME/.ssh/id_ed25519" ]] || ssh-keygen -t ed25519 -N "" -f "$HOME/.ssh/id_ed25519"
  ssh-copy-id "$HOST"
fi

ssh "${SSH_OPTS[@]}" "$HOST" true >/dev/null

###############################################################################
# 5) CREATE REMOTE DIRS
###############################################################################

ssh "${SSH_OPTS[@]}" "$HOST" \
  "mkdir -p \"$REMOTE_MAPI_HOME\" \"$REMOTE_SCRATCH/scratch\" /scratch/\$USER/mapi_logs" 2>/dev/null

###############################################################################
# 5a) SYNC EXCLUDES
###############################################################################

mkdir -p "$MAPI_ROOT/tools"
SYNC_EXC_FILE="$MAPI_ROOT/tools/.sync_excludes.txt"

if [[ ! -f "$SYNC_EXC_FILE" ]]; then
  cat >"$SYNC_EXC_FILE" <<'EOF'
.git
scratch
tools/miniconda3
tools/__HPC_NAME__/
EOF
fi

TMP_SYNC_EXC="$(mktemp)"
sed "s#__HPC_NAME__#${HPC_NAME}#g" "$SYNC_EXC_FILE" > "$TMP_SYNC_EXC"
trap 'rm -f "$TMP_SYNC_EXC"' EXIT

###############################################################################
# 6) MIRROR LOCAL → REMOTE HOME
###############################################################################

echo "[local] syncing MalariAPI → remote..."

tar -czf - \
  --exclude scratch \
  --exclude 'scratch/*' \
  --exclude envs \
  --exclude 'envs/*' \
  --exclude tools/miniconda3 \
  --exclude 'tools/miniconda3/*' \
  --exclude-from "$TMP_SYNC_EXC" \
  -C "$MAPI_ROOT" . \
  | ssh -o BatchMode=yes -o ConnectTimeout=8 "$HOST" \
      "env -u BASH_ENV bash --noprofile --norc -c 'set -euo pipefail; mkdir -p \"\$HOME/MalariAPI\"; tar -xzf - -C \"\$HOME/MalariAPI\"'"

###############################################################################
# 7) MIRROR LOCAL SCRATCH → REMOTE SCRATCH
###############################################################################

echo "[scratch] staging scratch → remote"

mkdir -p "$MAPI_ROOT/scratch"
tar -C "$MAPI_ROOT/scratch" --exclude tmp -czf - . \
  | ssh "${SSH_OPTS[@]}" "$HOST" \
      "bash --noprofile --norc -c 'set -euo pipefail; mkdir -p \"$REMOTE_SCRATCH/scratch\"; tar -xzf - -C \"$REMOTE_SCRATCH/scratch\"'"

###############################################################################
# 8) REMOTE SUBMIT WRAPPER
###############################################################################
# deprecated / removed

###############################################################################
# 9) INSTALL LOCAL HPC TOOL WRAPPERS
###############################################################################
install_local_hpc_tools() {
  local name="$1"
  local tmpl_root="$MAPI_ROOT/tools/hpc_templates"
  local dest_dir="$MAPI_ROOT/tools/$name"

  if [[ ! -d "$tmpl_root" ]]; then
    echo "[install] WARNING: no templates dir at $tmpl_root; skipping local HPC wrappers for '$name'"
    return 0
  fi

  echo "[install] Installing local HPC wrappers for '$name':"
  echo "  templates: $tmpl_root/*"
  echo "  dest:      $dest_dir/"

  mkdir -p "$dest_dir"

  rsync -a "$tmpl_root"/ "$dest_dir"/

  if ls "$dest_dir"/* >/dev/null 2>&1; then
    sed -i "s/__HPC_NAME__/$name/g" "$dest_dir"/* || true
    chmod +x "$dest_dir"/* || true
  fi
}

echo "[install] Setting up local HPC tool wrappers for $HPC_NAME..."
install_local_hpc_tools "$HPC_NAME"

###############################################################################
# 10) REMOTE MINICONDA INSTALL
###############################################################################

ssh "${SSH_OPTS[@]}" "$HOST" 'bash --noprofile --norc -s' <<'RMT'
set -euo pipefail

MAPI_HOME="$HOME/MalariAPI"
REMOTE_CONDA_DIR="$MAPI_HOME/tools/miniconda3"

# 0) Detect obviously broken Miniconda (e.g. conda shebang pointing elsewhere)
if [[ -x "$REMOTE_CONDA_DIR/bin/conda" ]]; then
  first_line="$(head -n 1 "$REMOTE_CONDA_DIR/bin/conda" 2>/dev/null || true)"
  if [[ "$first_line" != *"$REMOTE_CONDA_DIR/bin/python"* ]]; then
    echo "[remote] Existing Miniconda at $REMOTE_CONDA_DIR looks corrupted (shebang: $first_line)"
    echo "[remote] Removing and reinstalling..."
    rm -rf "$REMOTE_CONDA_DIR"
  fi
fi

# 1) install Miniconda if missing
if [[ ! -x "$REMOTE_CONDA_DIR/bin/conda" ]]; then
  echo "[remote] Installing Miniconda..."
  mkdir -p "$MAPI_HOME/tools"
  cd "$MAPI_HOME/tools"
  INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
  URL="https://repo.anaconda.com/miniconda/$INSTALLER"

  if command -v curl >/dev/null 2>&1; then
    curl -fsSL -o "$INSTALLER" "$URL"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$INSTALLER" "$URL"
  else
    echo "[remote] ERROR: neither curl nor wget is available to download Miniconda." >&2
    exit 3
  fi

  bash "$INSTALLER" -b -p "$REMOTE_CONDA_DIR"
  rm -f "$INSTALLER"
fi

echo "[remote] ensuring base env update (tools/yaml/base.yml if present)"
BASE=""
for p in \
  "$MAPI_HOME/tools/yaml/base.yml" \
  "$MAPI_HOME/tools/yaml/base.yaml"
do
  [[ -f "$p" ]] && BASE="$p" && break
done

if [[ -n "$BASE" ]]; then
  "$REMOTE_CONDA_DIR/bin/conda" env update -n base -f "$BASE" || true
fi
RMT

source "$MINICONDA_HOME/etc/profile.d/conda.sh"
conda activate base
###############################################################################
# 11) REMOTE env file only (NO ~/.bashrc modifications)
###############################################################################

echo "[remote] installing tools/mapi_remote_env.sh (no bashrc edits)"

ssh -o BatchMode=yes -o ConnectTimeout=8 "$HOST" 'bash --noprofile --norc -s' <<'RMT'
set -euo pipefail
MAPI_ROOT="$HOME/MalariAPI"
TOOLS_DIR="$MAPI_ROOT/tools"
mkdir -p "$TOOLS_DIR"

cat >"$TOOLS_DIR/mapi_remote_env.sh" <<'EOF_ENV'
# Auto-generated by MAPI installer (remote)
export MAPI_ROOT="$HOME/MalariAPI"
export REMOTE_MAPI_HOME="$HOME/MalariAPI"
export REMOTE_SCRATCH="${REMOTE_SCRATCH:-/scratch/$USER/MalariAPI}"
export MAPI_SCRATCH="$MAPI_ROOT/scratch"
EOF_ENV

chmod +x "$TOOLS_DIR/mapi_remote_env.sh"
RMT

###############################################################################
# DONE 🎉
###############################################################################

echo
echo "HPC init complete ✅"
echo "Local scratch:  $MAPI_ROOT/scratch/"
echo "Remote scratch: $REMOTE_SCRATCH/scratch/"
echo "Use:  mapi $HPC_NAME submit|status|cancel|look|peak|pull|push|sync"
