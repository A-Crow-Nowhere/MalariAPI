cat > ~/tools/mapi/lib.sh <<'BASH'
#!/usr/bin/env bash
set -euo pipefail

: "${ENV_YAML_DIR:=$HOME/envs/yaml}"

_mamba_or_conda() { command -v mamba >/dev/null 2>&1 && echo mamba || echo conda; }

ensure_env() {
  local env_name="$1" yaml_base="$2" tool="$(_mamba_or_conda)"
  if ! "$tool" env list | awk '{print $1}' | grep -qx "$env_name"; then
    echo "[mapi] Creating env '$env_name' from $ENV_YAML_DIR/${yaml_base}.yaml"
    "$tool" env create -n "$env_name" -f "$ENV_YAML_DIR/${yaml_base}.yaml"
  fi
}

run_in_env() { local env_name="$1"; shift; conda run -n "$env_name" "$@"; }
BASH
chmod +x ~/tools/mapi/lib.sh
