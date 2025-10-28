#!/usr/bin/env bash
set -euo pipefail
: "${MAPI_CONDA_ROOT:="$HOME/tools/miniconda3"}"
export PATH="$MAPI_CONDA_ROOT/bin:$PATH"
export CONDA_ENVS_DIRS="$MAPI_CONDA_ROOT/envs"
if [[ -f "$MAPI_CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1090
  source "$MAPI_CONDA_ROOT/etc/profile.d/conda.sh"
fi
