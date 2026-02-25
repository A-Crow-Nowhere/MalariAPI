#!/usr/bin/env bash
set -euo pipefail

dir="${1:-}"
list="${2:-}"

if [[ -z "$dir" ]]; then
  echo "Usage: count_files <dir> [-l]" >&2
  exit 2
fi
if [[ ! -d "$dir" ]]; then
  echo "ERROR: not a directory: $dir" >&2
  exit 2
fi

n="$(find "$dir" -maxdepth 1 -type f | wc -l | tr -d ' ')"
echo "N_FILES=$n"

if [[ "$list" == "-l" ]]; then
  # Prefer portable du usage; no GNU-only flags unless you know you have them
  find "$dir" -maxdepth 1 -type f -exec du -h {} +
fi

exit 0
