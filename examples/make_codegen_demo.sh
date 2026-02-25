#!/usr/bin/env bash
set -euo pipefail

# Demo: generate a module called "filter_tsv" from examples/filter_tsv.py
# Run this from your MalariAPI root (it should contain modules/, templates/, tools/).

ROOT="$(pwd)"
if [[ ! -d "$ROOT/bin/modules" || ! -d "$ROOT/tools/templates" || ! -d "$ROOT/tools" ]]; then
  echo "[demo] Run this from your MalariAPI root (needs modules/, templates/, tools/)." >&2
  exit 2
fi

echo "[demo] Generating module 'filter_tsv' (offline backend)"
python3 tools/ai/mapi_codegen.py \
  --in examples/filter_tsv.py \
  --name filter_tsv \
  --description "Filter TSV rows by numeric minimum across all numeric columns" \
  --out-root "$ROOT" \
  --backend offline \
  --print-spec

echo
echo "[demo] Generated:"
echo "  modules/filter_tsv/run"
echo "  modules/filter_tsv/module.yml"
echo "  tools/yaml/filter_tsv.yml"
echo
echo "[demo] Next steps (typical):"
echo "  1) Build env: conda env create -p envs/filter_tsv -f tools/yaml/filter_tsv.yml"
echo "  2) Edit modules/filter_tsv/run entrypoint if needed (look for {SCRIPT} placeholder)"
echo "  3) Try: mapi filter_tsv --help"
