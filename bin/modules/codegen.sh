#!/usr/bin/env bash
# =========================================================
# MAPI MODULE: codegen
# =========================================================
# Description:
#   Generate MAPI modules/ from scripts or code text (static analysis + optional LLM).
#
# Usage:
#   mapi codegen --in path/to/script.py --name my_module --out-root ~/MalariAPI
#
# Notes:
#   - This is the user-facing wrapper. The implementation lives in tools/ai/.
# =========================================================
set -euo pipefail

log(){ echo "[mapi codegen] $*"; }
die(){ echo "[mapi codegen] ERROR: $*" >&2; exit 2; }

# Resolve MAPI_ROOT from: <root>/bin/modules/codegen.sh
SELF="$(readlink -f "$0" 2>/dev/null || python3 -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' "$0")"
MAPI_ROOT="$(cd "$(dirname "$SELF")/../.." && pwd)"
export MAPI_ROOT

log "MAPI_ROOT=$MAPI_ROOT"

TOOL="$MAPI_ROOT/tools/ai/mapi_codegen.py"
[[ -f "$TOOL" ]] || die "Missing tool: $TOOL"

log "Dispatching to tools/ai/mapi_codegen.py"
exec python3 "$TOOL" "$@"
