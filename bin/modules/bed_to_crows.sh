#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# bed_to_svcrows.sh
# Convert BED/TSV-like SV calls into SVCROWS TSV schema
# ============================================================

say() { echo "[bed_to_svcrows] $*"; }
die() { say "ERROR: $*"; exit 1; }

usage() {
  cat >&2 <<'EOF'
Usage:
  bed_to_svcrows.sh --in <calls.tsv[.gz]> --out <out.tsv[.gz]> [options]

Required:
  --in FILE                Input TSV/BED-like file (optionally gz)
  --out FILE               Output SVCROWS TSV (optionally gz)

Options:
  --sample STR             Value to write into Var3 (default: Sample)
  --source STR             Source/caller tag (added into Var1) (default: bed)
  --type-col INT           1-based column holding SV type (default: 4)
  --id-col INT             1-based column holding a pre-existing ID/name (default: 0 = none)
  --su-col INT             1-based column holding support/reads/SU (default: 0 = none -> 0)
  --require-bounds         Require Start/End numeric and Start < End (default: on)
  --keep-types STR         Comma list of types to keep after normalization (default: DEL,DUP,INV)
  --header auto|yes|no     Input header handling (default: auto)
  --out-id-prefix STR      Prefix for numeric ID field (default: empty; ID is 1..N)
  --drop-nonstandard       Drop rows with unknown/unhandled types (default: on)

Notes:
- Normalizes type synonyms:
    duplication/dup -> DUP
    deletion/del   -> DEL
    inversion/inv  -> INV
- Length is computed as (End-Start). For DEL, output Length is negative.
- Var1 includes a VCF-like INFO string plus the original row fields where possible.
- Var2 is a lightweight genotype/support field: GT=./.;SU=<NumReads>
EOF
}

IN=""
OUT=""
SAMPLE="Sample"
SOURCE="bed"
TYPE_COL=4
ID_COL=0
SU_COL=0
REQUIRE_BOUNDS=1
KEEP_TYPES="DEL,DUP,INV"
HEADER_MODE="auto"
OUT_ID_PREFIX=""
DROP_NONSTANDARD=1

# -----------------------------
# Arg parse
# -----------------------------
[[ $# -eq 0 ]] && { usage; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in) IN="${2:-}"; shift 2;;
    --out) OUT="${2:-}"; shift 2;;
    --sample) SAMPLE="${2:-}"; shift 2;;
    --source) SOURCE="${2:-}"; shift 2;;
    --type-col) TYPE_COL="${2:-}"; shift 2;;
    --id-col) ID_COL="${2:-}"; shift 2;;
    --su-col) SU_COL="${2:-}"; shift 2;;
    --require-bounds) REQUIRE_BOUNDS=1; shift 1;;
    --keep-types) KEEP_TYPES="${2:-}"; shift 2;;
    --header) HEADER_MODE="${2:-}"; shift 2;;
    --out-id-prefix) OUT_ID_PREFIX="${2:-}"; shift 2;;
    --drop-nonstandard) DROP_NONSTANDARD=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

[[ -n "$IN" ]] || die "--in is required"
[[ -n "$OUT" ]] || die "--out is required"
[[ -f "$IN" ]] || die "Input not found: $IN"

# -----------------------------
# Resolve helper
# -----------------------------
THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER_DIR="$THIS_DIR/.svcrows"
PY="$HELPER_DIR/bed_to_svcrows.py"
[[ -f "$PY" ]] || die "Missing helper script: $PY"

# -----------------------------
# Run
# -----------------------------
say "IN: $IN"
say "OUT: $OUT"
say "sample: $SAMPLE | source: $SOURCE | type-col: $TYPE_COL | id-col: $ID_COL | su-col: $SU_COL | header: $HEADER_MODE"

python "$PY" \
  --in "$IN" \
  --out "$OUT" \
  --sample "$SAMPLE" \
  --source "$SOURCE" \
  --type-col "$TYPE_COL" \
  --id-col "$ID_COL" \
  --su-col "$SU_COL" \
  --keep-types "$KEEP_TYPES" \
  --header "$HEADER_MODE" \
  --out-id-prefix "$OUT_ID_PREFIX" \
  $( [[ "$REQUIRE_BOUNDS" -eq 1 ]] && echo --require-bounds ) \
  $( [[ "$DROP_NONSTANDARD" -eq 1 ]] && echo --drop-nonstandard )

say "Done."
