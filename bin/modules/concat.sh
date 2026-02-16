#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# concat_keep_header
# Concatenate TSV files in a directory, keeping one header
# ============================================================

say() { echo "[concat_keep_header] $*"; }
die() { say "ERROR: $*"; exit 1; }

usage() {
  cat >&2 <<'EOF'
Usage:
  concat_keep_header <DIR> [--out FILE]

Description:
  Concatenates all TSV/TSV.GZ files in DIR into a single file,
  keeping exactly one header (from the first file alphabetically).

Arguments:
  DIR                 Directory containing input files

Options:
  --out FILE          Output file (default: <DIR>/merged.tsv)
  -h, --help          Show this help

Notes:
  - Skips empty files
  - Supports .tsv and .tsv.gz
  - Header assumed to be the first line
EOF
}

# -----------------------------
# Parse args
# -----------------------------
[[ $# -eq 0 ]] && { usage; exit 1; }

DIR=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0;;
    --out) OUT="${2:-}"; shift 2;;
    *)
      [[ -z "$DIR" ]] && DIR="$1" || die "Unexpected argument: $1"
      shift
      ;;
  esac
done

[[ -d "$DIR" ]] || die "Not a directory: $DIR"

if [[ -z "$OUT" ]]; then
  OUT="$DIR/merged.tsv"
fi

say "Input dir: $DIR"
say "Output:    $OUT"

rm -f "$OUT"

# -----------------------------
# Find files
# -----------------------------
shopt -s nullglob
FILES=( "$DIR"/*.tsv "$DIR"/*.tsv.gz )

[[ ${#FILES[@]} -gt 0 ]] || die "No .tsv or .tsv.gz files found in $DIR"

# -----------------------------
# Header
# -----------------------------
FIRST="${FILES[0]}"
say "Using header from: $(basename "$FIRST")"

if [[ "$FIRST" == *.gz ]]; then
  zcat "$FIRST" | head -n 1 > "$OUT"
else
  head -n 1 "$FIRST" > "$OUT"
fi

# -----------------------------
# Body
# -----------------------------
for f in "${FILES[@]}"; do
  [[ -s "$f" ]] || continue

  if [[ "$f" == *.gz ]]; then
    zcat "$f" | tail -n +2 >> "$OUT"
  else
    tail -n +2 "$f" >> "$OUT"
  fi
done

say "Done. Wrote $(wc -l < "$OUT") lines."
