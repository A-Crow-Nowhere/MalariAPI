#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
bw_dir_to_bedgraph.sh - convert all BigWigs in a directory to bedGraph (BED4: chrom start end value)

Usage:
  bw_dir_to_bedgraph.sh --dir <BW_DIR> --out <OUT_DIR> [--gz] [--threads N]

Options:
  --dir        Input directory containing .bw/.bigWig files
  --out        Output directory for .bedGraph (or .bedGraph.gz if --gz)
  --gz         Gzip outputs
  --threads    Parallel jobs (default: 4)

Notes:
  - Requires: bigWigToBedGraph (UCSC tools), gzip (optional), xargs
  - Output format: chrom  start  end  value   (raw signal from BigWig)
EOF
}

DIR=""
OUT=""
GZ=0
THREADS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) DIR="${2:-}"; shift 2;;
    --out) OUT="${2:-}"; shift 2;;
    --gz)  GZ=1; shift;;
    --threads) THREADS="${2:-}"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

[[ -n "$DIR" ]] || { echo "ERROR: --dir is required" >&2; usage; exit 1; }
[[ -n "$OUT" ]] || { echo "ERROR: --out is required" >&2; usage; exit 1; }
[[ -d "$DIR" ]] || { echo "ERROR: not a directory: $DIR" >&2; exit 1; }

command -v bigWigToBedGraph >/dev/null 2>&1 || {
  echo "ERROR: bigWigToBedGraph not found in PATH" >&2
  echo "Install UCSC tools (bigWigToBedGraph) and try again." >&2
  exit 1
}

mkdir -p "$OUT"

# Find bigwigs
mapfile -t BWS < <(find "$DIR" -maxdepth 1 -type f \( -iname "*.bw" -o -iname "*.bigwig" -o -iname "*.bigWig" \) | sort)
(( ${#BWS[@]} > 0 )) || { echo "No .bw/.bigWig files found in: $DIR" >&2; exit 0; }

export OUT GZ
convert_one() {
  local bw="$1"
  local base
  base="$(basename "$bw")"
  base="${base%.*}"

  local out="$OUT/${base}.bedGraph"
  if [[ "$GZ" -eq 1 ]]; then
    out="${out}.gz"
    echo "[bw2bed] $bw -> $out"
    bigWigToBedGraph "$bw" /dev/stdout | gzip -c > "$out"
  else
    echo "[bw2bed] $bw -> $out"
    bigWigToBedGraph "$bw" "$out"
  fi
}
export -f convert_one

# Parallelize
printf "%s\0" "${BWS[@]}" | xargs -0 -n 1 -P "$THREADS" bash -lc 'convert_one "$@"' _

echo "[bw2bed] Done. Outputs in: $OUT"
