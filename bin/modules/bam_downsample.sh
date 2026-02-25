#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# bam_downsample (MAPI module)
# Downsample BAM(s) to a target number of alignment records.
# Samples ALL records present (primary/secondary/supp/split/disco/etc).
# ============================================================

VERSION="0.2.0"

# ----------------------------
# Defaults
# ----------------------------
IN_BAM=""
IN_DIR=""
PATTERN="*.bam"

OUT_BAM=""
OUT_DIR="."
TARGET_READS=""
SEED=42
THREADS=4
TOLERANCE=0.01
MAX_ITERS=4
INDEX=1
ALLOW_OVERWRITE=0
DRYRUN=0
KEEP_TMP=0
VERBOSE=0

# ----------------------------
# Helpers
# ----------------------------
die(){ echo "[bam_downsample] ERROR: $*" >&2; exit 1; }
log(){ echo "[bam_downsample] $*" >&2; }
dbg(){ [[ "$VERBOSE" -eq 1 ]] && echo "[bam_downsample] DEBUG: $*" >&2 || true; }

usage() {
  cat <<EOF
bam_downsample v$VERSION

Downsample BAM(s) to a target number of alignment records (samtools view -c).
Includes ALL record types present in the BAM.

Usage (single):
  mapi modules bam_downsample --bam in.bam --target 2000000 --out-dir out/

Usage (directory):
  mapi modules bam_downsample --dir bam_dir --target 2000000 --out-dir out/ [--pattern "*.bam"]

Required:
  --target INT          Target number of alignment records

Input (choose one):
  --bam PATH            Input BAM
  --dir DIR             Process all BAMs in DIR (non-recursive)
  --pattern GLOB        File glob within --dir (default: *.bam)

Output control:
  --out PATH            Output BAM path (single-mode only; overrides folder scheme)
  --out-dir DIR         Output base directory (default: .)
  --no-index            Do not write .bai index

Sampling control:
  --seed INT            Seed used by samtools (-s SEED.FRAC) (default: 42)
  --threads INT         Threads for samtools (default: 4)
  --tolerance FLOAT     Relative tolerance (default: 0.01 = 1%)
  --max-iters INT       Max auto-tune iterations (default: 4)

Behavior:
  --allow-overwrite     Overwrite existing outputs
  --keep-tmp            Keep temporary BAMs
  --dry-run             Print what would happen
  --verbose             Extra logging
  -h, --help            Show this help

Notes:
  - Counts/targets are alignment RECORDS, not read pairs.
  - Directory mode writes: OUT_DIR/<root>/<root>.downsampled.n<TARGET>.bam
EOF
}

bam_complete() {
  local bam="$1"
  [[ -s "$bam" ]] || return 1
  [[ -s "${bam}.bai" || -s "${bam%.bam}.bai" ]] || return 1
  return 0
}

ds_one_bam() {
  local in_bam="$1"
  local root="$2"
  local out_bam="$3"
  local out_dir_for_tmp="$4"

  [[ -s "$in_bam" ]] || { log "Skip (missing/empty): $in_bam"; return 0; }

  mkdir -p "$(dirname "$out_bam")"

  if [[ -s "$out_bam" && "$ALLOW_OVERWRITE" -eq 0 ]]; then
    if [[ "$INDEX" -eq 0 ]]; then
      log "Exists (skip): $out_bam"
      return 0
    fi
    if bam_complete "$out_bam"; then
      log "Exists (skip): $out_bam (and index)"
      return 0
    fi
  fi

  log "==> [$root] Counting records: $in_bam"
  local total
  total="$(samtools view -@ "$THREADS" -c "$in_bam")" || die "samtools count failed: $in_bam"
  [[ "$total" =~ ^[0-9]+$ ]] || die "Bad count for $in_bam"
  log "    input records: $total"

  if [[ "$total" -eq 0 ]]; then
    log "    skip: 0 records"
    return 0
  fi

  if [[ "$TARGET_READS" -ge "$total" ]]; then
    log "    target ($TARGET_READS) >= total ($total): copying"
    if [[ "$DRYRUN" -eq 1 ]]; then
      log "    [dry-run] cp $in_bam $out_bam ; index=$INDEX"
      return 0
    fi
    cp -f "$in_bam" "$out_bam"
    if [[ "$INDEX" -eq 1 ]]; then
      # Older samtools may not support -f. Remove existing index then re-index.
      if [[ "$ALLOW_OVERWRITE" -eq 1 ]]; then
        rm -f "${out_bam}.bai" "${out_bam%.bam}.bai" "${out_bam}.csi" "${out_bam%.bam}.csi" || true
      fi
      samtools index -@ "$THREADS" "$out_bam"
    fi
    return 0
  fi

  local frac
  frac="$(awk -v t="$TARGET_READS" -v n="$total" 'BEGIN{f=t/n; if(f>1)f=1; if(f<0)f=0; printf "%.6f", f}')"
  dbg "    initial fraction: $frac"

  local tmpdir
  tmpdir="$(mktemp -d "${out_dir_for_tmp%/}/.bam_downsample.${root}.XXXXXX")"
  local best_bam="" best_diff_abs="" best_count=""

  local iter=1
  while [[ "$iter" -le "$MAX_ITERS" ]]; do
    local tmpbam="${tmpdir}/${root}.iter${iter}.bam"


    if [[ "$DRYRUN" -eq 1 ]]; then
      log "    [dry-run] samtools view -@ $THREADS -b -s ${SEED}.${frac} '$in_bam' -o '$tmpbam'"
      log "    [dry-run] samtools view -@ $THREADS -c '$tmpbam'"
      break
    fi

	# samtools -s expects SEED.FRACDIGITS (e.g., 42.763197), not 42.0.763197
	frac_digits="$frac"
	frac_digits="${frac_digits#0.}"   # if "0.763197" -> "763197"
	frac_digits="${frac_digits#.}"    # if ".763197"  -> "763197" (extra safety)
	sarg="${SEED}.${frac_digits}"
	
	samtools view -@ "$THREADS" -b -s "$sarg" "$in_bam" -o "$tmpbam"

    local outn
    outn="$(samtools view -@ "$THREADS" -c "$tmpbam")"
    log "      got: $outn (target: $TARGET_READS)"

    local diff_abs
    diff_abs="$(awk -v a="$outn" -v b="$TARGET_READS" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.0f", d}')"

    if [[ -z "$best_diff_abs" || "$diff_abs" -lt "$best_diff_abs" ]]; then
      best_diff_abs="$diff_abs"
      best_bam="$tmpbam"
      best_count="$outn"
    fi

    local rel_ok
    rel_ok="$(awk -v a="$outn" -v b="$TARGET_READS" -v tol="$TOLERANCE" 'BEGIN{
      d=a-b; if(d<0)d=-d;
      r=d/b;
      if(r<=tol) print 1; else print 0
    }')"

    if [[ "$rel_ok" -eq 1 ]]; then
      log "      within tolerance (�$TOLERANCE): accept"
      break
    fi

    frac="$(awk -v f="$frac" -v t="$TARGET_READS" -v o="$outn" 'BEGIN{
      nf=f*(t/o);
      if(nf>1) nf=1;
      if(nf<0) nf=0;
      printf "%.6f", nf
    }')"
    dbg "      updated fraction -> $frac"

    iter=$((iter+1))
  done

  if [[ "$DRYRUN" -eq 1 ]]; then
    [[ "$KEEP_TMP" -eq 0 ]] && rm -rf "$tmpdir" || true
    return 0
  fi

  [[ -n "$best_bam" && -s "$best_bam" ]] || die "No sampled BAM produced for $in_bam"

  log "    best: $best_count (diff=$best_diff_abs). Writing: $out_bam"
  cp -f "$best_bam" "$out_bam"

  if [[ "$INDEX" -eq 1 ]]; then
    # Older samtools may not support -f. Remove existing index then re-index.
    if [[ "$ALLOW_OVERWRITE" -eq 1 ]]; then
      rm -f "${out_bam}.bai" "${out_bam%.bam}.bai" "${out_bam}.csi" "${out_bam%.bam}.csi" || true
    fi
    samtools index -@ "$THREADS" "$out_bam"
  fi

  if [[ "$KEEP_TMP" -eq 0 ]]; then
    rm -rf "$tmpdir" || true
  else
    log "    kept tmpdir: $tmpdir"
  fi

  log "    done."
}

# ----------------------------
# Arg parsing
# ----------------------------
if [[ $# -eq 0 ]]; then usage; exit 0; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) IN_BAM="${2:-}"; shift 2;;
    --dir) IN_DIR="${2:-}"; shift 2;;
    --pattern) PATTERN="${2:-}"; shift 2;;
    --target) TARGET_READS="${2:-}"; shift 2;;
    --out) OUT_BAM="${2:-}"; shift 2;;
    --out-dir) OUT_DIR="${2:-}"; shift 2;;
    --seed) SEED="${2:-}"; shift 2;;
    --threads) THREADS="${2:-}"; shift 2;;
    --tolerance) TOLERANCE="${2:-}"; shift 2;;
    --max-iters) MAX_ITERS="${2:-}"; shift 2;;
    --no-index) INDEX=0; shift 1;;
    --allow-overwrite) ALLOW_OVERWRITE=1; shift 1;;
    --keep-tmp) KEEP_TMP=1; shift 1;;
    --dry-run) DRYRUN=1; shift 1;;
    --verbose) VERBOSE=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# ----------------------------
# Validate
# ----------------------------
[[ -n "$TARGET_READS" ]] || die "--target is required"
[[ "$TARGET_READS" =~ ^[0-9]+$ ]] || die "--target must be an integer"
mkdir -p "$OUT_DIR"

# input mode checks
if [[ -n "$IN_BAM" && -n "$IN_DIR" ]]; then
  die "Use only one of --bam or --dir"
fi
if [[ -z "$IN_BAM" && -z "$IN_DIR" ]]; then
  die "Provide --bam or --dir"
fi

# ----------------------------
# Run single mode
# ----------------------------
if [[ -n "$IN_BAM" ]]; then
  [[ -s "$IN_BAM" ]] || die "Input BAM not found or empty: $IN_BAM"

  if [[ -z "$OUT_BAM" ]]; then
    root="$(basename "$IN_BAM")"
    root="${root%.bam}"
    out_sub="${OUT_DIR%/}/${root}"
    OUT_BAM="${out_sub}/${root}.downsampled.n${TARGET_READS}.bam"
  fi

  root="$(basename "$IN_BAM")"
  root="${root%.bam}"
  ds_one_bam "$IN_BAM" "$root" "$OUT_BAM" "$OUT_DIR"
  exit 0
fi

# ----------------------------
# Run directory mode
# ----------------------------
[[ -d "$IN_DIR" ]] || die "Directory not found: $IN_DIR"

shopt -s nullglob
bams=( "$IN_DIR"/$PATTERN )
shopt -u nullglob

if [[ "${#bams[@]}" -eq 0 ]]; then
  die "No files matched: $IN_DIR/$PATTERN"
fi

log "Directory mode: ${#bams[@]} file(s) matched '$PATTERN' in $IN_DIR"

for bam in "${bams[@]}"; do
  bn="$(basename "$bam")"
  root="${bn%.bam}"
  out_sub="${OUT_DIR%/}/${root}"
  out_bam="${out_sub}/${root}.downsampled.n${TARGET_READS}.bam"
  ds_one_bam "$bam" "$root" "$out_bam" "$OUT_DIR"
done

log "All done."
