#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# bwa_align - MAPI module (clean defaults + audit + guardrails)
#
# Design goals (dummy-proof):
#   1) The MAIN BAM is a true superset:
#        <sample>.<outkey>.aligned.sorted.bam includes mapped + unmapped records.
#      - This module DOES NOT offer any "samtools view -q/-f/-F" filtering upstream.
#      - Any quality/MAPQ filtering is only allowed on *derived* BAMs.
#
#   2) Derived BAMs are explicit subsets of ALIGNED_SORTED (unmapped/supp/primary/secondary),
#      and optional samblaster side-files (splitters/discordant) are produced ONLY if requested.
#
#   3) Defaults are "clean":
#      - bwa mem: only -t and auto RG (unless --no-rg). No hidden -a/-T/-k/-Y/-M.
#      - samblaster: OFF by default; when ON, defaults are conservative and do NOT drop reads.
#
#   4) --audit-counts writes a per-sample TSV that makes read accounting obvious.
#
# ============================================================

MAPI_ROOT="${MAPI_ROOT:-$HOME/MalariAPI}"
GENOMES_DIR="${GENOMES_DIR:-$MAPI_ROOT/genomes}"
THREADS="${THREADS:-4}"

# Clean defaults: no implicit bwa tuning.
BWA_ARGS_DEFAULT=()

# Samblaster defaults (ONLY used if samblaster is engaged).
# IMPORTANT: no --ignoreUnmated by default (that can drop records in-stream).
SAMBLASTER_ARGS_DEFAULT=( -M --maxSplitCount 10 --maxUnmappedBases 100 --minIndelSize 1 --minNonOverlap 1 --minClipSize 5 )

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
have() { command -v "$1" >/dev/null 2>&1; }
ensure_dir() { mkdir -p "$1"; }

print_help() {
  cat <<'EOF'
Usage:
  mapi modules bwa_align --r1 R1.fq.gz --r2 R2.fq.gz (--ref REF.fa | --genome TOKEN) [output opts] [flags]

Required:
  --r1 PATH
  --r2 PATH
  One of:
    --ref PATH
    --genome TOKEN

Output options:
  --out-dir DIR            Directory to write outputs (NO filename). If omitted:
                           defaults to: <dirname(R1)>/bwa_align_out
  --sample NAME            Override auto sample name inferred from R1
  --out-key KEY            Label added into filenames (default: bwa)

Back-compat (deprecated):
  --out-prefix PREFIX      Old behavior: use exact prefix (path+name). Still supported.

Core options:
  --threads N
  --bwa-extra "STRING"            Advanced passthrough to bwa mem (optional).
  --samblaster-extra "STRING"     Advanced passthrough to samblaster (optional).
  --samblaster-ignore-unmated     Add samblaster --ignoreUnmated (NOT default).

Read group (auto by default; safe on Slurm):
  --no-rg
  --rg-force
  --rg-id STR              (default: --sample)
  --rg-sm STR              (default: --sample)
  --rg-pl STR              (default: ILLUMINA)
  --rg-lb STR              (optional)
  --rg-pu STR              (optional)

Samblaster (OFF by default; enabled if any of these are requested):
  --use-samblaster
  --emit-discordant
  --emit-splitters

Derived splits from aligned BAM (sorted+indexed):
  --emit-unmapped
  --emit-supplementary
  --emit-primary
  --emit-secondary

MAPQ filtering (DERIVED outputs only; main aligned BAM is NEVER MAPQ-filtered):
  --min-mapq-primary N
  --min-mapq-secondary N
  --min-mapq-supplementary N
  --min-mapq-discordant N          (post-hoc filter of discordant output)
  --min-mapq-splitters N           (post-hoc filter of splitters output)

Reference/index behavior:
  --genomes-dir DIR
  --no-index

Audit / accounting:
  --audit-counts            Write <sample>.<outkey>.audit.tsv in OUT_DIR with:
                            FASTQ reads (R1/R2),
                            aligned BAM totals and key flag counts,
                            expected primary-like count (2*read_pairs),
                            and basic consistency checks.

Resume / recompute:
  --resume                  (default) Skip completed outputs; reuse unsorted intermediates if present
  --no-resume               Recompute everything from scratch (overwrite)
  --keep-tmp                Do not delete tmp intermediates (default keeps them; present for clarity)

Sorting controls:
  --sort-threads N          Threads for samtools sort (default: min(--threads,4))
  --sort-mem 256M           Memory per thread for samtools sort -m (default: 256M)
  --sort-tmp DIR            Directory for sort temp (default: OUTDIR/.tmp_bwa_align/sort_tmp)
  --sort-extra "STRING"     Extra args appended to samtools sort (optional)

EOF
}

# shell-style split of an arg string into array
split_args() {
  local s="${1:-}"
  if [[ -z "$s" ]]; then
    _SPLIT_ARGS=()
    return 0
  fi
  # shellcheck disable=SC2206
  _SPLIT_ARGS=($s)
}

infer_sample_from_r1() {
  local r1="$1"
  local b
  b="$(basename "$r1")"
  b="${b%.gz}"
  b="${b%.bgz}"
  b="${b%.fastq}"
  b="${b%.fq}"
  b="${b%_R1}"
  b="${b%_r1}"
  b="${b%_1}"
  b="${b%.R1}"
  b="${b%.r1}"
  echo "$b"
}

ref_lookup_first() {
  local token="$1"
  local hit=""
  hit="$(find "$GENOMES_DIR" -maxdepth 4 -type f \
      \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" -o -iname "*.fa.gz" -o -iname "*.fasta.gz" -o -iname "*.fna.gz" \) \
      2>/dev/null | grep -i "$token" | head -n 1 || true)"
  [[ -n "$hit" ]] || return 1
  echo "$hit"
}

prepare_reference() {
  local ref_in="$1"
  local outdir="$2"
  [[ -f "$ref_in" ]] || die "Reference not found: $ref_in"
  ensure_dir "$outdir"
  local cache_dir="$outdir/.ref_cache"
  ensure_dir "$cache_dir"

  if [[ "$ref_in" == *.gz ]]; then
    local base
    base="$(basename "$ref_in")"
    base="${base%.gz}"
    local ref_out="$cache_dir/$base"
    if [[ ! -s "$ref_out" ]]; then
      echo "==> Decompressing reference to cache: $ref_out"
      gzip -cd "$ref_in" > "$ref_out"
    fi
    echo "$ref_out"
  else
    echo "$ref_in"
  fi
}

ensure_indexes() {
  local ref="$1"
  local build_ok="$2"

  local bwa_ok=1
  [[ -s "${ref}.bwt" || -s "${ref}.0123.bwt" ]] || bwa_ok=0

  local fai_ok=1
  [[ -s "${ref}.fai" ]] || fai_ok=0

  if [[ "$bwa_ok" -eq 1 && "$fai_ok" -eq 1 ]]; then
    return 0
  fi
  [[ "$build_ok" -eq 1 ]] || die "Missing indexes for reference ($ref) and --no-index was set."

  if [[ "$bwa_ok" -eq 0 ]]; then
    echo "==> Building BWA index for: $ref"
    bwa index "$ref"
  fi
  if [[ "$fai_ok" -eq 0 ]]; then
    echo "==> Building samtools faidx for: $ref"
    samtools faidx "$ref"
  fi
}

bam_complete() {
  local bam="$1"
  [[ -s "$bam" && -s "${bam}.bai" ]]
}

# globals set later:
SORT_THREADS=""
SORT_MEM="256M"
SORT_TMP=""
SORT_EXTRA=""

sort_and_index_bam() {
  local in_bam="$1"
  local out_bam="$2"

  [[ -s "$in_bam" ]] || die "Expected BAM not found or empty: $in_bam"

  if bam_complete "$out_bam"; then
    echo "==> Exists (skip): $out_bam"
    return 0
  fi

  ensure_dir "$(dirname "$out_bam")"
  ensure_dir "$SORT_TMP"

  local sort_args=(-@ "$SORT_THREADS" -m "$SORT_MEM" -T "$SORT_TMP/$(basename "$out_bam").tmp")
  split_args "$SORT_EXTRA"
  ((${#_SPLIT_ARGS[@]})) && sort_args+=("${_SPLIT_ARGS[@]}")

  echo "==> Sorting+indexing: $out_bam  (threads=$SORT_THREADS mem=$SORT_MEM tmp=$SORT_TMP)"
  samtools sort "${sort_args[@]}" -o "$out_bam" "$in_bam"
  samtools index -@ "$SORT_THREADS" "$out_bam"
}

# Helper: MAPQ args for derived views
view_mapq_args() {
  local minq="${1:-}"
  VIEW_MAPQ_ARGS=()
  if [[ -n "$minq" ]]; then
    [[ "$minq" =~ ^[0-9]+$ ]] || die "MAPQ must be integer, got: $minq"
    VIEW_MAPQ_ARGS=( -q "$minq" )
  fi
}

# FASTQ read counting (gz or plain). Returns reads = lines/4.
count_fastq_reads() {
  local fq="$1"
  [[ -s "$fq" ]] || { echo "NA"; return 0; }

  if [[ "$fq" == *.gz || "$fq" == *.bgz ]]; then
    if have zcat; then
      zcat -f "$fq" | awk 'END{printf "%.0f\n", NR/4}'
    else
      gzip -cd "$fq" | awk 'END{printf "%.0f\n", NR/4}'
    fi
  else
    awk 'END{printf "%.0f\n", NR/4}' "$fq"
  fi
}

# Audit TSV writer (single-row, per-sample)
write_audit() {
  local audit_tsv="$1"
  local r1_count="$2"
  local r2_count="$3"
  local bam="$4"

  local total unmapped supp secondary primary_like
  total="$(samtools view -c "$bam")"
  unmapped="$(samtools view -c -f 4 "$bam")"
  supp="$(samtools view -c -f 2048 "$bam")"
  secondary="$(samtools view -c -f 256 "$bam")"
  primary_like="$(samtools view -c -F 2304 "$bam")"   # excludes 0x100 + 0x800

  local expected_primary_like="NA"
  local audit_note="ok"

  if [[ "$r1_count" != "NA" && "$r2_count" != "NA" ]]; then
    if [[ "$r1_count" != "$r2_count" ]]; then
      audit_note="WARN: R1!=R2 (paired input mismatch?)"
    else
      # paired-end: expected primary-like records = 2 * read_pairs = 2 * R1
      expected_primary_like=$(( 2 * r1_count ))
      if [[ "$primary_like" -ne "$expected_primary_like" ]]; then
        audit_note="WARN: primary_like != 2*read_pairs (possible stream loss or unusual input)"
      fi
    fi
  fi

  {
    echo -e "sample\tout_key\tr1_reads\tr2_reads\texpected_primary_like\tbam_records_total\tbam_unmapped_f4\tbam_supp_f2048\tbam_secondary_f256\tbam_not_secondary_or_supp_F2304\tnote"
    echo -e "${SAMPLE}\t${OUT_KEY}\t${r1_count}\t${r2_count}\t${expected_primary_like}\t${total}\t${unmapped}\t${supp}\t${secondary}\t${primary_like}\t${audit_note}"
  } > "$audit_tsv"
}

# Guardrails: warn about bwa-extra that commonly confuses accounting
warn_if_suspicious_bwa_extra() {
  local s="${1:-}"
  [[ -n "$s" ]] || return 0

  if [[ "$s" =~ (^|[[:space:]])-a([[:space:]]|$) ]]; then
    warn "--bwa-extra contains -a (output all alignments). This can greatly increase BAM records and confuse accounting."
  fi
  if [[ "$s" =~ (^|[[:space:]])-T([[:space:]]|$) ]]; then
    warn "--bwa-extra contains -T. This changes alignment reporting threshold (still should not drop reads, but will change mapped/unmapped balance)."
  fi
  if [[ "$s" =~ (^|[[:space:]])-Y([[:space:]]|$) ]]; then
    warn "--bwa-extra contains -Y. This changes how split alignments are flagged (supplementary vs secondary)."
  fi
  if [[ "$s" =~ (^|[[:space:]])-M([[:space:]]|$) ]]; then
    warn "--bwa-extra contains -M. This changes marking of shorter split hits; fine, but affects downstream interpretation."
  fi
}

# ----------------------------
# Args
# ----------------------------
R1=""
R2=""
REF_PATH=""
GENOME_TOKEN=""

OUT_DIR=""
SAMPLE=""
OUT_KEY="bwa"

# Read group behavior
AUTO_RG=1
RG_FORCE=0
RG_ID=""
RG_SM=""
RG_PL="ILLUMINA"
RG_LB=""
RG_PU=""

OUT_PREFIX=""

BWA_EXTRA=""
SAMBLASTER_EXTRA=""
SAMBLASTER_IGNORE_UNMATED=0

USE_SAMBLASTER=0
EMIT_DISCORDANT=0
EMIT_SPLITTERS=0
EMIT_UNMAPPED=0
EMIT_SUPPLEMENTARY=0
EMIT_PRIMARY=0
EMIT_SECONDARY=0

# MAPQ derived filters
MIN_MAPQ_PRIMARY=""
MIN_MAPQ_SECONDARY=""
MIN_MAPQ_SUPP=""
MIN_MAPQ_DISCORDANT=""
MIN_MAPQ_SPLITTERS=""

GENOMES_DIR_OVERRIDE=""
NO_INDEX=0

AUDIT_COUNTS=0

RESUME=1
KEEP_TMP=1

SORT_THREADS_CLI=""
SORT_MEM="256M"
SORT_TMP_CLI=""
SORT_EXTRA=""

[[ $# -eq 0 ]] && { print_help; exit 0; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) print_help; exit 0 ;;

    --r1) R1="${2:-}"; shift 2 ;;
    --r2) R2="${2:-}"; shift 2 ;;
    --ref) REF_PATH="${2:-}"; shift 2 ;;
    --genome) GENOME_TOKEN="${2:-}"; shift 2 ;;

    --out-dir) OUT_DIR="${2:-}"; shift 2 ;;
    --sample) SAMPLE="${2:-}"; shift 2 ;;
    --out-key) OUT_KEY="${2:-}"; shift 2 ;;

    --out-prefix) OUT_PREFIX="${2:-}"; shift 2 ;;

    --threads) THREADS="${2:-}"; shift 2 ;;
    --bwa-extra) BWA_EXTRA="${2:-}"; shift 2 ;;
    --samblaster-extra) SAMBLASTER_EXTRA="${2:-}"; shift 2 ;;
    --samblaster-ignore-unmated) SAMBLASTER_IGNORE_UNMATED=1; shift 1 ;;

    --use-samblaster) USE_SAMBLASTER=1; shift 1 ;;
    --emit-discordant) EMIT_DISCORDANT=1; shift 1 ;;
    --emit-splitters) EMIT_SPLITTERS=1; shift 1 ;;
    --emit-unmapped) EMIT_UNMAPPED=1; shift 1 ;;
    --emit-supplementary) EMIT_SUPPLEMENTARY=1; shift 1 ;;
    --emit-primary) EMIT_PRIMARY=1; shift 1 ;;
    --emit-secondary) EMIT_SECONDARY=1; shift 1 ;;

    --min-mapq-primary) MIN_MAPQ_PRIMARY="${2:-}"; shift 2 ;;
    --min-mapq-secondary) MIN_MAPQ_SECONDARY="${2:-}"; shift 2 ;;
    --min-mapq-supplementary) MIN_MAPQ_SUPP="${2:-}"; shift 2 ;;
    --min-mapq-discordant) MIN_MAPQ_DISCORDANT="${2:-}"; shift 2 ;;
    --min-mapq-splitters) MIN_MAPQ_SPLITTERS="${2:-}"; shift 2 ;;

    --audit-counts) AUDIT_COUNTS=1; shift 1 ;;

    --genomes-dir) GENOMES_DIR_OVERRIDE="${2:-}"; shift 2 ;;
    --no-index) NO_INDEX=1; shift 1 ;;

    --no-rg) AUTO_RG=0; shift 1 ;;
    --rg-force) RG_FORCE=1; shift 1 ;;
    --rg-id) RG_ID="${2:-}"; shift 2 ;;
    --rg-sm) RG_SM="${2:-}"; shift 2 ;;
    --rg-pl) RG_PL="${2:-}"; shift 2 ;;
    --rg-lb) RG_LB="${2:-}"; shift 2 ;;
    --rg-pu) RG_PU="${2:-}"; shift 2 ;;

    --resume) RESUME=1; shift 1 ;;
    --no-resume) RESUME=0; shift 1 ;;
    --keep-tmp) KEEP_TMP=1; shift 1 ;;

    --sort-threads) SORT_THREADS_CLI="${2:-}"; shift 2 ;;
    --sort-mem) SORT_MEM="${2:-}"; shift 2 ;;
    --sort-tmp) SORT_TMP_CLI="${2:-}"; shift 2 ;;
    --sort-extra) SORT_EXTRA="${2:-}"; shift 2 ;;

    *) die "Unknown argument: $1 (try --help)" ;;
  esac
done

[[ -n "$GENOMES_DIR_OVERRIDE" ]] && GENOMES_DIR="$GENOMES_DIR_OVERRIDE"

[[ -n "$R1" && -f "$R1" ]] || die "--r1 missing or not found: $R1"
[[ -n "$R2" && -f "$R2" ]] || die "--r2 missing or not found: $R2"

if [[ -n "$REF_PATH" && -n "$GENOME_TOKEN" ]]; then
  die "Use only one of --ref or --genome"
fi
if [[ -z "$REF_PATH" && -z "$GENOME_TOKEN" ]]; then
  die "You must provide --ref PATH or --genome TOKEN"
fi

have bwa || die "bwa not found in PATH"
have samtools || die "samtools not found in PATH"

# Decide output directory + sample base
if [[ -n "$OUT_PREFIX" ]]; then
  OUT_DIR="$(dirname "$OUT_PREFIX")"
  SAMPLE="$(basename "$OUT_PREFIX")"
  [[ -n "$OUT_DIR" ]] || die "--out-prefix produced empty out dir?"
  [[ -n "$SAMPLE" ]] || die "--out-prefix produced empty sample base?"
else
  [[ -n "$OUT_DIR" ]] || OUT_DIR="$(dirname "$R1")/bwa_align_out"
  [[ -n "$SAMPLE" ]] || SAMPLE="$(infer_sample_from_r1 "$R1")"
fi

[[ -n "$OUT_KEY" ]] || die "--out-key cannot be empty"
ensure_dir "$OUT_DIR"

TMPDIR="$OUT_DIR/.tmp_bwa_align"
ensure_dir "$TMPDIR"

# Sorting defaults
if [[ -n "$SORT_THREADS_CLI" ]]; then
  SORT_THREADS="$SORT_THREADS_CLI"
else
  if [[ "$THREADS" -le 4 ]]; then SORT_THREADS="$THREADS"; else SORT_THREADS=4; fi
fi
[[ -n "$SORT_TMP_CLI" ]] && SORT_TMP="$SORT_TMP_CLI" || SORT_TMP="$TMPDIR/sort_tmp"

# Determine if we need samblaster
NEED_SAMBLASTER=0
if [[ "$USE_SAMBLASTER" -eq 1 || "$EMIT_DISCORDANT" -eq 1 || "$EMIT_SPLITTERS" -eq 1 ]]; then
  NEED_SAMBLASTER=1
  have samblaster || die "samblaster not found in PATH (needed for requested outputs)"
fi

# Resolve reference
REF_IN=""
if [[ -n "$REF_PATH" ]]; then
  REF_IN="$REF_PATH"
else
  [[ -d "$GENOMES_DIR" ]] || die "Genomes dir not found: $GENOMES_DIR"
  REF_IN="$(ref_lookup_first "$GENOME_TOKEN" || true)"
  [[ -n "$REF_IN" ]] || die "No reference match for token '$GENOME_TOKEN' under: $GENOMES_DIR"
fi

REF="$(prepare_reference "$REF_IN" "$OUT_DIR")"
BUILD_OK=1
[[ "$NO_INDEX" -eq 1 ]] && BUILD_OK=0
ensure_indexes "$REF" "$BUILD_OK"

# Output names (sorted+indexed)
base="${OUT_DIR}/${SAMPLE}.${OUT_KEY}"

ALIGNED_SORTED="${base}.aligned.sorted.bam"
DISCORDANT_SORTED="${base}.discordant.sorted.bam"
SPLITTERS_SORTED="${base}.splitters.sorted.bam"
UNMAPPED_SORTED="${base}.unmapped.sorted.bam"
SUPP_SORTED="${base}.supplementary.sorted.bam"
PRIMARY_SORTED="${base}.primary.sorted.bam"
SECONDARY_SORTED="${base}.secondary.sorted.bam"

AUDIT_TSV="${base}.audit.tsv"

# Temp/unsorted names
ALIGNED_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.aligned.unsorted.bam"
DISCORDANT_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.discordant.unsorted.bam"
SPLITTERS_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.splitters.unsorted.bam"

UNMAPPED_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.unmapped.unsorted.bam"
SUPP_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.supplementary.unsorted.bam"
PRIMARY_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.primary.unsorted.bam"
SECONDARY_UNSORTED="$TMPDIR/${SAMPLE}.${OUT_KEY}.secondary.unsorted.bam"

# Build arg arrays
BWA_ARGS=( "${BWA_ARGS_DEFAULT[@]}" )
split_args "$BWA_EXTRA"
((${#_SPLIT_ARGS[@]})) && BWA_ARGS+=( "${_SPLIT_ARGS[@]}" )

warn_if_suspicious_bwa_extra "$BWA_EXTRA"

# Auto RG injection (escaped \t style)
BWA_EXTRA_HAS_R=0
if [[ -n "$BWA_EXTRA" ]]; then
  if [[ "$BWA_EXTRA" =~ (^|[[:space:]])-R([[:space:]]|$) ]]; then
    BWA_EXTRA_HAS_R=1
  fi
fi

if [[ "$AUTO_RG" -eq 1 ]]; then
  if [[ "$BWA_EXTRA_HAS_R" -eq 1 && "$RG_FORCE" -ne 1 ]]; then
    warn "--bwa-extra contains -R; skipping auto RG (use --rg-force to override)"
  else
    [[ -n "$RG_ID" ]] || RG_ID="$SAMPLE"
    [[ -n "$RG_SM" ]] || RG_SM="$SAMPLE"
    [[ -n "$RG_PL" ]] || RG_PL="ILLUMINA"

    RG_LINE="@RG\\tID:${RG_ID}\\tSM:${RG_SM}\\tPL:${RG_PL}"
    [[ -n "$RG_LB" ]] && RG_LINE+="\\tLB:${RG_LB}"
    [[ -n "$RG_PU" ]] && RG_LINE+="\\tPU:${RG_PU}"

    BWA_ARGS+=( -R "$RG_LINE" )
  fi
fi

SAMBLASTER_ARGS=( "${SAMBLASTER_ARGS_DEFAULT[@]}" )
split_args "$SAMBLASTER_EXTRA"
((${#_SPLIT_ARGS[@]})) && SAMBLASTER_ARGS+=( "${_SPLIT_ARGS[@]}" )
if [[ "$SAMBLASTER_IGNORE_UNMATED" -eq 1 ]]; then
  SAMBLASTER_ARGS+=( --ignoreUnmated )
  warn "samblaster --ignoreUnmated enabled; this may drop records if input pairing is irregular."
fi

echo "==> Reference: $REF"
echo "==> R1: $R1"
echo "==> R2: $R2"
echo "==> Threads: $THREADS"
echo "==> Sort threads: $SORT_THREADS"
echo "==> Sort mem: $SORT_MEM"
echo "==> OUT_DIR: $OUT_DIR"
echo "==> SAMPLE: $SAMPLE"
echo "==> OUT_KEY: $OUT_KEY"
echo "==> Main output: $ALIGNED_SORTED"
echo "==> NOTE: Main aligned BAM is NEVER MAPQ-filtered."

# ----------------------------
# Step 1: Produce aligned.unsorted if needed
# ----------------------------
need_aln=1
if [[ "$RESUME" -eq 1 ]]; then
  if bam_complete "$ALIGNED_SORTED"; then
    echo "==> Aligned sorted exists; skipping alignment+sort for main."
    need_aln=0
  elif [[ -s "$ALIGNED_UNSORTED" ]]; then
    echo "==> Found existing aligned unsorted; will only sort/index."
    need_aln=0
  fi
fi

# Audit: count FASTQs early (can be slow; only when requested)
R1_COUNT="NA"
R2_COUNT="NA"
if [[ "$AUDIT_COUNTS" -eq 1 ]]; then
  echo "==> Audit: counting FASTQ reads (R1/R2) ..."
  R1_COUNT="$(count_fastq_reads "$R1" || echo "NA")"
  R2_COUNT="$(count_fastq_reads "$R2" || echo "NA")"
fi

if [[ "$need_aln" -eq 1 ]]; then
  if [[ "$NEED_SAMBLASTER" -eq 1 ]]; then
    SB_OUT_FLAGS=()
    [[ "$EMIT_SPLITTERS" -eq 1 ]] && SB_OUT_FLAGS+=( --splitterFile "$SPLITTERS_UNSORTED" )
    [[ "$EMIT_DISCORDANT" -eq 1 ]] && SB_OUT_FLAGS+=( --discordantFile "$DISCORDANT_UNSORTED" )

    set -o pipefail
    bwa mem -t "$THREADS" "${BWA_ARGS[@]}" "$REF" "$R1" "$R2" \
      | samblaster "${SAMBLASTER_ARGS[@]}" "${SB_OUT_FLAGS[@]}" \
      | samtools view -Sb - \
      > "$ALIGNED_UNSORTED"
  else
    set -o pipefail
    bwa mem -t "$THREADS" "${BWA_ARGS[@]}" "$REF" "$R1" "$R2" \
      | samtools view -Sb - \
      > "$ALIGNED_UNSORTED"
  fi
fi

# Sort+index main if needed
if ! bam_complete "$ALIGNED_SORTED"; then
  if [[ -s "$ALIGNED_UNSORTED" ]]; then
    sort_and_index_bam "$ALIGNED_UNSORTED" "$ALIGNED_SORTED"
  else
    die "Neither $ALIGNED_SORTED nor $ALIGNED_UNSORTED exists; cannot proceed."
  fi
fi

# Write audit after aligned BAM exists
if [[ "$AUDIT_COUNTS" -eq 1 ]]; then
  echo "==> Audit: writing $AUDIT_TSV"
  write_audit "$AUDIT_TSV" "$R1_COUNT" "$R2_COUNT" "$ALIGNED_SORTED"
fi

# ----------------------------
# Step 2: samblaster outputs (discordant/splitters) sort+index if requested
# ----------------------------
if [[ "$EMIT_DISCORDANT" -eq 1 ]]; then
  if ! bam_complete "$DISCORDANT_SORTED"; then
    if [[ -s "$DISCORDANT_UNSORTED" ]]; then
      sort_and_index_bam "$DISCORDANT_UNSORTED" "$DISCORDANT_SORTED"
    else
      warn "discordant requested but unsorted missing."
    fi
  else
    echo "==> Exists (skip): $DISCORDANT_SORTED"
  fi
fi

if [[ "$EMIT_SPLITTERS" -eq 1 ]]; then
  if ! bam_complete "$SPLITTERS_SORTED"; then
    if [[ -s "$SPLITTERS_UNSORTED" ]]; then
      sort_and_index_bam "$SPLITTERS_UNSORTED" "$SPLITTERS_SORTED"
    else
      warn "splitters requested but unsorted missing."
    fi
  else
    echo "==> Exists (skip): $SPLITTERS_SORTED"
  fi
fi

# ----------------------------
# Step 3: Derived BAMs from aligned.sorted
# ----------------------------
if [[ "$EMIT_UNMAPPED" -eq 1 ]]; then
  if ! bam_complete "$UNMAPPED_SORTED"; then
    if [[ "$RESUME" -eq 1 && -s "$UNMAPPED_UNSORTED" ]]; then
      sort_and_index_bam "$UNMAPPED_UNSORTED" "$UNMAPPED_SORTED"
    else
      echo "==> Writing unmapped (0x4): $UNMAPPED_SORTED"
      # No MAPQ filter here (unmapped MAPQ is typically 0).
      samtools view -@ "$SORT_THREADS" -b -f 4 "$ALIGNED_SORTED" > "$UNMAPPED_UNSORTED"
      sort_and_index_bam "$UNMAPPED_UNSORTED" "$UNMAPPED_SORTED"
    fi
  else
    echo "==> Exists (skip): $UNMAPPED_SORTED"
  fi
fi

if [[ "$EMIT_SUPPLEMENTARY" -eq 1 ]]; then
  if ! bam_complete "$SUPP_SORTED"; then
    if [[ "$RESUME" -eq 1 && -s "$SUPP_UNSORTED" ]]; then
      sort_and_index_bam "$SUPP_UNSORTED" "$SUPP_SORTED"
    else
      echo "==> Writing supplementary (0x800): $SUPP_SORTED"
      view_mapq_args "$MIN_MAPQ_SUPP"
      samtools view -@ "$SORT_THREADS" -b -f 2048 "${VIEW_MAPQ_ARGS[@]}" "$ALIGNED_SORTED" > "$SUPP_UNSORTED"
      sort_and_index_bam "$SUPP_UNSORTED" "$SUPP_SORTED"
    fi
  else
    echo "==> Exists (skip): $SUPP_SORTED"
  fi
fi

if [[ "$EMIT_PRIMARY" -eq 1 ]]; then
  if ! bam_complete "$PRIMARY_SORTED"; then
    if [[ "$RESUME" -eq 1 && -s "$PRIMARY_UNSORTED" ]]; then
      sort_and_index_bam "$PRIMARY_UNSORTED" "$PRIMARY_SORTED"
    else
      echo "==> Writing primary-only (exclude 0x100+0x800): $PRIMARY_SORTED"
      view_mapq_args "$MIN_MAPQ_PRIMARY"
      samtools view -@ "$SORT_THREADS" -b -F 2304 "${VIEW_MAPQ_ARGS[@]}" "$ALIGNED_SORTED" > "$PRIMARY_UNSORTED"
      sort_and_index_bam "$PRIMARY_UNSORTED" "$PRIMARY_SORTED"
    fi
  else
    echo "==> Exists (skip): $PRIMARY_SORTED"
  fi
fi

if [[ "$EMIT_SECONDARY" -eq 1 ]]; then
  if ! bam_complete "$SECONDARY_SORTED"; then
    if [[ "$RESUME" -eq 1 && -s "$SECONDARY_UNSORTED" ]]; then
      sort_and_index_bam "$SECONDARY_UNSORTED" "$SECONDARY_SORTED"
    else
      echo "==> Writing secondary-only (0x100): $SECONDARY_SORTED"
      view_mapq_args "$MIN_MAPQ_SECONDARY"
      samtools view -@ "$SORT_THREADS" -b -f 256 "${VIEW_MAPQ_ARGS[@]}" "$ALIGNED_SORTED" > "$SECONDARY_UNSORTED"
      sort_and_index_bam "$SECONDARY_UNSORTED" "$SECONDARY_SORTED"
    fi
  else
    echo "==> Exists (skip): $SECONDARY_SORTED"
  fi
fi

# ----------------------------
# Step 4 (optional): Post-hoc MAPQ filtering on samblaster-derived BAMs
# ----------------------------
if [[ "$EMIT_DISCORDANT" -eq 1 && -n "$MIN_MAPQ_DISCORDANT" ]]; then
  if bam_complete "$DISCORDANT_SORTED"; then
    tmpf="$TMPDIR/${SAMPLE}.${OUT_KEY}.discordant.mapqtmp.bam"
    view_mapq_args "$MIN_MAPQ_DISCORDANT"
    echo "==> Filtering discordant by MAPQ >= $MIN_MAPQ_DISCORDANT (post-hoc)"
    samtools view -@ "$SORT_THREADS" -b "${VIEW_MAPQ_ARGS[@]}" "$DISCORDANT_SORTED" > "$tmpf"
    samtools index -@ "$SORT_THREADS" "$tmpf" || true
    mv -f "$tmpf" "$DISCORDANT_SORTED"
    mv -f "${tmpf}.bai" "${DISCORDANT_SORTED}.bai" 2>/dev/null || true
  fi
fi

if [[ "$EMIT_SPLITTERS" -eq 1 && -n "$MIN_MAPQ_SPLITTERS" ]]; then
  if bam_complete "$SPLITTERS_SORTED"; then
    tmpf="$TMPDIR/${SAMPLE}.${OUT_KEY}.splitters.mapqtmp.bam"
    view_mapq_args "$MIN_MAPQ_SPLITTERS"
    echo "==> Filtering splitters by MAPQ >= $MIN_MAPQ_SPLITTERS (post-hoc)"
    samtools view -@ "$SORT_THREADS" -b "${VIEW_MAPQ_ARGS[@]}" "$SPLITTERS_SORTED" > "$tmpf"
    samtools index -@ "$SORT_THREADS" "$tmpf" || true
    mv -f "$tmpf" "$SPLITTERS_SORTED"
    mv -f "${tmpf}.bai" "${SPLITTERS_SORTED}.bai" 2>/dev/null || true
  fi
fi

echo "==> Done."
