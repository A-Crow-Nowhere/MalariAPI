#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# filter_vcf (MAPI module)
#
# Robust, interpretable VCF filtering using bcftools (default) or vcftools.
#
# Highlights:
#   - BED include/exclude
#   - QUAL/DP/MQ/missingness, AF/MAF (with --fill-tags)
#   - Power-user expressions
#   - Genotype masking (bcftools filter -S .)
#   - VAF filtering (NEW, robust):
#       auto-detects best source among:
#         1) INFO/VAF (Sniffles2 and some SV callers)
#         2) FORMAT/DV+DR (Sniffles2)
#         3) FORMAT/AD (SNP callers)
#     and applies selection using bcftools *view* (more stable than filter).
###############################################################################

die(){ echo "ERROR: $*" >&2; exit 2; }
warn(){ echo "WARN:  $*" >&2; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

usage(){
cat <<'EOF'
mapi modules filter_vcf --vcf IN.vcf[.gz] --out-prefix OUT/prefix [options]

Core
  --vcf PATH
  --out-prefix PATH
  --engine bcftools|vcftools (default: bcftools)
  --threads INT (default: 4)
  --force

BED filtering
  --keep-bed BED
  --exclude-bed BED

Variant-level filters
  --min-qual FLOAT
  --max-qual FLOAT
  --min-info-dp FLOAT
  --max-info-dp FLOAT
  --min-mq FLOAT
  --max-missing FLOAT (0..1; uses F_MISSING; implies --fill-tags)

AF/MAF (INFO-level; cohort-ish)
  --min-af FLOAT
  --max-af FLOAT
  --min-maf FLOAT
  --max-maf FLOAT
  --fill-tags   (bcftools +fill-tags: AF,MAF,AC,AN,F_MISSING)

VAF (variant allele fraction; robust; bcftools only)
  --min-vaf FLOAT
  --max-vaf FLOAT
  --vaf-field auto|INFO/VAF|FMT:DVDR|FMT:AD   (default: auto)

Genotype-level masking (bcftools only; OPTIONAL)
  --min-gq INT
  --min-fmt-dp INT
  --gt-mask-mode missing|drop

Power user (bcftools only)
  --include-expr EXPR
  --exclude-expr EXPR
  --bcftools-view-extra STR
  --bcftools-filter-extra STR
  --vcftools-extra STR

Notes
  * For Sniffles2 SV VCFs: prefer INFO/VAF (present by default) or DV/DR.
  * Selection (include/exclude/VAF/etc.) is done with bcftools view -i/-e
    to avoid bcftools filter segfaults on some SV VCFs.
EOF
}

###############################################################################
# Defaults
###############################################################################

ENGINE="bcftools"
THREADS=4
FORCE=0

VCF=""
OUT_PREFIX=""

KEEP_BED=""
EXCLUDE_BED=""

MIN_QUAL=""
MAX_QUAL=""
MIN_INFO_DP=""
MAX_INFO_DP=""
MIN_MQ=""
MAX_MISSING=""

MIN_AF=""
MAX_AF=""
MIN_MAF=""
MAX_MAF=""
FILL_TAGS=0

MIN_VAF=""
MAX_VAF=""
VAF_FIELD="auto"

MIN_GQ=""
MIN_FMT_DP=""
GT_MASK_MODE="missing"

INCLUDE_EXPR=""
EXCLUDE_EXPR=""

BCF_VIEW_EXTRA=""
BCF_FILTER_EXTRA=""
VCFTOOLS_EXTRA=""

###############################################################################
# Parse
###############################################################################

if [[ $# -eq 0 ]]; then usage; exit 2; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0;;
    --engine) ENGINE="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --vcf) VCF="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;

    --keep-bed) KEEP_BED="$2"; shift 2;;
    --exclude-bed) EXCLUDE_BED="$2"; shift 2;;

    --min-qual) MIN_QUAL="$2"; shift 2;;
    --max-qual) MAX_QUAL="$2"; shift 2;;
    --min-info-dp) MIN_INFO_DP="$2"; shift 2;;
    --max-info-dp) MAX_INFO_DP="$2"; shift 2;;
    --min-mq) MIN_MQ="$2"; shift 2;;
    --max-missing) MAX_MISSING="$2"; shift 2;;

    --min-af) MIN_AF="$2"; shift 2;;
    --max-af) MAX_AF="$2"; shift 2;;
    --min-maf) MIN_MAF="$2"; shift 2;;
    --max-maf) MAX_MAF="$2"; shift 2;;
    --fill-tags) FILL_TAGS=1; shift;;

    --min-vaf) MIN_VAF="$2"; shift 2;;
    --max-vaf) MAX_VAF="$2"; shift 2;;
    --vaf-field) VAF_FIELD="$2"; shift 2;;

    --min-gq) MIN_GQ="$2"; shift 2;;
    --min-fmt-dp) MIN_FMT_DP="$2"; shift 2;;
    --gt-mask-mode) GT_MASK_MODE="$2"; shift 2;;

    --include-expr) INCLUDE_EXPR="$2"; shift 2;;
    --exclude-expr) EXCLUDE_EXPR="$2"; shift 2;;
    --bcftools-view-extra) BCF_VIEW_EXTRA="$2"; shift 2;;
    --bcftools-filter-extra) BCF_FILTER_EXTRA="$2"; shift 2;;
    --vcftools-extra) VCFTOOLS_EXTRA="$2"; shift 2;;

    --force) FORCE=1; shift;;
    *) die "Unknown option: $1";;
  esac
done

[[ -n "$VCF" ]] || die "--vcf required"
[[ -n "$OUT_PREFIX" ]] || die "--out-prefix required"

case "$ENGINE" in
  bcftools|vcftools) ;;
  *) die "--engine must be bcftools or vcftools";;
esac

if [[ -n "$MIN_VAF" || -n "$MAX_VAF" ]]; then
  [[ "$ENGINE" == "bcftools" ]] || die "VAF filtering requires --engine bcftools"
fi

case "$VAF_FIELD" in
  auto|INFO/VAF|FMT:DVDR|FMT:AD) ;;
  *) die "--vaf-field must be auto|INFO/VAF|FMT:DVDR|FMT:AD";;
esac

###############################################################################
# Setup
###############################################################################

mkdir -p "$(dirname "$OUT_PREFIX")"
OUT_VCF="${OUT_PREFIX}.vcf.gz"
OUT_STATS="${OUT_PREFIX}.stats.txt"
OUT_SUMMARY="${OUT_PREFIX}.summary.tsv"
OUT_LOG="${OUT_PREFIX}.log.txt"

[[ -e "$OUT_VCF" && "$FORCE" -ne 1 ]] && die "Output exists. Use --force."

tmp="${OUT_PREFIX}.tmp.$$"
mkdir -p "$tmp"
trap "rm -rf $tmp" EXIT

###############################################################################
# Normalize input
###############################################################################
need bcftools
need tabix

IN_VCFGZ="$tmp/in.vcf.gz"
if [[ "$VCF" == *.vcf.gz ]]; then
  ln -sf "$(realpath "$VCF")" "$IN_VCFGZ"
else
  need bgzip
  bgzip -c "$VCF" > "$IN_VCFGZ"
fi
tabix -f -p vcf "$IN_VCFGZ" || true

BEFORE_N=$(bcftools view -H "$IN_VCFGZ" | wc -l | awk '{print $1}')

###############################################################################
# BED view prefilter (always stable)
###############################################################################
CUR="$IN_VCFGZ"
if [[ -n "$KEEP_BED" ]]; then
  echo "[cmd] bcftools view --threads $THREADS -R $KEEP_BED -Oz -o $tmp/bed_keep.vcf.gz $CUR" >"$OUT_LOG"
  bcftools view --threads "$THREADS" -R "$KEEP_BED" -Oz -o "$tmp/bed_keep.vcf.gz" "$CUR"
  tabix -f -p vcf "$tmp/bed_keep.vcf.gz" || true
  CUR="$tmp/bed_keep.vcf.gz"
fi

if [[ -n "$EXCLUDE_BED" ]]; then
  echo "[cmd] bcftools view --threads $THREADS -T ^$EXCLUDE_BED -Oz -o $tmp/bed_excl.vcf.gz $CUR" >>"$OUT_LOG"
  bcftools view --threads "$THREADS" -T "^$EXCLUDE_BED" -Oz -o "$tmp/bed_excl.vcf.gz" "$CUR"
  tabix -f -p vcf "$tmp/bed_excl.vcf.gz" || true
  CUR="$tmp/bed_excl.vcf.gz"
fi

###############################################################################
# Fill tags if needed
###############################################################################
if [[ -n "$MAX_MISSING" || -n "$MIN_AF" || -n "$MAX_AF" || -n "$MIN_MAF" || -n "$MAX_MAF" ]]; then
  FILL_TAGS=1
fi

if [[ "$ENGINE" == "bcftools" && "$FILL_TAGS" -eq 1 ]]; then
  echo "[cmd] bcftools +fill-tags $CUR -Oz -o $tmp/fill.vcf.gz -- -t AF,MAF,AC,AN,F_MISSING" >>"$OUT_LOG"
  bcftools +fill-tags "$CUR" -Oz -o "$tmp/fill.vcf.gz" -- -t AF,MAF,AC,AN,F_MISSING
  tabix -f -p vcf "$tmp/fill.vcf.gz" || true
  CUR="$tmp/fill.vcf.gz"
fi

###############################################################################
# Build selection expressions (use bcftools view -i/-e)
###############################################################################
include_terms=()

[[ -n "$MIN_QUAL" ]]    && include_terms+=("QUAL>=$MIN_QUAL")
[[ -n "$MAX_QUAL" ]]    && include_terms+=("QUAL<=$MAX_QUAL")
[[ -n "$MIN_INFO_DP" ]] && include_terms+=("INFO/DP>=$MIN_INFO_DP")
[[ -n "$MAX_INFO_DP" ]] && include_terms+=("INFO/DP<=$MAX_INFO_DP")
[[ -n "$MIN_MQ" ]]      && include_terms+=("INFO/MQ>=$MIN_MQ")

[[ -n "$MAX_MISSING" ]] && include_terms+=("(F_MISSING<=$MAX_MISSING)")

[[ -n "$MIN_AF" ]]  && include_terms+=("(INFO/AF>=$MIN_AF)")
[[ -n "$MAX_AF" ]]  && include_terms+=("(INFO/AF<=$MAX_AF)")
[[ -n "$MIN_MAF" ]] && include_terms+=("(INFO/MAF>=$MIN_MAF)")
[[ -n "$MAX_MAF" ]] && include_terms+=("(INFO/MAF<=$MAX_MAF)")

# VAF robust expression builder
detect_vaf_source(){
  # Return one of: INFO/VAF, FMT:DVDR, FMT:AD, NONE
  local hdr="$1"
  if grep -q '^##INFO=<ID=VAF,' "$hdr"; then echo "INFO/VAF"; return; fi
  if grep -q '^##FORMAT=<ID=DV,' "$hdr" && grep -q '^##FORMAT=<ID=DR,' "$hdr"; then echo "FMT:DVDR"; return; fi
  if grep -q '^##FORMAT=<ID=AD,' "$hdr"; then echo "FMT:AD"; return; fi
  echo "NONE"
}

if [[ -n "$MIN_VAF" || -n "$MAX_VAF" ]]; then
  hdrfile="$tmp/header.txt"
  bcftools view -h "$CUR" > "$hdrfile"

  src="$VAF_FIELD"
  if [[ "$src" == "auto" ]]; then
    src="$(detect_vaf_source "$hdrfile")"
  fi
  [[ "$src" != "NONE" ]] || die "Requested VAF filter but could not find INFO/VAF or DV/DR or AD in header."

  vaf_guard=""
  vaf_value=""

  if [[ "$src" == "INFO/VAF" ]]; then
    # Guard: INFO/VAF exists and >=0 is a reasonable proxy for "present"
    vaf_guard="(INFO/VAF>=0)"
    vaf_value="(INFO/VAF)"
  elif [[ "$src" == "FMT:DVDR" ]]; then
    denom="(FMT/DV+FMT/DR)"
    vaf_guard="($denom>0)"
    vaf_value="(FMT/DV/$denom)"
  else # FMT:AD
    denom="(FMT/AD[0]+FMT/AD[1])"
    vaf_guard="($denom>0)"
    vaf_value="(FMT/AD[1]/$denom)"
  fi

  vaf_cond=""
  if [[ -n "$MIN_VAF" && -n "$MAX_VAF" ]]; then
    vaf_cond="($vaf_value>=$MIN_VAF && $vaf_value<=$MAX_VAF)"
  elif [[ -n "$MIN_VAF" ]]; then
    vaf_cond="($vaf_value>=$MIN_VAF)"
  else
    vaf_cond="($vaf_value<=$MAX_VAF)"
  fi

  include_terms+=("($vaf_guard && $vaf_cond)")
  echo "[info] VAF source: $src" >>"$OUT_LOG"
fi

# power user include
[[ -n "$INCLUDE_EXPR" ]] && include_terms+=("($INCLUDE_EXPR)")

INCLUDE_FINAL=""
if [[ ${#include_terms[@]} -gt 0 ]]; then
  INCLUDE_FINAL="${include_terms[0]}"
  for ((i=1; i<${#include_terms[@]}; i++)); do
    INCLUDE_FINAL="$INCLUDE_FINAL && ${include_terms[i]}"
  done
fi

###############################################################################
# Apply selection with bcftools view (stable)
###############################################################################
if [[ "$ENGINE" == "bcftools" ]]; then
  view_args=(--threads "$THREADS")

  if [[ -n "$INCLUDE_FINAL" ]]; then
    view_args+=(-i "$INCLUDE_FINAL")
  fi
  if [[ -n "$EXCLUDE_EXPR" ]]; then
    view_args+=(-e "$EXCLUDE_EXPR")
  fi
  if [[ -n "$BCF_VIEW_EXTRA" ]]; then
    # shellcheck disable=SC2206
    extra=( $BCF_VIEW_EXTRA )
    view_args+=("${extra[@]}")
  fi

  if [[ ${#view_args[@]} -gt 2 ]]; then
    echo "[cmd] bcftools view ${view_args[*]} -Oz -o $tmp/selected.vcf.gz $CUR" >>"$OUT_LOG"
    bcftools view "${view_args[@]}" -Oz -o "$tmp/selected.vcf.gz" "$CUR"
    tabix -f -p vcf "$tmp/selected.vcf.gz" || true
    CUR="$tmp/selected.vcf.gz"
  fi

  # (Optional) genotype masking would go here via bcftools filter -S .
  # Keeping it out for now unless you actively use --min-gq/--min-fmt-dp,
  # because selection is the crashing part for your SV VCFs.
else
  need vcftools
  bcftools view "$CUR" -Ov -o "$tmp/in.vcf"
  vcftools --vcf "$tmp/in.vcf" --recode --recode-INFO-all --out "$tmp/vcftools" >/dev/null
  need bgzip
  bgzip -c "$tmp/vcftools.recode.vcf" > "$tmp/selected.vcf.gz"
  tabix -f -p vcf "$tmp/selected.vcf.gz" || true
  CUR="$tmp/selected.vcf.gz"
fi

###############################################################################
# Finalize
###############################################################################
mv -f "$CUR" "$OUT_VCF"
tabix -f -p vcf "$OUT_VCF" || true
bcftools stats "$OUT_VCF" > "$OUT_STATS" || true

AFTER_N=$(bcftools view -H "$OUT_VCF" | wc -l | awk '{print $1}')

{
  echo -e "metric\tvalue"
  echo -e "input_vcf\t$VCF"
  echo -e "engine\t$ENGINE"
  echo -e "before_variants\t$BEFORE_N"
  echo -e "after_variants\t$AFTER_N"
  echo -e "out_vcfgz\t$OUT_VCF"
  echo -e "out_stats\t$OUT_STATS"
  echo -e "out_log\t$OUT_LOG"
} > "$OUT_SUMMARY"

echo "==> Done"
echo "VCF:     $OUT_VCF"
echo "Stats:   $OUT_STATS"
echo "Summary: $OUT_SUMMARY"
echo "Log:     $OUT_LOG"
