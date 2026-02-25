#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# sniffles2_call: Sniffles2 SV calling (+ optional CNV BED from DEL/DUP)
#
# Outputs:
#   - {sample}.sniffles2.vcf.gz (+ .tbi)
#   - {sample}.sniffles2.snf              (optional; for multisample workflows)
#   - {sample}.sniffles2.cnv.bed          (DEL/DUP extracted as intervals)
#
# Notes:
#   - Sniffles2 supports BAM/CRAM input. For DEL sequence output, provide --reference.
#   - Mosaic mode: --mosaic
#   - Force-genotype a known SV VCF: --genotype-vcf <vcf>
###############################################################################

# ----------------------------
# Defaults
# ----------------------------
BAM=""
OUT_DIR=""
SAMPLE=""
THREADS=4
REFERENCE=""
TANDEM_REPEATS=""
MOSAIC=0
MINSUPPORT=""
EMIT_SNF=0
GENOTYPE_VCF=""
EXTRA_ARGS=""

# ----------------------------
# Help
# ----------------------------
usage() {
  cat <<'EOF'
Usage:
  mapi modules sniffles2_call --bam <bam|cram> --out-dir <dir> --sample <name> [options]

Required:
  --bam FILE            Input long-read alignment (.bam or .cram)
  --out-dir DIR         Output directory
  --sample NAME         Sample name/prefix for outputs

Options:
  --threads INT         Threads (default: 4)
  --reference FASTA     Reference FASTA (recommended; required for DEL sequence output;
                        also strongly recommended for CRAM)
  --tandem-repeats BED  Tandem repeat annotations BED (improves calling in repeats)
  --mosaic              Enable mosaic/non-germline mode
  --minsupport INT      Minimum read support (passed to sniffles as --minsupport)
  --emit-snf            Also emit .snf (for population/multisample mode)
  --genotype-vcf VCF    Force-genotype a known SV set (passed as --genotype-vcf)
  --extra "ARGS"        Extra args passed directly to sniffles (quoted)

Outputs (in --out-dir):
  <sample>.sniffles2.vcf.gz
  <sample>.sniffles2.vcf.gz.tbi
  <sample>.sniffles2.cnv.bed
  <sample>.sniffles2.snf          (if --emit-snf)

Examples:
  mapi modules sniffles2_call \
    --bam reads.minimap2.bam \
    --out-dir /scratch/$USER/MalariAPI/scratch/test/sniffles \
    --sample C01_1 \
    --threads 8 \
    --reference ref.fa \
    --tandem-repeats repeats.bed

  # Mosaic mode:
  mapi modules sniffles2_call --bam in.bam --out-dir out --sample S1 --mosaic

  # Force genotyping:
  mapi modules sniffles2_call --bam in.bam --out-dir out --sample S1 \
    --genotype-vcf known_svs.vcf.gz
EOF
}

die(){ echo "ERROR: $*" >&2; exit 2; }

# ----------------------------
# Arg parse
# ----------------------------
if [[ $# -eq 0 ]]; then usage; exit 0; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --bam) BAM="${2:-}"; shift 2 ;;
    --out-dir) OUT_DIR="${2:-}"; shift 2 ;;
    --sample) SAMPLE="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --reference) REFERENCE="${2:-}"; shift 2 ;;
    --tandem-repeats) TANDEM_REPEATS="${2:-}"; shift 2 ;;
    --mosaic) MOSAIC=1; shift 1 ;;
    --minsupport) MINSUPPORT="${2:-}"; shift 2 ;;
    --emit-snf) EMIT_SNF=1; shift 1 ;;
    --genotype-vcf) GENOTYPE_VCF="${2:-}"; shift 2 ;;
    --extra) EXTRA_ARGS="${2:-}"; shift 2 ;;
    *) die "Unknown option: $1 (see --help)" ;;
  esac
done

# ----------------------------
# Validate
# ----------------------------
[[ -n "$BAM" ]] || die "--bam is required"
[[ -n "$OUT_DIR" ]] || die "--out-dir is required"
[[ -n "$SAMPLE" ]] || die "--sample is required"
[[ -s "$BAM" ]] || die "Input not found or empty: $BAM"

mkdir -p "$OUT_DIR"

# If CRAM, reference is generally needed by htslib/pysam
if [[ "$BAM" == *.cram && -z "$REFERENCE" ]]; then
  echo "WARNING: input is CRAM but --reference not provided; Sniffles/htslib may fail." >&2
fi

if [[ -n "$REFERENCE" && ! -s "$REFERENCE" ]]; then
  die "Reference FASTA not found: $REFERENCE"
fi
if [[ -n "$TANDEM_REPEATS" && ! -s "$TANDEM_REPEATS" ]]; then
  die "Tandem repeats BED not found: $TANDEM_REPEATS"
fi
if [[ -n "$GENOTYPE_VCF" && ! -s "$GENOTYPE_VCF" ]]; then
  die "Genotype VCF not found: $GENOTYPE_VCF"
fi

# ----------------------------
# Outputs
# ----------------------------
VCF_GZ="$OUT_DIR/${SAMPLE}.sniffles2.vcf.gz"
CNV_BED="$OUT_DIR/${SAMPLE}.sniffles2.cnv.bed"
SNF="$OUT_DIR/${SAMPLE}.sniffles2.snf"

# ----------------------------
# Run Sniffles2
# ----------------------------
cmd=(sniffles
  --input "$BAM"
  --vcf "$VCF_GZ"
  --threads "$THREADS"
)

# Reference helps for DEL sequences and often for CRAM workflows
if [[ -n "$REFERENCE" ]]; then
  cmd+=(--reference "$REFERENCE")
fi

if [[ -n "$TANDEM_REPEATS" ]]; then
  cmd+=(--tandem-repeats "$TANDEM_REPEATS")
fi

if [[ "$MOSAIC" -eq 1 ]]; then
  cmd+=(--mosaic)
fi

if [[ -n "$MINSUPPORT" ]]; then
  cmd+=(--minsupport "$MINSUPPORT")
fi

if [[ -n "$GENOTYPE_VCF" ]]; then
  cmd+=(--genotype-vcf "$GENOTYPE_VCF")
fi

if [[ "$EMIT_SNF" -eq 1 ]]; then
  cmd+=(--snf "$SNF")
fi

# Extra args (split on spaces intentionally)
if [[ -n "$EXTRA_ARGS" ]]; then
  # shellcheck disable=SC2206
  extra_arr=($EXTRA_ARGS)
  cmd+=("${extra_arr[@]}")
fi

echo "[sniffles2_call] Running:"
printf '  %q' "${cmd[@]}"; echo
"${cmd[@]}"

# Ensure index exists (Sniffles2 supports .vcf.gz output)
if [[ -s "$VCF_GZ" ]]; then
  if [[ ! -s "${VCF_GZ}.tbi" ]]; then
    echo "[sniffles2_call] Indexing VCF..."
    tabix -f -p vcf "$VCF_GZ"
  fi
else
  die "Sniffles did not produce VCF: $VCF_GZ"
fi

# ----------------------------
# Derive CNV BED from DEL/DUP records
# ----------------------------
# CNV-ish track: treat DEL/DUP SVs as CNV intervals (common convention).
# We use bcftools query: CHROM POS INFO/END INFO/SVTYPE INFO/SVLEN
echo "[sniffles2_call] Writing CNV BED (DEL/DUP) -> $CNV_BED"
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$VCF_GZ" \
  | awk 'BEGIN{OFS="\t"} $4=="DEL" || $4=="DUP" { 
      # BED is 0-based start; VCF POS is 1-based
      s=$2-1; e=$3;
      if (s<0) s=0;
      print $1, s, e, $4, $5
    }' \
  > "$CNV_BED"

echo "[sniffles2_call] Done."
echo "  VCF : $VCF_GZ"
echo "  CNV : $CNV_BED"
if [[ "$EMIT_SNF" -eq 1 ]]; then
  echo "  SNF : $SNF"
fi
