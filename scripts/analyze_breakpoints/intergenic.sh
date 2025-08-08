#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <dg_peaks.tsv> <annotation.gff>" >&2
  exit 1
fi

DG_PEAKS="$1"
GFF="$2"
OUTPUT="${DG_PEAKS%.tsv}.with_region.tsv"

# Extract 'gene' features only, convert to BED (0-based, half-open)
awk -F'\t' '$0 !~ /^#/ && $3=="gene" {print $1, $4-1, $5}' OFS='\t' "$GFF" > genes.bed

# Convert DG peaks to BED format (chrom, start=pos-1, end=pos, rest of fields)
awk -F'\t' -v OFS='\t' 'NR>1 {
    chrom = $3;
    pos = $8;
    gsub(/[^0-9]/,"",pos);
    start = pos - 1;
    if (start < 0) start = 0;
    end = pos;
    print chrom, start, end, $0
}' "$DG_PEAKS" > dg_peaks.bed.tmp

# Find overlaps with genes
bedtools intersect -a dg_peaks.bed.tmp -b genes.bed -wa -u > peaks_in_genes.bed

# Extract flank_name (column 2 in DG peaks = field 6 in dg_peaks.bed.tmp)
awk -F'\t' '{print $6}' peaks_in_genes.bed | sort | uniq > genic_flanks.txt

# Annotate DG peaks with genic/intergenic status
{
  read header
  echo -e "${header}\tregion_type"
  while IFS= read -r line; do
    flank=$(echo "$line" | cut -f2)
    if grep -qxF "$flank" genic_flanks.txt; then
      echo -e "${line}\tgenic"
    else
      echo -e "${line}\tintergenic"
    fi
  done
} < "$DG_PEAKS" > "$OUTPUT"

# Clean up temporary files
rm genes.bed dg_peaks.bed.tmp peaks_in_genes.bed genic_flanks.txt

echo "✅ Annotation with gene features done. Output saved to $OUTPUT"
