#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <final_summary.tsv> <annotation.gff>" >&2
  exit 1
fi

PEAKS="$1"
GFF="$2"
OUTPUT="${PEAKS%.tsv}.annotated.tsv"

echo "[INFO] Extracting CDS regions from GFF and converting to BED..."
awk -F'\t' 'BEGIN {OFS="\t"}
  $0 !~ /^#/ && $3 == "CDS" {
    start = $4 - 1
    if (start < 0) start = 0
    print $1, start, $5
  }
' "$GFF" | sort -k1,1 -k2,2n > cds.bed

echo "[DEBUG] Preview of cds.bed (first 5 lines):"
head -n 5 cds.bed | cat -t

echo "[INFO] Preparing peaks BED for bedtools intersect..."
awk -F'\t' 'BEGIN {OFS="\t"}
  NR > 1 {
    chrom = $2
    pos = $7
    start = pos - 1
    if (start < 0) start = 0
    end = pos
    print chrom, start, end
  }
' "$PEAKS" | sort -k1,1 -k2,2n > peaks.bed

echo "[DEBUG] Preview of peaks.bed (first 5 lines):"
head -n 5 peaks.bed | cat -t

echo "[DEBUG] File existence and type checks:"
ls -l cds.bed peaks.bed
file cds.bed peaks.bed

echo "[INFO] Running bedtools intersect to find peaks overlapping CDS regions..."
bedtools intersect -a peaks.bed -b cds.bed -wa -u > peaks_in_cds.bed

echo "[DEBUG] peaks_in_cds.bed preview (first 5 lines):"
head -n 5 peaks_in_cds.bed | cat -t

echo "[INFO] Generating genic peaks coordinates list..."
awk '{print $1, $2, $3}' peaks_in_cds.bed | tr ' ' '\t' | sort -k1,1 -k2,2n > genic_coords.tsv

echo "[INFO] Annotating original peaks using fast awk join..."

awk -F'\t' '
  NR==FNR { # genic coords file
    key = $1 FS $2 FS $3
    genic[key] = 1
    next
  }
  NR!=FNR { # peaks file
    if (FNR == 1) {
      print $0 "\tregion_type"
      next
    }
    chrom = $2
    abs_pos = $7
    start = abs_pos - 1
    key = chrom FS start FS abs_pos
    region = (key in genic) ? "genic" : "intergenic"
    print $0 "\t" region
  }
' genic_coords.tsv "$PEAKS" > "$OUTPUT"

echo "[INFO] Cleaning up temporary files..."
rm cds.bed peaks.bed peaks_in_cds.bed genic_coords.tsv

echo "✅ Annotation complete. Output saved to $OUTPUT"
