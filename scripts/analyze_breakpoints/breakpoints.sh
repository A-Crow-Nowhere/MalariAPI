#!/usr/bin/env bash
set -euo pipefail

# Activate conda environment (make sure conda is initialized in your shell)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate breakpointenv   # change to your actual env if needed

BED="$1"
FASTA="$2"
WINDOW_SIZE=50
FLANK=1000
DG_THRESHOLD=-5.8

OUTPUT="${BED%.bed}.dg_peaks.tsv"

# Write header with region_id column
echo -e "region_id\tflank_name\tchromosome\tpeak_call\tstatus\tdg_of_peak\tgc_content_around_peak\tabsolute_peak_position" > "$OUTPUT"

# Load chromosome lengths from FASTA index (.fai)
declare -A chrom_lengths
if [[ ! -f "${FASTA}.fai" ]]; then
    echo "[!] FASTA index file (${FASTA}.fai) not found. Creating it now..." >&2
    samtools faidx "$FASTA"
fi
while read -r chr len _; do
    chrom_lengths["$chr"]=$len
done < "${FASTA}.fai"

# Function to calculate DG using RNAfold
calc_dg() {
    local seq="$1"
    echo "$seq" | RNAfold --noPS 2>/dev/null | awk 'NR==2 {print $NF}' | tr -d '()'
}

# Function to calculate GC content
calc_gc() {
    local seq="$1"
    echo "$seq" | awk '{gc=gsub(/[GCgc]/,""); printf "%.2f", (gc/length($0))*100}'
}

# Function to find best DG peak in the flank region
find_peak() {
    local chrom="$1" start="$2" end="$3" flank_name="$4" orig_coord="$5"

    # Extract flanking sequence from FASTA using samtools faidx
    seq=$(samtools faidx "$FASTA" "${chrom}:${start}-${end}" | tail -n +2 | tr -d '\n')

    best_dg=999
    best_pos="NA"
    best_gc="NA"
    found_status="no_peak"
    best_dist=999999999

    seq_len=${#seq}
    for ((i=0; i<=seq_len-WINDOW_SIZE; i++)); do
        window_seq=${seq:i:WINDOW_SIZE}
        dg=$(calc_dg "$window_seq")

        # Check if DG is a valid number and below threshold
        if [[ "$dg" != "nan" && "$dg" != "" ]] && (( $(echo "$dg < $DG_THRESHOLD" | bc -l) )); then
            # Absolute genomic coordinate of this window's center
            center_offset=$((i + WINDOW_SIZE/2))
            abs_pos=$((start + center_offset))
            dist=$(( abs_pos > orig_coord ? abs_pos - orig_coord : orig_coord - abs_pos ))

            # Pick peak closest to original coordinate
            if [[ "$best_pos" == "NA" ]] || (( dist < best_dist )); then
                best_dg="$dg"
                best_pos="$abs_pos"
                best_gc=$(calc_gc "$window_seq")
                best_dist="$dist"
                found_status="peak_found"
            fi
        fi
    done

    # Extract region_id by stripping trailing _5prime or _3prime from flank_name
    region_id="${flank_name%_*}"

    echo -e "${region_id}\t${flank_name}\t${chrom}\t${best_pos}\t${found_status}\t${best_dg}\t${best_gc}\t${best_pos}" >> "$OUTPUT"
}

# Main loop over BED file; skip header line if present
tail -n +2 "$BED" | while IFS=$'\t' read -r chrom start end; do
    # Trim whitespace
    chrom=$(echo "$chrom" | tr -d '[:space:]')
    start=$(echo "$start" | tr -d '[:space:]')
    end=$(echo "$end" | tr -d '[:space:]')

    if [[ -z "$chrom" || -z "$start" || -z "$end" ]]; then
        echo "Skipping malformed or empty line" >&2
        continue
    fi

    # Get chromosome length; if not found, skip or warn
    chr_len=${chrom_lengths["$chrom"]:-0}
    if (( chr_len == 0 )); then
        echo "Warning: chromosome length not found for $chrom, skipping" >&2
        continue
    fi

    echo "DEBUG: chrom=$chrom, start=$start, end=$end, chr_len=$chr_len" >&2

    # 5' flank region (start - FLANK to start + FLANK)
    region_start=$((start - FLANK))
    (( region_start < 1 )) && region_start=1
    region_end=$((start + FLANK))
    (( region_end > chr_len )) && region_end=$chr_len
    echo "DEBUG 5prime region: $region_start - $region_end" >&2
    find_peak "$chrom" "$region_start" "$region_end" "${chrom}:${start}-${end}_5prime" "$start"

    # 3' flank region (end - FLANK to end + FLANK)
    region_start=$((end - FLANK))
    (( region_start < 1 )) && region_start=1
    region_end=$((end + FLANK))
    (( region_end > chr_len )) && region_end=$chr_len
    echo "DEBUG 3prime region: $region_start - $region_end" >&2
    find_peak "$chrom" "$region_start" "$region_end" "${chrom}:${start}-${end}_3prime" "$end"

done

echo "✅ DG peak search completed. Results saved to $OUTPUT"
