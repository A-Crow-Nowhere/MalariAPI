#!/bin/bash
set -euo pipefail

FASTA_FILE="accession.fa"     # Your input FASTA file
OUTPUT_FILE="taxid_map.txt"   # Output file for accession-to-taxid map

echo "Generating taxid map from $FASTA_FILE..."

# Extract unique accessions from FASTA headers (assumes first word after '>' is accession)
grep '^>' "$FASTA_FILE" | sed 's/^>//' | awk '{print $1}' | sort -u > accessions.txt

# Empty output file or create if doesn't exist
> "$OUTPUT_FILE"

# Function to get taxid for accession
get_taxid() {
    local acc=$1
    taxid=$(esearch -db nucleotide -query "$acc" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId 2>/dev/null || echo "0")
    echo -e "${acc}\t${taxid}"
}

# Iterate over each accession and write accession + taxid to output file
while read -r acc; do
    taxid_line=$(get_taxid "$acc")
    echo "$taxid_line" >> "$OUTPUT_FILE"
done < accessions.txt

echo "Taxid map saved to $OUTPUT_FILE"
