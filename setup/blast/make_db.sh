#!/usr/bin/env bash

# contam_db_setup.sh
# Downloads curated reference genomes for cell culture contamination screening
# and builds a BLAST database
set -euo pipefail
set -x
### 1. OUTPUT DIRECTORY
DB_DIR="$HOME/tools/blast/contamdb/out"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

### 2. GENOME LIST (accession <tab> label)
cat <<EOF > accessions.txt
GCF_000001405.39        Human
GCF_000001635.27        Mouse
GCF_000146045.2 Yeast
GCF_000002765.6 Plasmodium
GCF_000006985.1 Mycoplasma
GCF_000012865.1 Mycoplasma
GCF_009858895.2 Ecoli
GCF_000007565.1 Pseudomonas
GCF_000008865.1 Bacillus
GCF_000002655.1 Candida
GCF_000149845.1 Aspergillus
GCF_000002335.4 Trypanosoma
GCF_000002875.2 Leishmania
GCF_000006805.1 Halobacterium
GCF_000010525.1 Methanobrevibacter
GCF_000819615.1 Arabidopsis
GCF_001433935.1 Oryza
GCF_000819615.1 Plant
GCF_000819025.1 PhiX
GCF_000316655.1 pUC19
EOF

# Download assembly summary if not present
ASSEMBLY_SUMMARY="assembly_summary_refseq.txt"
if [ ! -f "$ASSEMBLY_SUMMARY" ]; then
  echo "📄 Downloading RefSeq assembly summary..."
  wget -q https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O "$ASSEMBLY_SUMMARY"
fi

# Read each accession and label
while IFS=$'\t' read -r acc label; do
  echo "📥 Downloading $acc ($label)"

  # Find line in assembly_summary with matching accession
  match=$(grep "$acc" "$ASSEMBLY_SUMMARY" || true)
  if [ -z "$match" ]; then
    echo "❌ Accession $acc not found in assembly summary"
    continue
  fi

  # Extract FTP path and construct filename
  ftp_path=$(echo "$match" | cut -f20)
  file_name=$(basename "$ftp_path")
  fasta_file="${file_name}_genomic.fna.gz"

  # Download
  wget -q -c "$ftp_path/$fasta_file" -O "$fasta_file" || {
    echo "❌ Failed to download $fasta_file"
    continue
  }

done < accessions.txt


### 4. CONCATENATE FASTA FILES
cat *.fna.gz > all_contaminants.fna.gz
gunzip -c all_contaminants.fna.gz > all_contaminants.fna

### 5. MAKE BLAST DATABASE
makeblastdb -in all_contaminants.fna \
  -dbtype nucl \
  -parse_seqids \
  -title "CuratedContamDB" \
  -out contaminants_db

### 6. CREATE CATEGORY MAP FILE
awk -F '\t' '{ print $1 "\t" $2 }' accessions.txt > tax_id_map.txt

### ✅ DONE
cat <<EOF

✅ Curated contamination BLAST database created at: $DB_DIR
- Database prefix: contaminants_db
- Accession-to-organism mapping: tax_id_map.txt
- Combined FASTA: all_contaminants.fna

Use with:
blastn -query sample.fasta -db $DB_DIR/contaminants_db [...]
EOF
