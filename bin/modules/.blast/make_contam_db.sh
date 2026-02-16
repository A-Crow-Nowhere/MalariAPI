#!/usr/bin/env bash
set -euo pipefail

say(){ echo "[make_contam_db] $*" >&2; }
die(){ say "ERROR: $*"; exit 1; }

usage(){
  cat <<'EOF'
Usage:
  make_contam_db.sh --out-dir DIR [--prefix contaminants_db] [--assembly-summary FILE]

Build curated contaminant DB (minimal-change port of Blast_Reads/setup/make.db.sh)

Outputs in --out-dir:
  accessions.txt
  assembly_summary_refseq.txt (if downloaded)
  *.fna.gz downloads
  all_contaminants.fna(.gz)
  <prefix>.n* (BLAST db files)
  tax_id_map.txt
EOF
}

OUT_DIR=""
PREFIX="contaminants_db"
ASSEMBLY_SUMMARY="assembly_summary_refseq.txt"

[[ $# -eq 0 ]] && { usage; exit 0; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --out-dir) OUT_DIR="${2:-}"; shift 2 ;;
    --prefix) PREFIX="${2:-}"; shift 2 ;;
    --assembly-summary) ASSEMBLY_SUMMARY="${2:-}"; shift 2 ;;
    *) die "Unknown option: $1" ;;
  esac
done

[[ -n "$OUT_DIR" ]] || die "--out-dir is required"

command -v wget >/dev/null 2>&1 || die "wget not found"
command -v makeblastdb >/dev/null 2>&1 || die "makeblastdb not found"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Same accession list as your original
cat <<EOF > accessions.txt
GCF_000001405.39	Human
GCF_000001635.27	Mouse
GCF_000146045.2	Yeast
GCF_000002765.6	Plasmodium
GCF_000006985.1	Mycoplasma
GCF_000012865.1	Mycoplasma
GCF_009858895.2	Ecoli
GCF_000007565.1	Pseudomonas
GCF_000008865.1	Bacillus
GCF_000002655.1	Candida
GCF_000149845.1	Aspergillus
GCF_000002335.4	Trypanosoma
GCF_000002875.2	Leishmania
GCF_000006805.1	Halobacterium
GCF_000010525.1	Methanobrevibacter
GCF_000819615.1	Arabidopsis
GCF_001433935.1	Oryza
GCF_000819615.1	Plant
GCF_000819025.1	PhiX
GCF_000316655.1	pUC19
EOF

if [[ ! -f "$ASSEMBLY_SUMMARY" ]]; then
  say "Downloading RefSeq assembly summary -> $ASSEMBLY_SUMMARY"
  wget -q https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O "$ASSEMBLY_SUMMARY"
fi

while IFS=$'\t' read -r acc label; do
  say "Downloading $acc ($label)"
  match=$(grep "$acc" "$ASSEMBLY_SUMMARY" || true)
  if [[ -z "$match" ]]; then
    say "Accession $acc not found in assembly summary (skipping)"
    continue
  fi

  ftp_path=$(echo "$match" | cut -f20)
  file_name=$(basename "$ftp_path")
  fasta_file="${file_name}_genomic.fna.gz"

  wget -q -c "$ftp_path/$fasta_file" -O "$fasta_file" || {
    say "Failed to download $fasta_file (skipping)"
    continue
  }
done < accessions.txt

say "Concatenating *.fna.gz -> all_contaminants.fna(.gz)"
cat *.fna.gz > all_contaminants.fna.gz
gunzip -c all_contaminants.fna.gz > all_contaminants.fna

say "makeblastdb -> ${PREFIX}"
makeblastdb -in all_contaminants.fna -dbtype nucl -parse_seqids -title "CuratedContamDB" -out "$PREFIX"

awk -F '\t' '{ print $1 "\t" $2 }' accessions.txt > tax_id_map.txt

say "DB prefix: $OUT_DIR/$PREFIX"
