# Contamination Detection Pipeline

This repository provides a workflow to identify and quantify contamination in single-cell samples by BLASTing reads against a curated contaminant database and summarizing genus-level proportions.

---

## Setup

### 1. Create Conda Environment

This project requires BLAST+ and Python packages including Biopython and requests.

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate blastenv
```

Alternatively, use pip for Python packages:

```bash
pip install -r requirements.txt
```

---

### 2. Prepare the Contaminant Database

1. Collect FASTA files for contaminants including:
   - All prokaryotes
   - All fungi
   - All protozoans
   - Human genome
   - 
#### The standard database construction:
      Prokaryotes (Bacteria & Archaea)
      Bacteria (RefSeq genomes)
      ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
      (You then download all relevant bacterial genome fasta files listed here)
      
      Archaea (RefSeq genomes)
      ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
      (Same as above for archaea)
      
      Fungi
      Fungi (RefSeq genomes)
      ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
      
      Protozoa (including Plasmodium)
      Protozoa (RefSeq genomes)
      ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt
      
      Or specific Plasmodium genome (example):
      ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.5_ASM276v2/GCF_000002765.5_ASM276v2_genomic.fna.gz
      
      Human Genome
      Human reference genome (GRCh38, primary assembly)
      ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
      
      How to get these
      Download the assembly_summary.txt for each group above.
      
      Parse the summary to get the FTP links to the actual fasta genome files.
      
      Use wget or curl to download those fasta files.
      
      Concatenate all downloaded fasta files into a single contaminants_all.fa.
      
      Build your BLAST DB from that combined fasta.



2. Concatenate all FASTA files into a single file, e.g.:

```bash
cat *.fa > contaminants_all.fa
```

3. Build the BLAST nucleotide database:

```bash
makeblastdb -in contaminants_all.fa -dbtype nucl -out contaminants_db -parse_seqids -taxid_map taxid_map.txt
```

> Note: Generating the `taxid_map.txt` file properly is recommended but optional for now.

---

### 3. Run BLAST and Summarize Genus Proportions

Run the provided bash script `blast.sh` to BLAST each FASTA sample against the contaminant database.
The python code for `summarize_top_genus_proportions.py` is also provided and needs to be in the same 
directory as `blast.sh`

This produces BLAST output files in `blast_nt_out/`. Logs for the the blastn will be made 
(ignore warnings), as well as the raw blast hits, and the summary of genuses hit.

run as:

`./bash.sh`

This script uses Entrez to fetch genus names from accession IDs.

---

## Files

- `blast.sh`: Bash wrapper to run BLAST on samples and keep top hits.
- `summarize_top_genus_proportions.py`: Python script to parse BLAST results and report genus proportions.
- `environment.yml`: Conda environment file.
- `requirements.txt`: Python dependencies.

---

## Notes

- Make sure you have an internet connection for Entrez queries.
- Replace the email in the Python script with your own email for Entrez usage.
- To speed up, you can cache genus mappings locally or build a local taxonomy DB.

---

## Contact

For questions or issues, open an issue or contact the maintainer.
njb8sg@virginia.edu

---

