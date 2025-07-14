# Contamination Detection Pipeline

This repository provides a workflow to identify and quantify contamination in single-cell samples by BLASTing reads against a curated contaminant database and summarizing genus-level proportions.

The code is a reactiation and automated workflow of the steps described in Liu et al.; [Found here](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00889-9#availability-of-data-and-materials)
---

## Setup

### 1. Create Conda Environment

This project requires BLAST+ and Python packages including Biopython and requests.

Create the conda environment with:

```bash
cd~
cd envs/yaml
wget https://raw.githubusercontent.com/a-crow-nowhere/MalariAPI/yaml/blastenv.yml
conda env create -f environment.yml
conda activate blastenv

cd~
mkdir -p /tools/blast/contamdb/out
cd /tools/blast/condtamdb
wget https://raw.githubusercontent.com/a-crow-nowhere/MalariAPI/setup/make_db.sh
./make_db.sh


```

---

### 2. Prepare the Contaminant Database



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

