# Contamination Detection Pipeline

This repository provides a workflow to identify and quantify contamination in single-cell samples by BLASTing reads against a curated contaminant database and summarizing genus-level proportions.

The code is a reactiation and automated workflow of the steps described in Liu et al.; [Found here](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00889-9#availability-of-data-and-materials)
---

## Copy and paste chunk by chunk into bash:

### Create environment and construct database (skip if done prevoiusly, takes ~5-10 minutes). Further, consider adding/subtracting the custom database genomes depeneding on what you need. Those can be found in the make_db.sh file. By default: 

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
```bash
cd~
cd envs/yaml
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/yaml/blastenv.yml
conda env create -f environment.yml
conda activate blastenv

cd~
mkdir -p /tools/blast/contamdb/out
cd /tools/blast/condtamdb
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/setup/make_db.sh
./make_db.sh
```

### Create in and out directories. Fastas are analized in this exact directory and must end with .fasta
```bash
cd~
mkdir -p /tools/blast/fastas #MOVE FASTA FILE(S) HERE
mkdir -p /tools/blast/blastOut
```

### Runs the code - the larger the fasta, the longer it takes. log files are produced mid run for you to check progress
```bash
cd /tools/blast
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/blast/blast.sh
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/blast/summarize_top_genus_proportions.py
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/blast/summarize.sh

chmod +x blast.sh
chmod +x summarize.sh

./blast.sh
./summarize.sh
```

### Raw blast outputs will be in the /tools/blast/blastOut dir
### summarized outputs will be in the /tools/blast dir
---


