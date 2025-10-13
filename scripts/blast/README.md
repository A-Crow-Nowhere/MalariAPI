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
conda env create -f blastenv.yml
conda activate blastenv

cd~
mkdir -p /bin/blast_tool/bin
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/blast/blast.sh
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/setup/blast/make_db.sh
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/blast/summarize_top_genus_proportions.py
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/blast/summarize.sh

chmod +x summarize_top_genus_proportions.py
chmod +x blast.sh
chmod +x summarize.sh
chmod +x make_db.sh

mkdir -p your/path/fastas #MOVE FASTA FILE(S) HERE
mkdir -p your/path/blastOut

./make_db.sh #will take a while to construct customDB

./blast.sh
./summarize.sh

### Raw blast outputs and log files will be in the /your/path/blastOut dir
### summarized outputs will be in the the working dir
---


