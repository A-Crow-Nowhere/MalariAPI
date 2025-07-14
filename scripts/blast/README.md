# Contamination Detection Pipeline

This repository provides a workflow to identify and quantify contamination in single-cell samples by BLASTing reads against a curated contaminant database and summarizing genus-level proportions.

The code is a reactiation and automated workflow of the steps described in Liu et al.; [Found here](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00889-9#availability-of-data-and-materials)
---

##Copy and paste chunk by chunk into bash:

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

```bash
cd~
mkdir -p /tools/blast/fastas #MOVE FASTA FILE(S) HERE
mkdir -p /tools/blast/blastOut
```

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

Raw blast outputs will be in the /tools/blast/blastOut dir
summarized outputs will be in the /tools/blast dir
---


