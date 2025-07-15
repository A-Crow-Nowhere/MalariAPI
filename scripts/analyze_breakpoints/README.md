```bash
cd~
cd envs/yaml
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/yaml/breakpoints.yml
conda env create -f breakpoints.yml
conda activate breakpointenv

cd~
mkdir -p tools/breakpoints/
cd tools/breakpoints/
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/breakpoints.sh
chmod +x breakpoints.sh

./breakpoints.sh <regions.bed> <reference.fasta>
```
Inputs for breakponts.sh need to be a normal bed file in the format: chr, start, end (no header)
and a reference genome fasta file with an index in the same directory.

Outputs:
dg_results.tsv, with columns:
SeqID	Start	WindowSeq	RNA_DG	Approx_DNA_DG

dg_peaks.tsv, with columns:
SeqID	PeakPos	PeakDG	Status	MeanDG	SD_DG	CI_Lower	CI_Upper	GC100bp

