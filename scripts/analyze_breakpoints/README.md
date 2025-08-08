```bash
cd~
cd envs/yaml
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/yaml/breakpoints.yml
conda env create -f breakpoints.yml
conda activate breakpointenv

cd~
mkdir -p tools/breakpoints/
cd tools/breakpoints/
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/analyze_breakpoints/breakpoints.sh
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/analyze_breakpoints/intergenic.sh
chmod +x breakpoints.sh

./breakpoints.sh <regions.bed> <reference.fasta> => final_summary.tsv
./intergenic.sh <final_summary.clean.tsv> 




```
Inputs for breakponts.sh need to be a normal bed file in the format: chr, start, end (no header)
and a reference genome fasta file with an index in the same directory.

Outputs:
dg_peaks.tsv (dg peak calling)
dg_peaks.with_regions.tsv (calls if the determined breakpoins are in an intergenic region). 



