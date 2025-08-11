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

#Call DG minima in 2000bp range around teh 5'&3' region of a bed region in sliding 50bp windows.
#The peak closest to the input breakpoints are retained. Also calculates GC content in a 100bp around that breakpoint.
./breakpoints.sh <regionstoanalyze.bed> <reference.fasta> => final_summary.tsv

#Determine if peak position called in final_summary.tsv is in an intergenic region based on CDS regions in a .gff file. 
./intergenic.sh <final_summary.clean.tsv> <reference.gff> => finall_summary.annotated.tsv




```
Inputs for breakponts.sh need to be a normal bed file in the format: chr, start, end (no header)
and a reference genome fasta file with an index in the same directory.

Inputs for intergenic.sh are the output of breakpoints.sh, and a .gff style file. 

Outputs:
\[temporary directory\] located in the local tmp dir, is kept temporarily, and retains intermediate files including viennaRNA raw outputs.
final_summary_raw.tsv (raw, unorganized dg peak calling)
final_summary.tsv (cleaned version)
final_summary.annotated.tsv (calls if the determined breakpoins are in an intergenic region)



