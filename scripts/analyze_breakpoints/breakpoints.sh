#!/bin/bash

# call_dg_peaks.sh - Extract sequences from BED, calculate ΔG in 50bp windows, call peaks closest to center,
# and output cleaned summary with absolute genomic peak positions

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.bed> <reference.fa>" >&2
  exit 1
fi

BED_IN="$1"
FASTA="$2"
SCRIPT_DIR=$(dirname "$(realpath "$0")")
TMPDIR=$(mktemp -d -t dgpeak_XXXXXX)

echo "[+] Temp dir: $TMPDIR"
trap 'echo "[!] Error occurred. Temp files preserved in $TMPDIR"; exit 1' ERR

FASTA_OUT="$TMPDIR/flanks.fa"
DG_TSV="$TMPDIR/dg_windows.tsv"
GC_TSV="$TMPDIR/gc.tsv"
PEAKS_TSV="$TMPDIR/peaks.tsv"
BED_FLANKS="$TMPDIR/flanks.bed"
SEQ_TSV="$TMPDIR/seq.tsv"
PARAMFILE="./dna_mathews2004.par"

SUMMARY_RAW="$PWD/final_summary_raw.tsv"
SUMMARY_OUT="$PWD/final_summary.tsv"

# Extract 2000bp regions centered on 5' and 3' ends (±1000bp)
sed 's/\r//' "$BED_IN" | awk -v OFS="\t" 'BEGIN{flank=1000} {
  if($2-flank >= 0) print $1, $2-flank, $2+flank, $1":"$2"-"$3"_flank_5";
  print $1, $3-flank, $3+flank, $1":"$2"-"$3"_flank_3";
}' > "$BED_FLANKS"

bedtools getfasta -fi "$FASTA" -bed "$BED_FLANKS" -name -fo "$FASTA_OUT"

# Clean FASTA output and extract only the short flank_name
awk 'BEGIN{RS=">"; ORS=""} NR>1 {
  n=split($0, a, "\n"); header=a[1]; sub(/::.*/, "", header);
  seq=""; for(i=2;i<=n;i++) seq=seq a[i];
  gsub("\r", "", seq);
  print header"\t"seq"\n"
}' "$FASTA_OUT" > "$SEQ_TSV"

# Calculate ΔG
: > "$DG_TSV"
while IFS=$'\t' read -r region seq; do
  seqlen=${#seq}
  for ((i=0; i<=seqlen-50; i++)); do
    window=${seq:i:50}
    rna=${window//T/U}
    echo -e ">win\n$rna" > "$TMPDIR/tmp.fa"
    dg=$(RNAfold --noGU --noLP --noPS --paramFile="$PARAMFILE" < "$TMPDIR/tmp.fa" 2>/dev/null | \
         tail -n1 | sed -n 's/.*(\s*\([-0-9.]\+\)).*/\1/p')
    [[ -z "$dg" ]] && dg="NA"
    echo -e "$region\t$((i+1))\t$window\t$dg" >> "$DG_TSV"
  done
done < "$SEQ_TSV"

# Peak calling - choose local minima peak closest to position 1000 with dg ≤ -5.8
awk -v OFS="\t" -v dg_thresh=-5.8 -v OUT="$PEAKS_TSV" '
{
  reg=$1; pos=$2; dg=($4 == "NA" ? "NA" : $4+0);
  seqs[reg]=1;
  dgs[reg, NR]=$4;
  poss[reg, NR]=$2;
  ends[reg] = NR
}
END {
  for (r in seqs) {
    best_pos=-1; best_dg="NA"; min_dist=1e9; n=0; sum=0; sumsq=0;
    for (i=2; i<ends[r]; i++) {
      dgm=dgs[r,i]; dgp=dgs[r,i+1]; dgl=dgs[r,i-1];
      if (dgm != "NA") {
        sum += dgm; sumsq += dgm*dgm; n++;
        if (dgm <= dg_thresh && dgm < dgl && dgm < dgp) {
          dist = (poss[r,i] > 1000) ? poss[r,i] - 1000 : 1000 - poss[r,i];
          if (dist < min_dist) {
            min_dist = dist;
            best_pos = poss[r,i];
            best_dg = dgm;
          }
        }
      }
    }
    mean=(n>0?sum/n:"NA"); sd=(n>1?sqrt((sumsq-n*mean^2)/n):"NA");
    if (best_pos == -1) {
      peak_pos=1000; best_dg="NA"; status="unresolved";
    } else {
      peak_pos=best_pos; status="predicted";
    }
    split(r, parts, "::");
    print parts[1], peak_pos, status, best_dg > OUT;
  }
}' "$DG_TSV"

# GC content only around 100bp centered on peak
awk -v OFS="\t" '
{
  region=$1; seq=$2;
  center=1000; start=center-50; end=center+49;
  sub(/::.*/, "", region);
  gc=0; win=substr(seq, start+1, 100);
  for(i=1;i<=length(win);i++) {
    b=substr(win,i,1);
    if(b=="C"||b=="c"||b=="G"||b=="g") gc++
  }
  printf "%s\t%.6f\n", region, gc/100;
}' "$SEQ_TSV" > "$GC_TSV"

# Merge summary - output to raw summary first
echo -e "flank_name\tstart\tend\tpeak_call\tstatus\tdg_of_peak\tgc_content_around_peak" > "$SUMMARY_RAW"

awk '{gsub("\r",""); print}' "$BED_FLANKS" | \
  awk -v OFS="\t" '{sub(/::.*/, "", $4); print $4,$1,$2,$3}' | sort > "$TMPDIR/coords.tsv"

sort "$PEAKS_TSV" > "$TMPDIR/sorted_peaks.tsv"
sort "$GC_TSV" > "$TMPDIR/sorted_gc.tsv"
sort "$TMPDIR/coords.tsv" > "$TMPDIR/sorted_coords.tsv"

paste "$TMPDIR/sorted_peaks.tsv" "$TMPDIR/sorted_gc.tsv" "$TMPDIR/sorted_coords.tsv" | \
awk -F'\t' -v OFS='\t' '{
  print $1, $7, $8, $2, $3, $4, $6
}' >> "$SUMMARY_RAW"

# Clean and fix final output with absolute peak position
awk -F'\t' -v OFS='\t' '
BEGIN {
    print "flank_name", "chromosome", "peak_call", "status", "dg_of_peak", "gc_content_around_peak", "absolute_peak_position"
}
NR > 1 {
    flank = $1;
    chrom = $3;
    offset = $4 + 0;

    # Extract left and right coordinates from flank_name
    match(flank, /:([0-9]+)-([0-9]+)/, arr);
    left_coord = arr[1] + 0;
    right_coord = arr[2] + 0;

    # Determine flank type (3 or 5) from flank_name ending
    if (flank ~ /_flank_3$/) {
        abs_pos = right_coord - 1000 + offset;
    }
    else if (flank ~ /_flank_5$/) {
        abs_pos = left_coord - 1000 + offset;
    }
    else {
        abs_pos = left_coord + offset;
    }

    print flank, chrom, offset, $5, $6, $7, abs_pos;
}
' "$SUMMARY_RAW" > "$SUMMARY_OUT"

echo "[+] Done. Final cleaned summary: $SUMMARY_OUT"
echo "[+] Temp folder kept at: $TMPDIR"
