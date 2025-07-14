#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="./blastOut"
OUT_SUMMARY_BY_SAMPLE="summary_by_sample.tsv"
OUT_SUMMARY_TOTAL="summary_total.tsv"

# Clear output files and add headers
echo -e "Sample\tCategory\tPercentage" > "$OUT_SUMMARY_BY_SAMPLE"
echo -e "Category\tMeanPercentage\tMedianPercentage\tNumberOfSamples" > "$OUT_SUMMARY_TOTAL"

# Temporary file to collect all category percentages with sample name
TMP_ALL="all_samples_categories.tsv"
rm -f "$TMP_ALL"

# Loop over all TSV files in input dir
for f in "$INPUT_DIR"/*.tsv; do
    # Extract sample name (filename without directory and extension)
    sample=$(basename "$f" .tsv)
    
    # Read each line: Category \t Percentage (assuming two columns per file)
    # Append sample name as first column
    awk -v sample="$sample" -F'\t' 'NF>=2 {print sample"\t"$1"\t"$2}' "$f" >> "$TMP_ALL"
done

# Append all entries to summary_by_sample.tsv
cat "$TMP_ALL" >> "$OUT_SUMMARY_BY_SAMPLE"

# Now compute summary_total.tsv

# awk script to calculate mean, median, and count for each category
awk -F'\t' '
{
    cat=$2; perc=$3+0;
    vals[cat][counts[cat]++] = perc;
    sums[cat] += perc;
}
END {
    for (cat in sums) {
        n = counts[cat];
        # Compute mean
        mean = sums[cat] / n;

        # Compute median
        # sort vals[cat]
        for(i=0; i<n; i++) arr[i] = vals[cat][i];
        asort(arr);
        if (n % 2 == 1) {
            median = arr[int(n/2)];
        } else {
            median = (arr[n/2 -1] + arr[n/2]) / 2;
        }

        printf "%s\t%.3f\t%.3f\t%d\n", cat, mean, median, n;
    }
}
' "$TMP_ALL" | sort > "$OUT_SUMMARY_TOTAL"

# Clean up temp file
rm -f "$TMP_ALL"

echo "✅ Summaries written:"
echo " - $OUT_SUMMARY_BY_SAMPLE"
echo " - $OUT_SUMMARY_TOTAL"
