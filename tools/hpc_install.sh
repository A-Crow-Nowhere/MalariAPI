1. TO MAKE SVCROWS FILE FROM VCF FILES IN A DIRECTORY:

#!/bin/bash

for vcf in *.vcf; do
    # remove the .vcf extension
    base="${vcf%.vcf}"
    
    echo "Processing $vcf ..."
    vcf_to_crows.sh -V "$vcf" -O "$base" -s "$base"
done

2. FILTERS, DELS, DUPS, INS, AND INV INTO SEPARATE FOLDERS, COUNTS EACH VARIANT TYPE, AND SAVES IT TO A FILE

#!/bin/bash

# Name of summary output file
summary_file="summary_counts.tsv"

# Overwrite existing summary file and add header
echo -e "Sample\tDEL\tDEL_Reads\tDUP\tDUP_Reads\tINS\tINS_Reads\tINV\tINV_Reads" > "$summary_file"

# Create output folders for each SV type
mkdir -p dels_only dups_only ins_only inv_only

# Loop through each .crows.tsv file in the current directory
for file in *.crows.tsv; do
    base="${file%.crows.tsv}"
    echo "Processing $file ..."

    # Step 1 — Create a TEMP fixed version of the file with NumReads replaced using SUPPORT
    fixed="${base}_fixed.tmp"

    awk 'BEGIN{FS=OFS="\t"}
        NR==1 { print; next }
        {
            match($7, /SUPPORT=([0-9]+)/, m);
            if (m[1] != "") {
                $12 = m[1];
            }
            print;
        }' "$file" > "$fixed"

    # Step 2 — Filter SV types *from the fixed file*
    awk 'NR==1 || $5 == "DEL"' "$fixed" > "dels_only/${base}_DEL.tsv"
    awk 'NR==1 || $5 == "DUP"' "$fixed" > "dups_only/${base}_DUP.tsv"
    awk 'NR==1 || $5 == "INS"' "$fixed" > "ins_only/${base}_INS.tsv"
    awk 'NR==1 || $5 == "INV"' "$fixed" > "inv_only/${base}_INV.tsv"

    # Step 3 — Count events for each SV type
    dels=$(awk '$5 == "DEL"' "$fixed" | wc -l)
    dups=$(awk '$5 == "DUP"' "$fixed" | wc -l)
    ins=$(awk '$5 == "INS"' "$fixed" | wc -l)
    inv=$(awk '$5 == "INV"' "$fixed" | wc -l)

    # Step 4 — Sum NumReads for each SV type
    del_reads=$(awk '$5 == "DEL" {sum += $12} END{print sum+0}' "$fixed")
    dup_reads=$(awk '$5 == "DUP" {sum += $12} END{print sum+0}' "$fixed")
    ins_reads=$(awk '$5 == "INS" {sum += $12} END{print sum+0}' "$fixed")
    inv_reads=$(awk '$5 == "INV" {sum += $12} END{print sum+0}' "$fixed")

    # Print per-file summary
    printf "  DELs: %d (Reads: %d)\n" "$dels" "$del_reads"
    printf "  DUPs: %d (Reads: %d)\n" "$dups" "$dup_reads"
    printf "  INSs: %d (Reads: %d)\n" "$ins" "$ins_reads"
    printf "  INVs: %d (Reads: %d)\n" "$inv" "$inv_reads"
    echo "-------------------------------------"

    # Append summary row to TSV
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
        "$base" "$dels" "$del_reads" "$dups" "$dup_reads" "$ins" "$ins_reads" "$inv" "$inv_reads" >> "$summary_file"

    # Clean up temp file
    rm "$fixed"
done

echo "✅ All files processed successfully!"
echo "📊 Summary written to: $summary_file"


4. VAF

#!/bin/bash

make_vaf_summary() {
    folder=$1
    outfile="${folder%_only}_VAF.tsv"
    echo "Creating $outfile ..."

    tmpdir=$(mktemp -d)

    # Loop through all TSVs in the folder
    for file in ${folder}/*.tsv; do
        # Use the FULL filename, minus `.crows.tsv`
        sample=$(basename "$file" .tsv)
        sample=${sample%.crows}
        sample=${sample%.crows.tsv}

        # Extract VAF from column 7
        awk -F'\t' '
            NR>1 {
                v = "NA"
                if (match($7, /VAF=([0-9]*\.?[0-9]+)/, m)) {
                    v = m[1]
                }
                print v
            }
        ' "$file" > "$tmpdir/$sample.txt"
    done

    # Build header
    header=$(for f in "$tmpdir"/*.txt; do basename "$f" .txt; done | paste -sd'\t')
    echo -e "$header" > "$outfile"

    # Paste columns row-by-row
    paste "$tmpdir"/*.txt >> "$outfile"

    rm -r "$tmpdir"
    echo "Wrote $outfile"
}

make_vaf_summary dels_only
make_vaf_summary dups_only
make_vaf_summary ins_only
make_vaf_summary inv_only

echo "All VAF summaries created."



3. MAKES LENGTH SUMMARIES FOR ALL DELS DUPS IVS INS 

#!/bin/bash

# Function to build summary table for a given variant type folder
make_summary() {
    folder=$1
    outfile="${folder%_only}_lengths.tsv"
    echo "Creating $outfile ..."

    # Temporary directory for intermediate files
    tmpdir=$(mktemp -d)

    # Extract sample name and length column from each TSV
    for file in ${folder}/*.tsv; do
        sample=$(basename "$file" | sed 's/_.*//')   # sample name before first underscore
        awk 'NR>1 {print $4}' "$file" > "$tmpdir/$sample.txt"
    done

    # Paste columns together and add headers
    paste -d'\t' <(for f in "$tmpdir"/*.txt; do basename "$f" .txt; done) \
        <(printf "") > /dev/null 2>&1  # dummy call to avoid syntax highlighting issues

    # Build header row
    header=$(for f in "$tmpdir"/*.txt; do basename "$f" .txt; done | paste -sd'\t')
    echo -e "$header" > "$outfile"

    # Merge columns by row position
    paste "$tmpdir"/*.txt >> "$outfile"

    # Clean up
    rm -r "$tmpdir"
    echo "✅ Wrote $outfile"
}

# Generate summaries for each variant type
make_summary dels_only
make_summary dups_only
make_summary ins_only
make_summary inv_only

echo "🎉 All summary tables created!"
