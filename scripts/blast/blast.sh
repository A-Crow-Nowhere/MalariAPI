#!/bin/bash

OUTDIR="blastOut"
FILTERED_FASTAS="filtered"
mkdir -p "$OUTDIR" "$FILTERED_FASTAS"

DB_PATH="$HOME/tools/blast/contamdb/out/contaminants_db"
SKIP_EXISTING=true

for fasta in ./fastas/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    outfile="$OUTDIR/${sample}.nt.blast.out"
    log="$OUTDIR/${sample}.log"
    filtered_fasta="$FILTERED_FASTAS/${sample}.filtered.fasta"

    if [[ "$SKIP_EXISTING" == true && -s "$outfile" ]]; then
        echo "$sample already completed, skipping."
        continue
    fi

    echo "$(date) - Filtering $sample" | tee "$log"

    # Filtering: Remove short or mostly-N reads
    echo "$(date) - Filtering $sample" | tee "$log"

    passed=0
    failed=0

    awk -v out="$filtered_fasta" '
    BEGIN { RS = ">" ; ORS = "" }
    NR > 1 {
        header = substr($0, 1, index($0, "\n") - 1)
        seq = substr($0, index($0, "\n") + 1)
        gsub(/\n/, "", seq)
        total = length(seq)
        ncount = gsub(/[Nn]/, "", seq)

        if (total >= 30 && ncount / total <= 0.5) {
            print ">" header "\n" seq "\n" >> out
            passed++
        } else {
            failed++
        }
    }
    END {
        printf("🔍 Filtering complete: %d passed, %d failed\n", passed, failed) > "/dev/stderr"
    }' "$fasta" 2>>"$log"


    blastn -query "$filtered_fasta" \
        -db "$DB_PATH" \
        -outfmt '6 qseqid sseqid pident length mismatch evalue sscinames' \
        -max_target_seqs 10 \
        -num_threads 4 \
        -dust no \
        -qcov_hsp_perc 80 \
        -evalue 1e-10 \
        -out "$outfile.tmp" >> "$log" 2>&1

    # Keep only top hit per read
    sort -k1,1 -k6,6g "$outfile.tmp" | awk '!seen[$1]++' > "$outfile"
    rm "$outfile.tmp"

    if [[ $? -ne 0 ]]; then
        echo "$(date) - ERROR for $sample" | tee -a "$log"
        continue
    fi

    echo "$(date) - Finished BLAST for $sample" | tee -a "$log"

    echo "🧮 Summarizing genus proportions for $sample:"
    python3 summarize_top_genus_proportions.py "$outfile" | tee -a "$log"
    echo ""

done

echo "✅ All samples processed or attempted."
