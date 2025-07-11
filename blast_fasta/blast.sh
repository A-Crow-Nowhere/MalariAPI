#!/bin/bash

OUTDIR="blast_nt_out"
mkdir -p "$OUTDIR"

DB_PATH="$HOME/contamdb/contaminants_db"
SKIP_EXISTING=true

for fasta in *.fasta; do
    sample=$(basename "$fasta" .fasta)
    outfile="$OUTDIR/${sample}.nt.blast.out"
    log="$OUTDIR/${sample}.log"

    if [[ "$SKIP_EXISTING" == true && -s "$outfile" ]]; then
        echo "$sample already completed, skipping."
        continue
    fi

    echo "$(date) - BLASTing $sample" | tee "$log"

    # Run blastn with multiple targets to allow best match sorting
    blastn -query "$fasta" \
        -db "$DB_PATH" \
        -outfmt '6 qseqid sseqid pident length mismatch evalue sscinames' \
        -max_target_seqs 10 \
        -num_threads 4 \
        -dust no \
        -qcov_hsp_perc 80 \
        -evalue 1e-10 \
        -out "$outfile.tmp" >> "$log" 2>&1

    # Keep only the top hit per read (lowest evalue)
    sort -k1,1 -k6,6g "$outfile.tmp" | awk '!seen[$1]++' > "$outfile"
    rm "$outfile.tmp"

    if [[ $? -ne 0 ]]; then
        echo "$(date) - ERROR for $sample" | tee -a "$log"
        continue
    fi

    echo "$(date) - Finished BLAST for $sample" | tee -a "$log"

    # Run python genus summarization on the filtered output and capture
    echo "🧮 Summarizing genus proportions for $sample:"
    python3 summarize_top_genus_proportions.py "$outfile" | tee -a "$log"
    echo ""

done

echo "✅ All samples processed or attempted."
