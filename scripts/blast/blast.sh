#!/bin/bash

OUTDIR="blastOut"
mkdir -p "$OUTDIR"

DB_PATH="$HOME/tools/blast/contamdb/out/contaminants_db"
SKIP_EXISTING=true

MIN_MEM_KB=700000  # 700 MB threshold
SLEEP_AFTER_EXIT=30 # seconds to wait before restart on low memory

DEBUG=0
if [[ "$1" == "--debug" ]]; then
  DEBUG=1
fi

process_samples() {
  for fasta in ./fastas/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    outfile="$OUTDIR/${sample}.nt.blast.out"
    log="$OUTDIR/${sample}.log"
    filtered_fasta="$OUTDIR/${sample}.filtered.fasta"

    if [[ "$SKIP_EXISTING" == true && -s "$outfile" ]]; then
        echo "$sample already completed, skipping."
        continue
    fi

    echo "$(date) - Filtering $sample" | tee "$log"

    # Filter: remove reads <100 bp or >50% Ns
    awk -v out="$filtered_fasta" '
    BEGIN { RS = ">" ; ORS = "" }
    NR > 1 {
        header = substr($0, 1, index($0, "\n") - 1)
        seq = substr($0, index($0, "\n") + 1)
        gsub(/\n/, "", seq)
        total = length(seq)
        ncount = gsub(/[Nn]/, "", seq)

        if (total >= 100 && ncount / total <= 0.5) {
            print ">" header "\n" seq "\n" >> out
            passed++
        } else {
            failed++
        }
    }
    END {
        printf("✅ Filtering complete: %d passed, %d failed\n", passed, failed) > "/dev/stderr"
    }' "$fasta" 2>>"$log"

    if [[ ! -s "$filtered_fasta" ]]; then
        echo "⚠️ No reads passed filter for $sample. Skipping." | tee -a "$log"
        continue
    fi

    # Determine number of sequences to blast
    total_seqs=$(grep -c '^>' "$filtered_fasta")
    if [[ $DEBUG -eq 1 ]]; then
        num_seqs_to_blast=100
    else
        num_seqs_to_blast=5000
    fi

    blast_input="$filtered_fasta"
    if [[ $num_seqs_to_blast -lt $total_seqs ]]; then
        blast_input="$OUTDIR/${sample}.blast_input.fasta"
        awk -v n=$num_seqs_to_blast 'BEGIN{RS=">"; ORS=""} NR==1{next} NR<=n+1{print ">"$0}' "$filtered_fasta" > "$blast_input"
        echo "DEBUG mode: blasting only first $num_seqs_to_blast sequences for $sample" | tee -a "$log"
    fi

    echo "$(date) - BLASTing $sample" | tee -a "$log"
    echo "?? Running blastn on $sample with up to $num_seqs_to_blast sequences" | tee -a "$log"

    free_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    if [[ $free_kb -lt $MIN_MEM_KB ]]; then
        echo "⚠️ Low available memory (${free_kb} KB) before BLAST for $sample. Exiting to avoid crash." | tee -a "$log"
        rm -f "$blast_input"
        return 1  # signal to restart
    fi

    blastn -query "$blast_input" \
        -db "$DB_PATH" \
        -outfmt '6 qseqid sseqid pident length mismatch evalue sscinames' \
        -max_target_seqs 10 \
        -num_threads 8 \
        -dust no \
        -qcov_hsp_perc 90 \
        -evalue 1e-20 \
        -out "$outfile.tmp" >> "$log" 2>&1

    blast_exit=$?

    if [[ $blast_exit -ne 0 ]]; then
        echo "$(date) - ❌ BLAST failed with exit $blast_exit for $sample. Skipping." | tee -a "$log"
        rm -f "$outfile.tmp" "$blast_input"
        continue
    fi

    sort -k1,1 -k6,6g "$outfile.tmp" | awk '!seen[$1]++' > "$outfile"
    rm -f "$outfile.tmp" "$blast_input"

    echo "$(date) - ✅ Finished BLAST for $sample" | tee -a "$log"

    echo "?? Summarizing genus proportions for $sample:" | tee -a "$log"
    python3 summarize_top_genus_proportions.py "$outfile" | tee -a "$log"
    echo "" | tee -a "$log"

    echo "$(date) - Memory usage after $sample:" | tee -a "$log"
    free -h | tee -a "$log"

    free_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
    if [[ $free_kb -lt $MIN_MEM_KB ]]; then
        echo "⚠️ Low available memory (${free_kb} KB) after processing $sample. Exiting to avoid crash." | tee -a "$log"
        return 1  # signal to restart
    fi

    echo "Sleeping 2 seconds to free resources..." | tee -a "$log"
    sleep 2
  done

  return 0
}

echo "Starting BLAST loop. Debug mode: $DEBUG"
while true; do
  process_samples
  exitcode=$?
  if [[ $exitcode -eq 0 ]]; then
    echo "✅ All samples processed."
    break
  else
    echo "⚠️ Exiting early due to low memory. Sleeping $SLEEP_AFTER_EXIT seconds before retry..."
    sleep $SLEEP_AFTER_EXIT
    echo "Restarting BLAST loop..."
  fi
done
