for f in *.tsv; do
  sample="${f%%.*}"
  awk -v S="$sample" 'BEGIN{OFS="\t"}
  NR==1 { print $0, "sample_name"; next }
  { print $0, S }
  ' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done
