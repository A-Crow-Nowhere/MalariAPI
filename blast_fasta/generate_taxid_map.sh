#!/bin/bash

FASTA="contaminants_all.fa"
KEYWORDS="taxid_keywords.tsv"
OUT="taxid_map.txt"

# Make sure files exist
if [[ ! -f $FASTA || ! -f $KEYWORDS ]]; then
  echo "Missing $FASTA or $KEYWORDS"
  exit 1
fi

# Build awk-friendly lookup string
awk -v kfile="$KEYWORDS" '
  BEGIN {
    while ((getline < kfile) > 0) {
      kwd[tolower($1)] = $2
    }
  }
  /^>/ {
    header = substr($0, 2)
    split(header, a, " ")
    id = a[1]
    tax = 0
    for (k in kwd) {
      if (tolower($0) ~ k) {
        tax = kwd[k]
        break
      }
    }
    if (tax > 0)
      print id "\t" tax
    else
      print id "\t0"  # fallback if no match
  }
' "$FASTA" > "$OUT"

echo "✅ Generated $OUT"
