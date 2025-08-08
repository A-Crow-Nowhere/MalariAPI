#!/usr/bin/env bash
set -euo pipefail

# Default to 'blastOut' if no argument is provided
BLAST_DIR=${1:-blastOut}

python3 - "$BLAST_DIR" <<EOF
import os
import sys
import pandas as pd

blast_dir = sys.argv[1]
rows = []
all_samples = []

for filename in os.listdir(blast_dir):
    if filename.endswith(".genus_proportions.tsv"):
        sample_name = filename.replace(".genus_proportions.tsv", "")
        all_samples.append(sample_name)
        filepath = os.path.join(blast_dir, filename)

        try:
            df = pd.read_csv(filepath, sep="\\t")
        except Exception:
            df = pd.DataFrame()

        if df.empty or not {"genus", "percentage", "number"}.issubset(df.columns):
            # No hits or malformed — still include empty entry
            rows.append({
                "Sample": sample_name,
                "Category": "NoHits",
                "Percentage": 0,
                "Number": 0,
                "TotalReads": 0
            })
            continue

        total = df["number"].sum()

        for _, row in df.iterrows():
            category = row["genus"]
            percent = row["percentage"]
            count = row["number"]
            rows.append({
                "Sample": sample_name,
                "Category": category,
                "Percentage": percent,
                "Number": count,
                "TotalReads": total
            })

# Convert to dataframe
summary_df = pd.DataFrame(rows)

# Ensure all samples are included in pivot table, even if completely missing
if not summary_df.empty:
    summary_df.to_csv("summary_by_sample.tsv", sep="\\t", index=False)
    print("✅ Summary by sample saved to summary_by_sample.tsv")

    pivot_df = summary_df.pivot_table(index="Sample", columns="Category", values="Percentage", fill_value=0)
    pivot_df["TotalReads"] = summary_df.groupby("Sample")["TotalReads"].first().reindex(pivot_df.index, fill_value=0)
else:
    # Empty result fallback
    pivot_df = pd.DataFrame({"TotalReads": [0 for _ in all_samples]}, index=all_samples)

pivot_df.to_csv("summary_matrix.tsv", sep="\\t")
print("✅ Summary matrix saved to 'summary_matrix.tsv'")

# Summary stats (only for non-empty categories)
if not summary_df.empty:
    stats_df = summary_df.groupby("Category").agg({
        "Percentage": ["mean", "std", lambda x: ",".join(f"{v:.2f}" for v in x)],
        "Sample": "count"
    })
    stats_df.columns = ["Mean", "StdDev", "AllPercentages", "SampleCount"]
    stats_df.to_csv("summary_total.tsv", sep="\\t")
    print("✅ Summary statistics saved to 'summary_total.tsv'")
else:
    with open("summary_total.tsv", "w") as f:
        f.write("No data available\\n")
    print("⚠️ No hits found in any sample. Empty summary written.")
EOF
