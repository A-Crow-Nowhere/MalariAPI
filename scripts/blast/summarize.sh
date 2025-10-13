#!/usr/bin/env python3

import os
import glob
import pandas as pd
import numpy as np
import argparse

def summarize_blast_folder(blast_dir):
    # Look for files ending with .genus_proportions.tsv
    pattern = os.path.join(blast_dir, "*.genus_proportions.tsv")
    all_files = sorted(glob.glob(pattern))
    if not all_files:
        raise FileNotFoundError(f"No .genus_proportions.tsv files found in directory: {blast_dir}")

    data = {}

    for file_path in all_files:
        sample_name = os.path.basename(file_path).replace(".genus_proportions.tsv", "")
        df = pd.read_csv(file_path, sep="\t")

        # Columns are lowercase: genus, percentage, number
        label_col = 'genus'
        pct_col = 'percentage'

        if label_col not in df.columns or pct_col not in df.columns:
            raise ValueError(f"Expected columns '{label_col}' and '{pct_col}' in {file_path}")

        for _, row in df.iterrows():
            category = row[label_col]
            percentage = row[pct_col]
            if category not in data:
                data[category] = []
            data[category].append(percentage)

    summary_rows = []
    for category, proportions in data.items():
        mean = np.mean(proportions)
        std = np.std(proportions)
        sample_count = len(proportions)
        details = ','.join([f"{p:.2f}" for p in proportions])
        summary_rows.append([category, mean, std, details, sample_count])

    summary_df = pd.DataFrame(summary_rows, columns=["Category", "Mean%", "StdDev%", "Total", "N_samples"])
    summary_df.sort_values(by="Mean%", ascending=False, inplace=True)
    summary_df.to_csv("summary_total.tsv", sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize genus proportions from BLAST output.")
    parser.add_argument(
        "--blast-dir",
        type=str,
        default="blastOut",
        help="Directory containing .genus_proportions.tsv files (default: blastOut)"
    )
    args = parser.parse_args()
    summarize_blast_folder(args.blast_dir)
