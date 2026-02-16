#!/usr/bin/env python3

import os
import glob
import argparse
import pandas as pd
import numpy as np

def summarize_folder(blast_dir, out_path):
    files = sorted(glob.glob(os.path.join(blast_dir, "*.genus_proportions.tsv")))
    if not files:
        raise FileNotFoundError(f"No genus_proportions.tsv files found in {blast_dir}")

    data = {}

    for fp in files:
        df = pd.read_csv(fp, sep="\t")
        if "genus" not in df.columns or "percentage" not in df.columns:
            raise ValueError(f"Malformed file: {fp}")

        for _, row in df.iterrows():
            genus = row["genus"]
            pct = row["percentage"]
            data.setdefault(genus, []).append(pct)

    rows = []
    for genus, values in data.items():
        rows.append([
            genus,
            float(np.mean(values)),
            float(np.std(values)),
            ",".join(f"{v:.2f}" for v in values),
            len(values),
        ])

    out_df = pd.DataFrame(
        rows,
        columns=["Category", "Mean%", "StdDev%", "Total", "N_samples"]
    ).sort_values("Mean%", ascending=False)

    out_df.to_csv(out_path, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast-dir", required=True)
    ap.add_argument("--out", default="summary_total.tsv")
    args = ap.parse_args()

    summarize_folder(args.blast_dir, args.out)

if __name__ == "__main__":
    main()
