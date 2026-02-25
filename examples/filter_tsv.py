#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_tsv", required=True, help="Input TSV file")
    ap.add_argument("--out", dest="out_tsv", required=True, help="Output TSV file")
    ap.add_argument("--min", dest="min_val", type=float, default=0.0, help="Minimum value")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t")
    df = df[df.select_dtypes(include="number").ge(args.min_val).all(axis=1)]
    df.to_csv(args.out_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
