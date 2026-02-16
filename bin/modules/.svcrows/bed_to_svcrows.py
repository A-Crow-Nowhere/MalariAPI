#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
from typing import List, Optional, Tuple

SVCROWS_HEADER = [
    "Chr","Start","End","Length","Type","ID","Var1","Var2","Var3","IsKnown","QScore","NumReads"
]

def open_maybe_gz(path: str, mode: str):
    # mode: "rt" or "wt"
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, newline="")

def norm_type(raw: str) -> Optional[str]:
    if raw is None:
        return None
    s = raw.strip().lower()
    s = re.sub(r"[^a-z]", "", s)  # drop punctuation/underscores/etc
    if s in ("dup", "duplication", "dupl", "tandemdup", "tandemduplication"):
        return "DUP"
    if s in ("del", "deletion", "loss"):
        return "DEL"
    if s in ("inv", "inversion"):
        return "INV"
    # allow already-normalized
    if s in ("dup", "del", "inv"):
        return s.upper()
    if s in ("dup", "del", "inv"):
        return s.upper()
    if s in ("dup",):
        return "DUP"
    if s in ("del",):
        return "DEL"
    if s in ("inv",):
        return "INV"
    # also accept literal "DUP"/"DEL"/"INV" if preserved
    if raw.strip().upper() in ("DUP","DEL","INV"):
        return raw.strip().upper()
    return None

def parse_int(x: str) -> Optional[int]:
    try:
        # accept floats that are actually ints
        if x is None:
            return None
        x = x.strip()
        if x == "":
            return None
        if re.match(r"^-?\d+(\.0+)?$", x):
            return int(float(x))
        # if it's clearly an int
        return int(x)
    except Exception:
        return None

def parse_float(x: str) -> Optional[float]:
    try:
        if x is None:
            return None
        x = x.strip()
        if x == "":
            return None
        return float(x)
    except Exception:
        return None

def detect_header(first_row: List[str]) -> bool:
    # Heuristic: if col1 isn't a chrom-like token OR col2/col3 aren't ints
    if not first_row:
        return False
    c1 = first_row[0].strip()
    c2 = first_row[1].strip() if len(first_row) > 1 else ""
    c3 = first_row[2].strip() if len(first_row) > 2 else ""
    if parse_int(c2) is None or parse_int(c3) is None:
        return True
    # if it literally says "chr" or "chrom"
    if c1.lower() in ("chr","chrom","chromosome"):
        return True
    return False

def sanitize_kv_key(k: str) -> str:
    k = k.strip()
    k = re.sub(r"\s+", "_", k)
    k = re.sub(r"[^A-Za-z0-9_]", "", k)
    if k == "":
        k = "COL"
    return k

def build_var1(
    svtype: str,
    svlen: int,
    end: int,
    source: str,
    orig_header: Optional[List[str]],
    row: List[str],
) -> str:
    # VCF-like INFO base
    parts = [f"SVTYPE={svtype}", f"SVLEN={svlen}", f"END={end}", f"SOURCE={source}"]

    # Append original fields as FIELD=VALUE
    if orig_header and len(orig_header) == len(row):
        for k, v in zip(orig_header, row):
            kk = sanitize_kv_key(k)
            vv = (v or "").strip()
            if vv == "":
                continue
            # avoid duplicating base keys too noisily
            if kk.upper() in ("CHR","CHROM","START","BEGIN","END","SVTYPE","SVLEN","TYPE","SOURCE"):
                continue
            # keep values but strip semicolons
            vv = vv.replace(";", ",")
            parts.append(f"{kk}={vv}")
    else:
        # no header: keep positional cols
        for i, v in enumerate(row, start=1):
            vv = (v or "").strip()
            if vv == "":
                continue
            vv = vv.replace(";", ",")
            parts.append(f"COL{i}={vv}")

    return ";".join(parts)

def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--sample", default="Sample")
    ap.add_argument("--source", default="bed")
    ap.add_argument("--type-col", type=int, default=4)  # 1-based
    ap.add_argument("--id-col", type=int, default=0)     # 1-based; 0=none
    ap.add_argument("--su-col", type=int, default=0)     # 1-based; 0=none
    ap.add_argument("--keep-types", default="DEL,DUP,INV")
    ap.add_argument("--header", default="auto", choices=["auto","yes","no"])
    ap.add_argument("--out-id-prefix", default="")
    ap.add_argument("--require-bounds", action="store_true")
    ap.add_argument("--drop-nonstandard", action="store_true")
    ap.add_argument("-h","--help", action="help")
    args = ap.parse_args()

    keep = set([t.strip().upper() for t in args.keep_types.split(",") if t.strip()])

    type_idx = args.type_col - 1
    id_idx = args.id_col - 1 if args.id_col > 0 else None
    su_idx = args.su_col - 1 if args.su_col > 0 else None

    # Read first row to figure header
    with open_maybe_gz(args.inp, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        first = None
        for row in reader:
            if not row or all((c.strip()=="" for c in row)):
                continue
            first = row
            break

        if first is None:
            raise SystemExit("Empty input")

        has_header = False
        if args.header == "yes":
            has_header = True
        elif args.header == "no":
            has_header = False
        else:
            has_header = detect_header(first)

    # Now do full pass
    with open_maybe_gz(args.inp, "rt") as f_in, open_maybe_gz(args.out, "wt") as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")

        writer.writerow(SVCROWS_HEADER)

        orig_header = None
        wrote = 0
        kept = 0

        for row in reader:
            if not row or all((c.strip()=="" for c in row)):
                continue

            if orig_header is None and has_header:
                orig_header = row
                continue

            # Need at least 3 columns
            if len(row) < 3:
                continue

            chrom = row[0].strip()
            start = parse_int(row[1])
            end = parse_int(row[2])

            if args.require_bounds:
                if start is None or end is None:
                    continue
                if start >= end:
                    continue

            # Type
            raw_type = row[type_idx] if 0 <= type_idx < len(row) else ""
            svtype = norm_type(raw_type)
            if svtype is None:
                if args.drop_nonstandard:
                    continue
                else:
                    # keep but mark unknown
                    continue

            if svtype not in keep:
                continue

            length = (end - start) if (start is not None and end is not None) else 0
            if svtype == "DEL":
                svlen = -abs(length)
            else:
                svlen = abs(length)

            # Support / reads
            su_val = 0
            if su_idx is not None and 0 <= su_idx < len(row):
                # accept int or float, cast down
                su_parsed = parse_float(row[su_idx])
                if su_parsed is not None:
                    su_val = int(round(su_parsed))
            num_reads = su_val

            # ID
            if id_idx is not None and 0 <= id_idx < len(row):
                rid = row[id_idx].strip()
                if rid == "":
                    rid = str(kept + 1)
            else:
                rid = str(kept + 1)

            if args.out_id_prefix:
                rid_out = f"{args.out_id_prefix}{rid}"
            else:
                # SVCROWS examples show numeric-ish IDs; keep clean
                rid_out = rid

            var1 = build_var1(
                svtype=svtype,
                svlen=svlen,
                end=end if end is not None else 0,
                source=args.source,
                orig_header=orig_header,
                row=row
            )
            var2 = f"GT=./.;SU={num_reads}"

            out_row = [
                chrom,
                str(start if start is not None else ""),
                str(end if end is not None else ""),
                str(svlen),
                svtype,
                rid_out,
                var1,
                var2,
                args.sample,
                "FALSE",
                "0",
                str(num_reads),
            ]
            writer.writerow(out_row)
            wrote += 1
            kept += 1

    # done

if __name__ == "__main__":
    main()
