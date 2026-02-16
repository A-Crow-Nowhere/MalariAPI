#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import tempfile
from typing import Optional, List

import numpy as np
import pandas as pd
import pysam


# ============================================================
# IO helpers
# ============================================================

def _is_gzip_by_magic(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def read_tsv_any(path: str, header="infer") -> pd.DataFrame:
    """
    Read TSV. If file ends with .gz but isn't actually gzipped, detect by magic bytes.
    """
    is_gz = _is_gzip_by_magic(path)
    if is_gz:
        return pd.read_csv(path, sep="\t", header=header, compression="gzip")
    return pd.read_csv(path, sep="\t", header=header)


def write_tsv_any(df: pd.DataFrame, path: str) -> None:
    """
    Write TSV; if path ends with .gz, gzip it. Always write NA explicitly for missing.
    """
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            df.to_csv(f, sep="\t", index=False, na_rep="NA")
    else:
        df.to_csv(path, sep="\t", index=False, na_rep="NA")


def write_bedgraph_any(df_bg: pd.DataFrame, path: str) -> None:
    """
    Write bedGraph (no header). Always write NA explicitly for missing.
    """
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            df_bg.to_csv(f, sep="\t", header=False, index=False, na_rep="NA")
    else:
        df_bg.to_csv(path, sep="\t", header=False, index=False, na_rep="NA")


def _atomic_write_no_follow(write_fn, df: pd.DataFrame, out_path: str) -> None:
    """
    Write output without following an existing symlink at out_path.
    Writes to temp in same dir, then os.replace (replaces symlink itself if present).
    """
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)

    fd, tmp_path = tempfile.mkstemp(prefix=".tmp.", dir=out_dir)
    os.close(fd)

    try:
        write_fn(df, tmp_path)
        os.replace(tmp_path, out_path)
    finally:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass


# ============================================================
# Segment parsing
# ============================================================

def _maybe_infer_ratio_col(df: pd.DataFrame, y_col: int) -> Optional[int]:
    """
    If there is a column immediately after y_col that looks like ratio (positive, ~2^y),
    use it; else return None (we'll compute ratio).
    """
    if y_col + 1 >= df.shape[1]:
        return None

    y = pd.to_numeric(df.iloc[:, y_col], errors="coerce").to_numpy(dtype=float)
    r = pd.to_numeric(df.iloc[:, y_col + 1], errors="coerce").to_numpy(dtype=float)

    finite = np.isfinite(y) & np.isfinite(r)
    if finite.sum() < 50:
        return None

    yy = y[finite]
    rr = r[finite]

    # ratio should be >0 and roughly exp2(y)
    if np.mean(rr > 0) < 0.95:
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.nanmedian(np.abs(np.log2(rr) - yy))
    if not np.isfinite(err):
        return None

    # If the median absolute error is small, trust it
    if err <= 0.25:
        return y_col + 1

    return None


def _infer_y_col(df: pd.DataFrame, candidate_cols: List[int], n_check: int = 200) -> int:
    """
    Infer which column is y (log2 ratio) from a headerless segments BED-like table.

    IMPORTANT: We avoid BED score-like columns (often constant 0) by preferring columns
    that have both negative and positive values (y is log2 ratio and should cross 0).
    """
    n = min(n_check, len(df))
    if n == 0:
        raise RuntimeError("Empty segments file; cannot infer y column.")

    best_col = None
    best_score = -1e18

    for c in candidate_cols:
        vals = pd.to_numeric(df.iloc[:n, c], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(vals)
        if finite.sum() < max(10, int(0.5 * n)):
            continue

        v = vals[finite]

        has_neg = np.any(v < 0)
        has_pos = np.any(v > 0)
        sign_bonus = 2.0 if (has_neg and has_pos) else (-2.0)

        uniq = len(np.unique(np.round(v, 6)))
        uniq_bonus = min(1.0, uniq / 20.0)

        frac_in = np.mean((v >= -20.0) & (v <= 20.0))
        med_abs = np.median(np.abs(v))

        score = (2.0 * frac_in) + (1.0 / (1.0 + med_abs)) + sign_bonus + uniq_bonus

        if score > best_score:
            best_score = score
            best_col = c

    if best_col is None:
        raise RuntimeError(
            "Could not infer y (log2 ratio) column. "
            "Your segments.bed column layout is unusual—please inspect columns."
        )
    return int(best_col)


def read_segments_bed(path: str) -> pd.DataFrame:
    """
    Match the ORIGINAL contract first:

    segments.bed(.gz) expected columns (0-based):
      0 chrom
      1 start
      2 end
      3 name
      4 bed_score
      5 strand
      6 y (log2 ratio)
      7 ratio (optional; computed if missing)
      9 n_probes (optional; if present)

    If the file doesn't have >=7 columns, fall back to inference.
    """
    df = read_tsv_any(path, header=None)

    if df.shape[1] < 4:
        raise RuntimeError("segments.bed must have >= 4 columns (chrom,start,end,name).")

    out = pd.DataFrame({
        "chrom": df.iloc[:, 0].astype(str),
        "start": pd.to_numeric(df.iloc[:, 1], errors="coerce"),
        "end": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
        "name": df.iloc[:, 3].astype(str),
    })

    if out["start"].isna().any() or out["end"].isna().any():
        raise RuntimeError("segments.bed has non-numeric start/end in columns 1/2—unexpected format.")

    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)

    # ---- Prefer original layout: y at col 6 if present ----
    if df.shape[1] >= 7:
        y_col = 6
        out["y"] = pd.to_numeric(df.iloc[:, y_col], errors="coerce").astype(float)

        if df.shape[1] >= 8:
            out["ratio"] = pd.to_numeric(df.iloc[:, 7], errors="coerce").astype(float)
        else:
            out["ratio"] = np.power(2.0, out["y"].to_numpy(dtype=float))

        if df.shape[1] >= 10:
            out["n_probes"] = pd.to_numeric(df.iloc[:, 9], errors="coerce").astype("Int64")
        else:
            out["n_probes"] = pd.Series([pd.NA] * len(out), dtype="Int64")

        return out

    # ---- Fallback: inference ----
    candidate_cols = list(range(4, df.shape[1]))
    y_col = _infer_y_col(df, candidate_cols)
    out["y"] = pd.to_numeric(df.iloc[:, y_col], errors="coerce").astype(float)

    ratio_col = _maybe_infer_ratio_col(df, y_col)
    if ratio_col is not None:
        out["ratio"] = pd.to_numeric(df.iloc[:, ratio_col], errors="coerce").astype(float)
    else:
        out["ratio"] = np.power(2.0, out["y"].to_numpy(dtype=float))

    out["n_probes"] = pd.Series([pd.NA] * len(out), dtype="Int64")
    return out


# ============================================================
# Boundary support TSV
# ============================================================

def read_boundary_support(path: str) -> pd.DataFrame:
    """
    boundary_support.tsv expected at least:
      chrom, name, combined_z_chr
    """
    df = read_tsv_any(path, header="infer")
    required = {"chrom", "name", "combined_z_chr"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"boundary_support missing required columns: {sorted(missing)}")
    df["chrom"] = df["chrom"].astype(str)
    df["name"] = df["name"].astype(str)
    df["combined_z_chr"] = pd.to_numeric(df["combined_z_chr"], errors="coerce").astype(float)
    return df


# ============================================================
# BAM counting helpers
# ============================================================

def _open_bam(path: str) -> Optional[pysam.AlignmentFile]:
    if not path:
        return None
    if not os.path.exists(path):
        raise RuntimeError(f"BAM not found: {path}")
    return pysam.AlignmentFile(path, "rb")


def _bam_has_contig(bam: pysam.AlignmentFile, chrom: str) -> bool:
    return chrom in bam.references


def count_reads_in_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    read1_only: bool = True,
    drop_dup: bool = False,
    include_supp: bool = False,
    include_secondary: bool = False,
) -> int:
    if start < 0:
        start = 0
    if end < start:
        end = start
    if not _bam_has_contig(bam, chrom):
        return 0

    n = 0
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped:
            continue
        if (not include_secondary) and r.is_secondary:
            continue
        if (not include_supp) and r.is_supplementary:
            continue
        if drop_dup and r.is_duplicate:
            continue

        if read1_only:
            if r.is_paired and (not r.is_read1):
                continue
        n += 1
    return n


# ============================================================
# Core logic helpers
# ============================================================

def min2(a: Optional[float], b: Optional[float]) -> float:
    vals = []
    if a is not None and not (isinstance(a, float) and math.isnan(a)):
        vals.append(float(a))
    if b is not None and not (isinstance(b, float) and math.isnan(b)):
        vals.append(float(b))
    if not vals:
        return float("nan")
    return float(min(vals))


def z_tier(z: float, weak_z: float, strong_z: float) -> str:
    if not math.isfinite(z):
        return "NA"
    if z >= strong_z:
        return "strong"
    if z >= weak_z:
        return "weak"
    return "low"


def ratio_tier(r: float, weak_ratio: float, strong_ratio: float) -> str:
    if r <= 0 or not math.isfinite(r):
        return "NA"
    if r >= strong_ratio or r <= (1.0 / strong_ratio):
        return "strong"
    if r >= weak_ratio or r <= (1.0 / weak_ratio):
        return "weak"
    return "low"


def call_direction(r: float) -> str:
    if not math.isfinite(r) or r <= 0:
        return "NA"
    if r > 1.0:
        return "gain"
    if r < 1.0:
        return "loss"
    return "neutral"


def weighted_mean(vals, weights):
    v = np.asarray(vals, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not ok.any():
        return float("nan")
    return float(np.sum(v[ok] * w[ok]) / np.sum(w[ok]))


def _mean_int(vals):
    vv = []
    for x in vals:
        if x is None or x is pd.NA:
            continue
        if isinstance(x, float) and (not math.isfinite(x)):
            continue
        try:
            vv.append(float(x))
        except Exception:
            pass
    if not vv:
        return pd.NA
    return float(np.mean(vv))


def _sum_int(vals):
    s = 0
    any_ok = False
    for x in vals:
        if x is None or x is pd.NA:
            continue
        if isinstance(x, float) and (not math.isfinite(x)):
            continue
        s += int(x)
        any_ok = True
    return int(s) if any_ok else pd.NA


def fuse_segments_preserve_support(df: pd.DataFrame, max_gap: int = 0) -> pd.DataFrame:
    """
    Fuse adjacent CNV segments of the same direction (gain or loss) into larger blocks.

    Eligibility (per-segment):
      - ratio_tier != "low"   (i.e., at least weak effect)
      - direction in {"gain","loss"}

    Additional gating (per candidate block):
      - Pass if mean(ratio_z) > 0 OR strict majority of ratio_z > 0
        (If ratio_z missing, fall back to seg_conf_z.)

    Fusion rules:
      - Only fuse within each chromosome
      - Only fuse contiguous runs (by genomic order) where:
          * each member is eligible
          * direction is identical across the run
          * gap between adjacent segments <= max_gap
      - Output includes original rows (unchanged) PLUS fused rows (your downstream code
        should select fused-only if desired).

    Preserved support semantics:
      - Fused CN is probe-weighted (weights = n_probes if available else span length)
      - Fused seg_conf_z uses ONLY the outer boundary z's (leftmost left + rightmost right)
      - Internal boundaries fused-away are summarized into left/right "internal support" columns
        by splitting internal boundaries around midpoint.

    Returns:
      A dataframe including fused rows with extra fused_* and outer_* columns.
    """

    if df.empty:
        return df.copy()

    out = df.copy()

    # Ensure types / helpers
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype(int)
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype(int)
    out["ratio"] = pd.to_numeric(out["ratio"], errors="coerce")
    out["y"] = pd.to_numeric(out["y"], errors="coerce")
    out["cn"] = pd.to_numeric(out["cn"], errors="coerce")

    if "n_probes" in out.columns:
        # may be Int64 with NA; keep numeric view when needed
        nprobes_num = pd.to_numeric(out["n_probes"], errors="coerce")
    else:
        nprobes_num = pd.Series([np.nan] * len(out), index=out.index)

    # ---- Per-segment eligibility: weak+strong effect only; NO keep_suggested filter ----
    eligible = (
        (out.get("ratio_tier", "low") != "low") &
        (out.get("direction", "NA").isin(["gain", "loss"]))
    )



    def _safe_tier_ratio(ratio_val: float, weak_ratio: float, strong_ratio: float) -> str:
        return ratio_tier(float(ratio_val), float(weak_ratio), float(strong_ratio))

    def _safe_tier_z(zval: float, weak_z: float, strong_z: float) -> str:
        return z_tier(float(zval), float(weak_z), float(strong_z))

    def _keep_suggested(rt: str, zt: str, keep_mode: str) -> bool:
        ratio_keep = rt in ["weak", "strong"]
        z_keep = zt in ["weak", "strong"]
        if keep_mode == "ratio_only":
            return ratio_keep
        if keep_mode == "z_only":
            return z_keep
        return ratio_keep and z_keep

    def _call_label(direction: str, rt: str, zt: str) -> str:
        if rt not in ["weak", "strong"]:
            return "no_CNV_effect"
        return f"CNV_{direction}_{rt}Effect_{zt}Support"


    # ratio_z accessor + block gate
    def _get_ratio_z_series(df_block: pd.DataFrame) -> pd.Series:
        if "ratio_z" in df_block.columns:
            z = pd.to_numeric(df_block["ratio_z"], errors="coerce")
        else:
            z = pd.to_numeric(df_block.get("seg_conf_z", np.nan), errors="coerce")
        return z

    def _block_passes_ratio_z(df_block: pd.DataFrame) -> bool:
        z = _get_ratio_z_series(df_block)
        z = z[np.isfinite(z.to_numpy(dtype=float))]

        # Conservative: no usable z => don't fuse
        if len(z) == 0:
            return False

        mean_pass = float(z.mean()) > 0.0
        maj_pass = int((z > 0.0).sum()) > (len(z) / 2.0)  # strict majority
        return bool(mean_pass or maj_pass)

    # Helpers for weights and outer boundaries
    def _weights(df_block: pd.DataFrame) -> np.ndarray:
        wp = pd.to_numeric(df_block.get("n_probes", np.nan), errors="coerce").to_numpy(dtype=float)
        span = (pd.to_numeric(df_block["end"], errors="coerce") - pd.to_numeric(df_block["start"], errors="coerce")).to_numpy(dtype=float)
        w = np.where(np.isfinite(wp) & (wp > 0), wp, span)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        return w

    def _wavg(vals: np.ndarray, w: np.ndarray) -> float:
        m = np.isfinite(vals) & np.isfinite(w) & (w > 0)
        if m.sum() == 0:
            return float("nan")
        return float(np.sum(vals[m] * w[m]) / np.sum(w[m]))

    def _min2(a: float, b: float) -> float:
        aa = a if (a is not None and math.isfinite(a)) else float("nan")
        bb = b if (b is not None and math.isfinite(b)) else float("nan")
        if math.isfinite(aa) and math.isfinite(bb):
            return float(min(aa, bb))
        if math.isfinite(aa):
            return float(aa)
        if math.isfinite(bb):
            return float(bb)
        return float("nan")

    fused_rows = []
    fused_member_idx = set()  # original df indices that were consumed into a fused block

    # Process per chromosome, preserve original ordering
    for chrom, sub in out.groupby("chrom", sort=False):
        sub = sub.sort_values(["start", "end"]).copy()
        sub["_orig_idx"] = sub.index
        sub = sub.reset_index(drop=True)  # <--- critical: now i/j are positional
    
        sub["_eligible"] = (
            eligible.reindex(sub["_orig_idx"])
                    .fillna(False)
                    .to_numpy(dtype=bool)
        )
    
        i = 0
        while i < len(sub):
            if not bool(sub.iloc[i]["_eligible"]):   # <--- iloc everywhere
                i += 1
                continue
            direction = str(sub.iloc[i]["direction"])
            block_idx = [i]
            j = i + 1
            while j < len(sub):
                if not bool(sub.iloc[j]["_eligible"]):
                    break
                if str(sub.iloc[j]["direction"]) != direction:
                    break
                gap = int(sub.iloc[j]["start"]) - int(sub.iloc[j - 1]["end"])
                if gap > int(max_gap):
                    break
                block_idx.append(j)
                j += 1
    
            block = sub.iloc[block_idx].copy()

            # Fuse only if block length > 1 AND passes ratio_z gate
            if len(block) > 1 and _block_passes_ratio_z(block):
                fused_member_idx.update(block["_orig_idx"].tolist())
                w = _weights(block)

                # New interval
                new_start = int(block["start"].min())
                new_end = int(block["end"].max())

                # Weighted effect size
                ratio_vals = pd.to_numeric(block["ratio"], errors="coerce").to_numpy(dtype=float)
                y_vals = pd.to_numeric(block["y"], errors="coerce").to_numpy(dtype=float)
                cn_vals = pd.to_numeric(block["cn"], errors="coerce").to_numpy(dtype=float)
                nprobes_vals = pd.to_numeric(block.get("n_probes", np.nan), errors="coerce").to_numpy(dtype=float)

                fused_ratio = _wavg(ratio_vals, w)
                fused_y = _wavg(y_vals, w)
                fused_cn = _wavg(cn_vals, w)

                # n_probes: sum if available
                fused_nprobes = float(np.nansum(nprobes_vals)) if np.isfinite(nprobes_vals).any() else float("nan")

                # Outer boundaries (names + z)
                left_name = str(block.iloc[0].get("left_boundary_name", "NA"))
                right_name = str(block.iloc[-1].get("right_boundary_name", block.iloc[-1].get("name", "NA")))

                left_z = pd.to_numeric(block.iloc[0].get("left_boundary_z", np.nan), errors="coerce")
                right_z = pd.to_numeric(block.iloc[-1].get("right_boundary_z", np.nan), errors="coerce")

                fused_seg_conf_z = _min2(
                    float(left_z) if math.isfinite(float(left_z)) else float("nan"),
                    float(right_z) if math.isfinite(float(right_z)) else float("nan")
                )

                # ---- Internal boundary support summarization ----
                # Internal boundaries are the boundaries between segments in the block (k-1).
                # Split them around the midpoint; average support to left and right.
                mid = (new_start + new_end) / 2.0

                internal_support = []
                internal_pos = []

                if "right_support_reads" in block.columns:
                    for k in range(len(block) - 1):
                        # boundary between segment k and k+1 is at end of seg k (== start of seg k+1)
                        bpos = int(block.iloc[k]["end"])
                        sval = pd.to_numeric(block.iloc[k].get("right_support_reads", np.nan), errors="coerce")
                        sval_f = float(sval) if math.isfinite(float(sval)) else float("nan")
                        internal_support.append(sval_f)
                        internal_pos.append(bpos)

                # Split into left/right of midpoint
                left_int = [s for s, p in zip(internal_support, internal_pos) if math.isfinite(s) and p <= mid]
                right_int = [s for s, p in zip(internal_support, internal_pos) if math.isfinite(s) and p > mid]

                fused_internal_left_support_mean = float(np.mean(left_int)) if len(left_int) else float("nan")
                fused_internal_right_support_mean = float(np.mean(right_int)) if len(right_int) else float("nan")

                # Outer boundary support if present on endpoints
                outer_left_support = pd.to_numeric(block.iloc[0].get("left_support_reads", np.nan), errors="coerce")
                outer_right_support = pd.to_numeric(block.iloc[-1].get("right_support_reads", np.nan), errors="coerce")

                # Build fused row (start from first row for stable columns)
                base = block.iloc[0].to_dict()

                # ---- fused_segment_names (order-preserving, de-duplicated) ----
                names = [str(x) for x in block["name"].astype(str).tolist()]
                seen = set()
                names_uniq = []
                for nm in names:
                    if nm not in seen:
                        names_uniq.append(nm)
                        seen.add(nm)
                fused_names_str = ",".join(names_uniq)

                # -----------------------------
                # Recompute fusion-derived metrics
                # -----------------------------

                # Pull thresholds/mode from the first row (all rows should share these)
                weak_ratio = float(base.get("weak_ratio", np.nan))
                strong_ratio = float(base.get("strong_ratio", np.nan))
                weak_z = float(base.get("weak_z", np.nan))
                strong_z = float(base.get("strong_z", np.nan))
                keep_mode = str(base.get("keep_mode", "ratio_only"))

                # Determine seg_conf_z policy (do NOT hardcode)
                # If you already store it, use it. Otherwise default conservatively.
                conf_method = str(base.get("conf_method", base.get("confidence_method", "min"))).lower()
                if conf_method not in {"min", "mean", "max"}:
                    conf_method = "min"

                # Outer boundary z values
                left_z_val = float(left_z) if math.isfinite(float(left_z)) else float("nan")
                right_z_val = float(right_z) if math.isfinite(float(right_z)) else float("nan")

                # Apply seg_conf_z policy to the TWO outer boundaries
                z_pair = [z for z in [left_z_val, right_z_val] if math.isfinite(z)]
                if len(z_pair) == 0:
                    fused_seg_conf_z = float("nan")
                elif conf_method == "mean":
                    fused_seg_conf_z = float(np.mean(z_pair))
                elif conf_method == "max":
                    fused_seg_conf_z = float(np.max(z_pair))
                else:
                    fused_seg_conf_z = float(np.min(z_pair))

                # Recompute tiers + keep + call on fused values
                fused_rtier = ratio_tier(float(fused_ratio), weak_ratio, strong_ratio)
                fused_ztier = z_tier(float(fused_seg_conf_z), weak_z, strong_z)
                ratio_keep = fused_rtier in ["weak", "strong"]
                z_keep = fused_ztier in ["weak", "strong"]

                if keep_mode == "ratio_only":
                    fused_keep = ratio_keep
                elif keep_mode == "z_only":
                    fused_keep = z_keep
                else:
                    fused_keep = ratio_keep and z_keep

                if not ratio_keep:
                    fused_call = "no_CNV_effect"
                else:
                    fused_call = f"CNV_{direction}_{fused_rtier}Effect_{fused_ztier}Support"
                fused_call = fused_call + "_FUSED"

                # ratio_z: if your table uses it, define it consistently.
                # Most common: ratio_z == seg_conf_z (z support), or could be z-scored ratio.
                # Here we set ratio_z := fused_seg_conf_z unless you truly have something else.
                fused_ratio_z = fused_seg_conf_z

                # Outer boundary support (keep your existing outer_left_support, outer_right_support)
                out_left_support = float(outer_left_support) if math.isfinite(float(outer_left_support)) else float("nan")
                out_right_support = float(outer_right_support) if math.isfinite(float(outer_right_support)) else float("nan")

                # -----------------------------
                # Update base dict
                # -----------------------------
                base.update({
                    "chrom": str(chrom),
                    "start": int(new_start),
                    "end": int(new_end),

                    # Name: deterministic (you can change later if you want a new fused name)
                    "name": str(block.iloc[0].get("name", "")),

                    "y": float(fused_y),
                    "ratio": float(fused_ratio),
                    "cn": float(fused_cn),

                    "n_probes": int(fused_nprobes) if math.isfinite(fused_nprobes) else pd.NA,
                    "direction": direction,

                    "left_boundary_name": left_name,
                    "right_boundary_name": right_name,
                    "left_boundary_z": left_z_val,
                    "right_boundary_z": right_z_val,
                    "seg_conf_z": float(fused_seg_conf_z),

                    # Recomputed tiers/keep/call
                    "ratio_tier": fused_rtier,
                    "support_tier": fused_ztier,
                    "keep_suggested": bool(fused_keep),
                    "call": fused_call,

                    # If present in your schema
                    "ratio_z": float(fused_ratio_z) if math.isfinite(float(fused_ratio_z)) else float("nan"),

                    # Tags / bookkeeping
                    "fused": True,
                    "fused_n_segments": int(len(names_uniq)),   # de-duplicated count
                    "fused_internal_boundaries": int(max(0, len(names_uniq) - 1)),
                    "fused_segment_names": fused_names_str,

                    # Preserve internal boundary evidence
                    "outer_left_support_reads": out_left_support,
                    "outer_right_support_reads": out_right_support,
                    "internal_left_boundary_support_mean": float(fused_internal_left_support_mean) if math.isfinite(float(fused_internal_left_support_mean)) else float("nan"),
                    "internal_right_boundary_support_mean": float(fused_internal_right_support_mean) if math.isfinite(float(fused_internal_right_support_mean)) else float("nan"),
                })

                # Remove scratch columns so they never appear in output
                base.pop("_orig_idx", None)
                base.pop("_eligible", None)


                # --- recompute length for fused interval ---
                fused_seg_len_bp = int(new_end) - int(new_start)

                # --- recompute / aggregate BAM-derived columns for fused row ---
                # Primary reads: sum across all member segments (counts are per-segment over segment span)
                fused_seg_primary_reads = _sum_int(block.get("seg_primary_reads", [])) if "seg_primary_reads" in block.columns else pd.NA

                # For fused row, left/right boundary window counts should come from the OUTER boundaries only:
                # left boundary = left boundary of FIRST segment in block
                # right boundary = right boundary of LAST segment in block
                fused_left_split_reads  = block.iloc[0].get("left_split_reads", pd.NA)  if "left_split_reads"  in block.columns else pd.NA
                fused_left_disco_reads  = block.iloc[0].get("left_disco_reads", pd.NA)  if "left_disco_reads"  in block.columns else pd.NA
                fused_left_support_reads = (
                    int(fused_left_split_reads) + int(fused_left_disco_reads)
                    if (fused_left_split_reads is not pd.NA and fused_left_disco_reads is not pd.NA)
                    else pd.NA
                )

                fused_right_split_reads = block.iloc[-1].get("right_split_reads", pd.NA) if "right_split_reads" in block.columns else pd.NA
                fused_right_disco_reads = block.iloc[-1].get("right_disco_reads", pd.NA) if "right_disco_reads" in block.columns else pd.NA
                fused_right_support_reads = (
                    int(fused_right_split_reads) + int(fused_right_disco_reads)
                    if (fused_right_split_reads is not pd.NA and fused_right_disco_reads is not pd.NA)
                    else pd.NA
                )

                # outer_left/right_support_reads should mirror the fused left/right support totals
                # (these are your “preserved outer boundary evidence” columns)
                outer_left_support_reads = fused_left_support_reads
                outer_right_support_reads = fused_right_support_reads

                # Internal boundary count should match the *de-duplicated* segment list
                fused_internal_boundaries_n = int(max(0, len(names_uniq) - 1))

                base.update({
                    # length + primary coverage
                    "seg_len_bp": fused_seg_len_bp,
                    "seg_primary_reads": fused_seg_primary_reads,

                    # boundary-window counts at outer boundaries only
                    "left_split_reads": fused_left_split_reads,
                    "left_disco_reads": fused_left_disco_reads,
                    "left_support_reads": fused_left_support_reads,

                    "right_split_reads": fused_right_split_reads,
                    "right_disco_reads": fused_right_disco_reads,
                    "right_support_reads": fused_right_support_reads,

                    # preserve outer boundary support explicitly
                    "outer_left_support_reads": outer_left_support_reads,
                    "outer_right_support_reads": outer_right_support_reads,

                    # keep internal summary + count consistent
                    "fused_internal_boundaries": fused_internal_boundaries_n,
                    "internal_left_boundary_support_mean": float(fused_internal_left_support_mean) if math.isfinite(float(fused_internal_left_support_mean)) else float("nan"),
                    "internal_right_boundary_support_mean": float(fused_internal_right_support_mean) if math.isfinite(float(fused_internal_right_support_mean)) else float("nan"),
                })




                base["fused"] = True
                base["fused_n_segments"] = int(len(names_uniq))
                base["fused_internal_boundaries"] = int(max(0, len(names_uniq) - 1))



                

                fused_rows.append(base)               

            # Move to next chunk
            i = j

    if not fused_rows:
        # Nothing fused; return original df as-is
        out2 = out.copy()
        out2["fused"] = False
        return out2

    fused_df = pd.DataFrame(fused_rows)

    # Defensive de-dupe of fused rows (in case something upstream creates repeats)
    # Keyed on genomic interval + direction + fused segment composition when available.
    dedupe_keys = [c for c in ["chrom", "start", "end", "direction", "fused_segment_names"] if c in fused_df.columns]
    if dedupe_keys:
        fused_df = fused_df.drop_duplicates(subset=dedupe_keys, keep="last").reset_index(drop=True)

    # Mark originals: which ones were consumed into a fused row
    out2 = out.copy()
    if "fused" not in out2.columns:
        out2["fused"] = False

    out2["_fused_member"] = out2.index.map(lambda i: i in fused_member_idx)

    # **Key behavior change**:
    # Drop original segments that are part of any fused block, keep everything else.
    originals_kept = out2.loc[~out2["_fused_member"]].drop(columns=["_fused_member"])

    # Combine: non-member originals + fused rows
    combined = pd.concat([originals_kept, fused_df], ignore_index=True, sort=False)
    # Normalize fused flag: any missing -> False, enforce bool dtype
    if "fused" not in combined.columns:
        combined["fused"] = False
    combined["fused"] = combined["fused"].fillna(False).astype(bool)
    
    # Stable-ish ordering: chrom,start,end,fused (put fused after originals per region)
    combined["start"] = pd.to_numeric(combined["start"], errors="coerce").astype(int)
    combined["end"] = pd.to_numeric(combined["end"], errors="coerce").astype(int)

    if "fused" in combined.columns:
        combined["fused_sort"] = combined["fused"].astype(int)
        combined = combined.sort_values(["chrom", "start", "end", "fused_sort"], kind="mergesort") \
                           .drop(columns=["fused_sort"])
    else:
        combined = combined.sort_values(["chrom", "start", "end"], kind="mergesort")

    return combined




# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="Finalize CNV segments by combining segment ratio (effect size tiers) with boundary z support (confidence tiers), with optional BAM-based read support and optional fusion."
    )
    ap.add_argument("--segments-bed", required=True, help="segments.bed(.gz ok) from Step 2")
    ap.add_argument("--boundary-support-tsv", required=True, help="boundary_support.tsv(.gz ok) from Step 3")
    ap.add_argument("--out-tsv", required=True, help="Output TSV (plain .tsv recommended)")

    ap.add_argument("--out-bedgraph-y", default="", help="Optional bedGraph of y per segment (.gz ok)")
    ap.add_argument("--out-bedgraph-ratio", default="", help="Optional bedGraph of ratio per segment (.gz ok)")

    ap.add_argument("--weak-ratio", type=float, default=1.25)
    ap.add_argument("--strong-ratio", type=float, default=1.5)

    ap.add_argument("--weak-z", type=float, default=1.0)
    ap.add_argument("--strong-z", type=float, default=2.0)

    ap.add_argument("--keep-mode", choices=["ratio_only", "z_only", "both"], default="ratio_only")
    ap.add_argument("--cn-base", type=float, default=1.0)

    # Optional fusion
    ap.add_argument("--fuse", action="store_true", help="Fuse adjacent kept CNVs of same direction (allow weak+strong)")
    ap.add_argument("--fuse-max-gap", type=int, default=0, help="Maximum bp gap allowed between fused segments")

    # Optional BAM-based counts
    ap.add_argument("--primary-bam", default="", help="Primary/aligned BAM for segment read counting (.bai required)")
    ap.add_argument("--split-bam", default="", help="Splitters BAM for boundary read counting (.bai required)")
    ap.add_argument("--disco-bam", default="", help="Discordant BAM for boundary read counting (.bai required)")
    ap.add_argument("--count-flank", type=int, default=1250, help="Flank bp on each side of boundary for read counting")

    ap.add_argument("--count-all-alignments", action="store_true", default=False,
                    help="Count all alignments (otherwise counts read1 only for paired reads)")
    ap.add_argument("--count-drop-dup", action="store_true", default=False,
                    help="Exclude duplicate-marked reads (FLAG 1024) from counts")
    ap.add_argument("--count-include-supp", action="store_true", default=False,
                    help="Include supplementary alignments in counts")
    ap.add_argument("--count-include-secondary", action="store_true", default=False,
                    help="Include secondary alignments in counts")

    args = ap.parse_args()

    if args.weak_ratio <= 1.0 or args.strong_ratio <= 1.0:
        raise RuntimeError("--weak-ratio and --strong-ratio must be > 1.0")
    if args.strong_ratio < args.weak_ratio:
        raise RuntimeError("--strong-ratio must be >= --weak-ratio")
    if args.strong_z < args.weak_z:
        raise RuntimeError("--strong-z must be >= --weak-z")

    seg = read_segments_bed(args.segments_bed)
    b = read_boundary_support(args.boundary_support_tsv)

    # boundary keyed by LEFT segment name
    bz = (
        b.groupby(["chrom", "name"], as_index=False)["combined_z_chr"]
        .max()
        .rename(columns={"combined_z_chr": "boundary_z"})
    )
    zmap = {(r.chrom, r.name): float(r.boundary_z) for r in bz.itertuples(index=False)}

    # Open BAMs once
    primary_bam = _open_bam(args.primary_bam) if args.primary_bam else None
    split_bam = _open_bam(args.split_bam) if args.split_bam else None
    disco_bam = _open_bam(args.disco_bam) if args.disco_bam else None

    read1_only = not bool(args.count_all_alignments)
    drop_dup = bool(args.count_drop_dup)
    include_supp = bool(args.count_include_supp)
    include_secondary = bool(args.count_include_secondary)
    flank = int(args.count_flank)

    out_rows = []

    for chrom, sub in seg.groupby("chrom", sort=False):
        sub = sub.sort_values(["start", "end"]).reset_index(drop=True)

        for i in range(len(sub)):
            r = {k: sub.loc[i, k] for k in sub.columns}
            name = str(r["name"])

            # Right boundary keyed by this segment name
            right_z = zmap.get((chrom, name), float("nan"))

            # Left boundary keyed by previous segment name
            if i == 0:
                left_name = "NA"
                left_z = float("nan")
            else:
                left_name = str(sub.loc[i - 1, "name"])
                left_z = zmap.get((chrom, left_name), float("nan"))

            seg_conf_z = min2(left_z, right_z)

            ratio = float(r["ratio"])
            y = float(r["y"])
            cn = float(args.cn_base) * ratio

            rtier = ratio_tier(ratio, args.weak_ratio, args.strong_ratio)
            ztier = z_tier(seg_conf_z, args.weak_z, args.strong_z)
            direction = call_direction(ratio)

            ratio_keep = (rtier in ["weak", "strong"])
            z_keep = (ztier in ["weak", "strong"])

            if args.keep_mode == "ratio_only":
                keep_suggested = ratio_keep
            elif args.keep_mode == "z_only":
                keep_suggested = z_keep
            else:
                keep_suggested = ratio_keep and z_keep

            if not ratio_keep:
                call = "no_CNV_effect"
            else:
                call = f"CNV_{direction}_{rtier}Effect_{ztier}Support"

            # ---- BAM-derived counts (optional) ----
            seg_len = int(r["end"]) - int(r["start"])
            seg_primary_reads = pd.NA
            left_split_reads = pd.NA
            right_split_reads = pd.NA
            left_disco_reads = pd.NA
            right_disco_reads = pd.NA

            if primary_bam is not None:
                seg_primary_reads = count_reads_in_interval(
                    primary_bam, chrom, int(r["start"]), int(r["end"]),
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=include_supp, include_secondary=include_secondary
                )

            # boundary windows around segment edges
            left_pos = int(r["start"])
            right_pos = int(r["end"])
            l0, l1 = left_pos - flank, left_pos + flank
            r0, r1 = right_pos - flank, right_pos + flank

            if split_bam is not None:
                left_split_reads = count_reads_in_interval(
                    split_bam, chrom, l0, l1,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=include_supp, include_secondary=include_secondary
                )
                right_split_reads = count_reads_in_interval(
                    split_bam, chrom, r0, r1,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=include_supp, include_secondary=include_secondary
                )

            if disco_bam is not None:
                left_disco_reads = count_reads_in_interval(
                    disco_bam, chrom, l0, l1,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=include_supp, include_secondary=include_secondary
                )
                right_disco_reads = count_reads_in_interval(
                    disco_bam, chrom, r0, r1,
                    read1_only=read1_only, drop_dup=drop_dup,
                    include_supp=include_supp, include_secondary=include_secondary
                )

            # Totals if both present
            left_support_reads = pd.NA
            right_support_reads = pd.NA
            if pd.notna(left_split_reads) and pd.notna(left_disco_reads):
                left_support_reads = int(left_split_reads) + int(left_disco_reads)
            if pd.notna(right_split_reads) and pd.notna(right_disco_reads):
                right_support_reads = int(right_split_reads) + int(right_disco_reads)

            r.update(dict(
                cn_base=float(args.cn_base),
                cn=cn,
                direction=direction,

                left_boundary_name=left_name,
                right_boundary_name=name,
                left_boundary_z=left_z,
                right_boundary_z=right_z,
                seg_conf_z=seg_conf_z,

                ratio_tier=rtier,
                support_tier=ztier,
                weak_ratio=args.weak_ratio,
                strong_ratio=args.strong_ratio,
                weak_z=args.weak_z,
                strong_z=args.strong_z,

                keep_mode=args.keep_mode,
                keep_suggested=bool(keep_suggested),
                call=call,

                seg_len_bp=int(seg_len),
                seg_primary_reads=seg_primary_reads,

                left_split_reads=left_split_reads,
                right_split_reads=right_split_reads,
                left_disco_reads=left_disco_reads,
                right_disco_reads=right_disco_reads,
                left_support_reads=left_support_reads,
                right_support_reads=right_support_reads,
            ))

            out_rows.append(r)

    out = pd.DataFrame(out_rows)

    # ------------------------------------------------------------
    # NEW: compute ratio_z (magnitude-based) so fusion gating works
    #      for BOTH gains and losses.
    #
    # ratio_z = robust z-score of |y| per chromosome using MAD:
    #   ratio_z = 0.6745 * (|y| - median(|y|)) / MAD(|y|)
    # ------------------------------------------------------------
    y_abs = pd.to_numeric(out["y"], errors="coerce").abs()

    def _robust_mag_z(g: pd.Series) -> pd.Series:
        v = pd.to_numeric(g, errors="coerce").to_numpy(dtype=float)
        med = np.nanmedian(v)
        mad = np.nanmedian(np.abs(v - med))
        # avoid div-by-zero if chromosome is flat
        denom = mad if (np.isfinite(mad) and mad > 0) else np.nanstd(v)
        if not np.isfinite(denom) or denom <= 0:
            # totally flat: give zeros
            return pd.Series([0.0] * len(v), index=g.index, dtype=float)
        z = 0.6745 * (v - med) / denom
        # magnitude-based score; always >= 0
        z = np.abs(z)
        # replace nans with 0 for safety in gating
        z = np.where(np.isfinite(z), z, 0.0)
        return pd.Series(z, index=g.index, dtype=float)

    out["ratio_z"] = y_abs.groupby(out["chrom"], sort=False).apply(_robust_mag_z).reset_index(level=0, drop=True)

    # Optional fusion (runs AFTER we computed per-segment read columns)
    if args.fuse:
        out = fuse_segments_preserve_support(out, max_gap=int(args.fuse_max_gap))

    # ------------------------------------------------------------
    # NEW: also write fused-only TSV next to the main output if any
    # ------------------------------------------------------------
    out_dir = os.path.dirname(os.path.abspath(args.out_tsv)) or "."
    fused_only_path = os.path.join(out_dir, os.path.basename(args.out_tsv).replace(".tsv", ".fused.tsv").replace(".gz", ".fused.tsv.gz"))

    if "fused" in out.columns and out["fused"].astype(bool).any():
        out_fused_only = out[out["fused"].astype(bool)].copy()
        _atomic_write_no_follow(write_tsv_any, out_fused_only, fused_only_path)
    else:
        # still emit an empty fused file with headers (helps debugging downstream)
        out_empty = out.head(0).copy()
        out_empty["fused"] = pd.Series([], dtype=bool)
        _atomic_write_no_follow(write_tsv_any, out_empty, fused_only_path)

    
    # Optional fusion (runs AFTER we computed per-segment read columns)
    if args.fuse:
        out = fuse_segments_preserve_support(out, max_gap=int(args.fuse_max_gap))

    # Stable column order (preserve any extras after these)
    col_order = [
        "chrom", "start", "end", "name",
        "y", "ratio", "n_probes",
        "cn_base", "cn", "direction",
        "left_boundary_name", "right_boundary_name",
        "left_boundary_z", "right_boundary_z", "seg_conf_z",
        "ratio_tier", "support_tier",
        "weak_ratio", "strong_ratio", "weak_z", "strong_z",
        "keep_mode", "keep_suggested", "call",

        # BAM-derived support
        "seg_len_bp", "seg_primary_reads",
        "left_split_reads", "left_disco_reads", "left_support_reads",
        "right_split_reads", "right_disco_reads", "right_support_reads",

        # Fusion-preserved extra support (if present)
        "outer_left_split_reads", "outer_left_disco_reads", "outer_left_support_reads",
        "outer_right_split_reads", "outer_right_disco_reads", "outer_right_support_reads",
        "fused_internal_boundaries", "fused_left_boundaries_n", "fused_right_boundaries_n",
    ]
    out = out[[c for c in col_order if c in out.columns] + [c for c in out.columns if c not in col_order]]

    _atomic_write_no_follow(write_tsv_any, out, args.out_tsv)

    if args.out_bedgraph_y:
        bg = out[["chrom", "start", "end"]].copy()
        bg["value"] = pd.to_numeric(out["y"], errors="coerce").astype(float)
        _atomic_write_no_follow(write_bedgraph_any, bg, args.out_bedgraph_y)

    if args.out_bedgraph_ratio:
        bg = out[["chrom", "start", "end"]].copy()
        bg["value"] = pd.to_numeric(out["ratio"], errors="coerce").astype(float)
        _atomic_write_no_follow(write_bedgraph_any, bg, args.out_bedgraph_ratio)

    # Close BAM handles
    for bam in (primary_bam, split_bam, disco_bam):
        try:
            if bam is not None:
                bam.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
