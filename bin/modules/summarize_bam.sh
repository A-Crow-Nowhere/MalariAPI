#!/usr/bin/env bash
# MAPI module template (validator-friendly)
# Recommended final location: bin/modules/<name>.sh
# Expected metadata: bin/modules/yaml/<name>.yml
set -euo pipefail

# Identify self early so Usage can print even if lock isn’t reachable yet
self="$(basename "$0")"
mod="${self%.sh}"

# ------------------- Usage (designed to satisfy validator) -------------------
usage(){ cat <<EOF
Usage: $self --in <reads.fastq> [--in2 <reads2.fastq>] [--threads N] [-- ...tool flags...]
Aliases: --in1 ≡ --in ≡ -1 ; --in2 ≡ -2 ; --threads ≡ -t
Notes:
  • Outputs are ALWAYS written to: {sample_prefix}_mapi-out/
  • Filenames MUST be: {sample_prefix}.<TAG>.<ext>
  • Default TAG for this module: '${DEFAULT_TAG}'
  • This is a template; real modules should live under bin/modules/<name>.sh
  • Metadata YAML is expected at bin/modules/yaml/<name>.yml unless overridden by \$MAPI_MODULE_YAML

EOF
}

# Hardcoded module tag (documented in Usage above)
DEFAULT_TAG="${DEFAULT_TAG:-mapiTool}"

# Let --help/-h succeed before sourcing anything (so validator can read Usage)
if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then
  usage; exit 0
fi

# ---------------- Enforce MAPI miniconda (validator-friendly) ----------------
# The FIRST line below MUST remain exactly as shown to satisfy E105.
# shellcheck disable=SC1090
source "$(dirname "$0")/../_mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../bin/_mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../._mapi_conda_lock.sh" \
  || source "$(dirname "$0")/../bin/._mapi_conda_lock.sh" \
  || { echo "[ERROR] Could not source _mapi_conda_lock(.sh) from parent or bin/." >&2; exit 3; }

# ----------------------- Repo root & metadata discovery ----------------------
REPO_ROOT="${MAPI_REPO_ROOT:-"$(cd "$(dirname "$0")/../.." 2>/dev/null || pwd)"}"
meta="${MAPI_MODULE_YAML:-"$REPO_ROOT/bin/modules/yaml/$mod.yml"}"

if [[ ! -f "$meta" ]]; then
  echo "[warn] Missing metadata YAML (expected: $meta). Template will still run." >&2
fi

# Optionally activate env from YAML, if present (no-op when mapi already runs env)
if [[ -f "$meta" ]]; then
  env_yaml="$(awk -F: '/^env:/ {f=1} f&&/conda:/ {gsub(/[\"'\'' ]/,""); print $2; exit}' "$meta" || true)"
  if [[ -n "${env_yaml:-}" && -f "$REPO_ROOT/$env_yaml" ]]; then
    env_name="$(basename "${env_yaml%.yml}")"
    conda env list >/dev/null 2>&1 && conda activate "$env_name" || true
  fi
fi

# ----------------------------- Defaults & CLI --------------------------------
THREADS="${THREADS:-8}"
# OUTDIR is determined by sample_prefix; user-supplied --out is ignored (warned) to enforce the rule.
OUTDIR=""     # placeholder; set after parsing IN1
IN1="" IN2=""

EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in1|--in|-1)            IN1="$2"; shift 2;;
    --in2|-2)                 IN2="$2"; shift 2;;
    --out|--outdir|-o)        echo "[warn] --out/--outdir is ignored; outputs go to {sample_prefix}_mapi-out/ by spec" >&2; shift 2;;
    --threads|-t)             THREADS="$2"; shift 2;;
    --)                       shift; EXTRA+=("$@"); break;;
    -h|--help)                usage; exit 0;;
    *) echo "Unknown: $1" >&2; usage; exit 2;;
  esac
done

# Validator-visible guards (documented in Usage above)
[[ -n "$IN1" ]] || { echo "[error] missing --in/--in1" >&2; usage; exit 2; }

# Derive sample_prefix from IN1 (basename without first extension)
bn="$(basename "$IN1")"
sample_prefix="${bn%%.*}"

# Enforce outdir & filename conventions
OUTDIR="${sample_prefix}_mapi-out"
mkdir -p "$OUTDIR"

# Helper to produce tagged filenames
# Usage: out_path <ext>  -> prints  ${OUTDIR}/${sample_prefix}.${DEFAULT_TAG}.<ext>
out_path(){ printf "%s/%s.%s.%s\n" "$OUTDIR" "$sample_prefix" "$DEFAULT_TAG" "$1"; }

# ------------------------- IMPLEMENT YOUR TOOL HERE --------------------------
# Example: run an inline Python tool that summarizes BAM/CRAM to a report.
# Mapping:
#   • Wrapper --in/--in1 → Python positional input
#   • Wrapper --out (directory) → Python -o <OUTDIR>/<mod>_summary.txt (default)
#   • Any flags after "--" are forwarded to Python (e.g., -r ref.fa --primary-only)

# Choose Python runner (honor MAPI_ENV if set)
if [[ -n "${MAPI_ENV:-}" ]]; then
  RUN_PY=(conda run -n "${MAPI_ENV}" python3)
else
  RUN_PY=(python3)
fi

# Default output file inside OUTDIR; users can override by passing "-o ..." after "--"
py_out="${OUTDIR}/${mod}_summary.txt"

# Execute embedded Python; order matters: our default "-o" comes first so user "-o" after "--" can override it.
exec "${RUN_PY[@]}" - "$IN1" -o "$py_out" "${EXTRA[@]}" <<'PYCODE'
#!/usr/bin/env python3
import argparse, collections, sys, os
try:
    import pysam
except ModuleNotFoundError:
    sys.stderr.write(
        "ERROR: pysam is required.\n"
        "Install via conda (recommended):   conda install -c bioconda pysam\n"
        "or via pip:                         python -m pip install pysam\n"
    )
    sys.exit(1)

def summarize(in_path, out_path, reference=None, primary_only=False, min_mapq=None):
    mode = "r" if in_path.endswith(".sam") else "rb"
    kwargs = {}
    if in_path.lower().endswith(".cram"):
        if not reference:
            sys.exit("ERROR: CRAM input requires --reference FASTA")
        kwargs["reference_filename"] = reference
    try:
        bam = pysam.AlignmentFile(in_path, mode, **kwargs)
    except Exception as e:
        sys.exit(f"ERROR: failed to open {in_path}: {e}")

    counts = collections.Counter()
    soft_bases = 0
    hard_bases = 0

    for aln in bam.fetch(until_eof=True):
        if primary_only and (aln.is_unmapped or aln.is_secondary or aln.is_supplementary):
            continue
        if (min_mapq is not None) and (not aln.is_unmapped) and (aln.mapping_quality < min_mapq):
            continue

        if aln.is_unmapped: counts["unmapped"] += 1
        else:
            counts["mapped"] += 1
            if not aln.is_secondary and not aln.is_supplementary:
                counts["primary_mapped"] += 1
        if aln.mate_is_unmapped: counts["mate_unmapped"] += 1
        if aln.is_reverse:       counts["reverse_strand"] += 1
        if aln.is_secondary:     counts["secondary"] += 1
        if aln.is_qcfail:        counts["qc_fail"] += 1
        if aln.is_duplicate:     counts["duplicate"] += 1
        if aln.is_supplementary: counts["supplementary"] += 1
        if aln.is_paired:        counts["paired"] += 1
        if aln.is_proper_pair:   counts["proper_pair"] += 1
        if aln.is_read1:         counts["read1"] += 1
        if aln.is_read2:         counts["read2"] += 1

        hadS = hadH = False
        if aln.cigartuples:
            for op, length in aln.cigartuples:  # 4=S, 5=H
                if op == 4: soft_bases += length; hadS = True
                elif op == 5: hard_bases += length; hadH = True
        if hadS: counts["soft_clipped_any"] += 1
        if hadH: counts["hard_clipped_any"] += 1

        if aln.is_supplementary or aln.has_tag("SA"):
            counts["split_read_like"] += 1

    bam.close()

    lines = []
    lines.append(f"# BWA summary for: {os.path.abspath(in_path)}")
    lines.append("# Counts are read-level unless noted.")
    order = [
        "mapped","primary_mapped","unmapped","mate_unmapped",
        "paired","proper_pair","read1","read2",
        "reverse_strand","secondary","supplementary","qc_fail","duplicate",
        "soft_clipped_any","hard_clipped_any","split_read_like"
    ]
    for k in order:
        lines.append(f"{k}\t{counts.get(k,0)}")
    lines.append(f"soft_clipped_bases\t{soft_bases}")
    lines.append(f"hard_clipped_bases\t{hard_bases}")

    try:
        with open(out_path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
    except Exception as e:
        sys.exit(f"ERROR: failed to write {out_path}: {e}")
    return out_path

def build_parser():
    p = argparse.ArgumentParser(
        prog="bwa-summary",
        description="Summarize SAM/BAM/CRAM flags and clipping into a text report."
    )
    p.add_argument("in_bam", help="Input SAM/BAM/CRAM file")
    p.add_argument("-o","--out", default="bwa_summary.txt",
                   help="Output text file (default: bwa_summary.txt)")
    p.add_argument("-r","--reference", help="Reference FASTA (required for CRAM)")
    p.add_argument("--primary-only", action="store_true",
                   help="Count only primary mapped reads")
    p.add_argument("--min-mapq", type=int, default=None,
                   help="Minimum MAPQ to include")
    return p

def main(argv=None):
    args = build_parser().parse_args(argv)
    outp = summarize(args.in_bam, args.out, args.reference, args.primary_only, args.min_mapq)
    print(f"Wrote summary: {outp}")

if __name__ == "__main__":
    main()
PYCODE


# Safe demo outputs to illustrate the contract:
echo "[$mod] meta=${meta:-none} env=${env_name:-none} threads=$THREADS in1=$IN1 in2=${IN2:-none} out=$OUTDIR tag=$DEFAULT_TAG" \
  | tee "$(out_path log)"
echo "TEMPLATE OUTPUT for $mod" > "$(out_path out)"
