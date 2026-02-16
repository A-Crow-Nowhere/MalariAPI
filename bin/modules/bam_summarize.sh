#!/usr/bin/env bash
# MAPI module template (validator-friendly)
# Final location: bin/modules/bam_summarize.sh
set -euo pipefail

self="$(basename "$0")"
mod="${self%.sh}"

# ------------------- Usage (designed to satisfy validator) -------------------
usage(){ cat <<EOF
Usage:
  $self --dir <DIR_WITH_BAMS> [options]
  $self --bam <A.bam> [--bam <B.bam> ...] [options]

Notes:
  • This module writes ONE TSV row per BAM/CRAM.
  • TSV metrics are derived from samtools + mosdepth only.
  • Qualimap is visualization-only; it does NOT feed TSV metrics.
  • This module deliberately does NOT allow upstream filtering of the BAM stream.
    (No "samtools-extra", no MAPQ filtering knobs here.)

Inputs:
  --dir DIR
      Scan DIR for *.bam (and *.cram if present).

  --bam PATH
      Provide one or more BAM/CRAM paths explicitly.

Outputs:
  --tsv PATH
      Output TSV file. Default: bam_summarize.tsv (in CWD)

  --out-dir DIR
      Output directory for side artifacts (CV/qualimap). Default: bam_summarize_out/

Options:
  --threads N
      Threads for samtools/mosdepth where applicable. Default: 8

Coverage CV (optional):
  --cv
      Compute coefficient of variation (CV) over non-overlapping bins using mosdepth:
        cv_depth_bins = (sd(depth_bin_means) / mean(depth_bin_means)) * 100
      Results appear in TSV columns:
        depth_mean_bins, depth_sd_bins, cv_depth_bins
  --binsize N
      Bin size for mosdepth windows. Default: 20000

Qualimap (optional visualization):
  --qualimap
      Run qualimap bamqc per BAM. Outputs go under:
        <out-dir>/qualimap/<sample>/
  --qualimap-threads N
      Threads for qualimap (-nt). Default: same as --threads
  --qualimap-nreads N
      Run qualimap on an approximate subsample of N reads by creating a temp BAM.
      (Internally computes fraction = N / total_records and uses samtools view -s.)
  --qualimap-subsample FRACTION
      Alternative to --qualimap-nreads. FRACTION in (0,1].
  --qualimap-seed INT
      Seed used for subsampling (samtools -s). Default: 42

EOF
}

# Let --help/-h succeed before sourcing anything (validator-friendly)
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
  echo "[warn] Missing metadata YAML (expected: $meta). Module will still run." >&2
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
QUALI_THREADS=""
BINSIZE=20000
DO_CV=0
DO_QUALIMAP=0

TSV_OUT="bam_summarize.tsv"
OUTDIR="bam_summarize_out"

QUALI_NREADS=""
QUALI_FRAC=""
QUALI_SEED=42

DIR=""
BAMS=()

die(){ echo "[error] $*" >&2; exit 2; }
have(){ command -v "$1" >/dev/null 2>&1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) DIR="${2:-}"; shift 2;;
    --bam) BAMS+=("${2:-}"); shift 2;;

    --tsv) TSV_OUT="${2:-}"; shift 2;;
    --out-dir|--outdir|-o) OUTDIR="${2:-}"; shift 2;;

    --threads|-t) THREADS="${2:-}"; shift 2;;

    --cv) DO_CV=1; shift 1;;
    --binsize) BINSIZE="${2:-}"; shift 2;;

    --qualimap) DO_QUALIMAP=1; shift 1;;
    --qualimap-threads) QUALI_THREADS="${2:-}"; shift 2;;
    --qualimap-nreads) QUALI_NREADS="${2:-}"; shift 2;;
    --qualimap-subsample) QUALI_FRAC="${2:-}"; shift 2;;
    --qualimap-seed) QUALI_SEED="${2:-}"; shift 2;;

    -h|--help) usage; exit 0;;
    *)
      echo "Unknown: $1" >&2
      usage; exit 2;;
  esac
done

[[ -n "${QUALI_THREADS:-}" ]] || QUALI_THREADS="$THREADS"

# Guards
[[ -n "$DIR" || ${#BAMS[@]} -gt 0 ]] || { usage; die "Provide --dir or at least one --bam"; }
[[ "$BINSIZE" =~ ^[0-9]+$ ]] || die "--binsize must be integer"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be integer"
[[ "$QUALI_THREADS" =~ ^[0-9]+$ ]] || die "--qualimap-threads must be integer"
[[ "$QUALI_SEED" =~ ^[0-9]+$ ]] || die "--qualimap-seed must be integer"

if [[ -n "$QUALI_NREADS" && -n "$QUALI_FRAC" ]]; then
  die "Use only one of --qualimap-nreads or --qualimap-subsample"
fi
if [[ -n "$QUALI_NREADS" ]]; then
  [[ "$QUALI_NREADS" =~ ^[0-9]+$ ]] || die "--qualimap-nreads must be integer"
fi
if [[ -n "$QUALI_FRAC" ]]; then
  # accept formats like 0.25 or 1
  [[ "$QUALI_FRAC" =~ ^0\.[0-9]+$|^1(\.0+)?$ ]] || die "--qualimap-subsample must be in (0,1]"
fi

have samtools || die "samtools not found in PATH"
if [[ "$DO_CV" -eq 1 ]]; then
  have mosdepth || die "mosdepth not found in PATH (required for --cv)"
fi
if [[ "$DO_QUALIMAP" -eq 1 ]]; then
  have qualimap || die "qualimap not found in PATH (required for --qualimap)"
fi

mkdir -p "$OUTDIR"
mkdir -p "$(dirname "$TSV_OUT")" 2>/dev/null || true

# ---------------------------- Collect BAM list --------------------------------
if [[ -n "$DIR" ]]; then
  [[ -d "$DIR" ]] || die "--dir not found: $DIR"
  while IFS= read -r -d '' f; do BAMS+=("$f"); done < <(find "$DIR" -maxdepth 1 -type f \( -name "*.bam" -o -name "*.cram" \) -print0 | sort -z)
fi

# De-dupe while preserving order
if [[ ${#BAMS[@]} -gt 1 ]]; then
  mapfile -t BAMS < <(printf "%s\n" "${BAMS[@]}" | awk '!seen[$0]++')
fi
[[ ${#BAMS[@]} -gt 0 ]] || die "No BAM/CRAM files found."

# ---------------------------- Helpers -----------------------------------------
# Ensure index exists (bam.bai or cram.crai)
ensure_index() {
  local bam="$1"
  if [[ "$bam" == *.cram || "$bam" == *.CRAM ]]; then
    [[ -s "${bam}.crai" ]] || samtools index -@ "$THREADS" "$bam"
  else
    # samtools index may write bam.bai or <bam>.bai depending on version; accept either
    if [[ -s "${bam}.bai" ]]; then
      return 0
    fi
    if [[ -s "$(dirname "$bam")/$(basename "$bam" .bam).bai" ]]; then
      return 0
    fi
    samtools index -@ "$THREADS" "$bam"
  fi
}

# Parse samtools flagstat into simple scalar metrics
# Outputs tab-separated:
# total, secondary, supplementary, duplicates, mapped, paired, properly_paired, singletons, mate_mapped_different_chr
flagstat_metrics() {
  local bam="$1"
  local txt
  txt="$(samtools flagstat -@ "$THREADS" "$bam")"

  # shell-safe parsing with awk matching line starts
  awk -v RS='\n' '
    BEGIN{
      total=secondary=supp=dups=mapped=paired=proper=single=diffchr=0
    }
    $0 ~ / in total/ { total=$1 }
    $0 ~ / secondary/ { secondary=$1 }
    $0 ~ / supplementary/ { supp=$1 }
    $0 ~ / duplicates/ { dups=$1 }
    $0 ~ / mapped \(/ { mapped=$1 }
    $0 ~ / paired in sequencing/ { paired=$1 }
    $0 ~ / properly paired/ { proper=$1 }
    $0 ~ / singletons/ { single=$1 }
    $0 ~ / with mate mapped to a different chr/ && $0 !~ /mapQ/ { diffchr=$1 }
    END{
      printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", total,secondary,supp,dups,mapped,paired,proper,single,diffchr
    }
  ' <<<"$txt"
}

# Parse samtools stats "SN" lines for a few stable fields
# Outputs tab-separated:
# reads, bases, average_length, insert_size_avg (if present else NA)
stats_metrics() {
  local bam="$1"
  local txt
  txt="$(samtools stats -@ "$THREADS" "$bam" 2>/dev/null | awk -F'\t' '$1=="SN"{print $2"\t"$3}')"

  # helper to pull a key by exact label
  get_sn() {
    local key="$1"
    awk -v k="$key" -F'\t' '$1==k{print $2; found=1} END{if(!found) print "NA"}' <<<"$txt"
  }

  local reads bases avglen isize
  reads="$(get_sn "raw total sequences:")"
  bases="$(get_sn "total length:")"
  avglen="$(get_sn "average length:")"
  isize="$(get_sn "insert size average:")"

  printf "%s\t%s\t%s\t%s\n" "$reads" "$bases" "$avglen" "$isize"
}

# Compute CV over mosdepth bins: mean, sd, cv%
cv_from_mosdepth_bins() {
  local bam="$1"
  local prefix="$2"   # full path prefix inside OUTDIR
  local binsize="$3"

  # mosdepth regions output: <prefix>.regions.bed.gz with columns:
  # chrom  start  end  mean_depth
  mosdepth --threads "$THREADS" --by "$binsize" --no-per-base "$prefix" "$bam" >/dev/null 2>&1

  local regions_gz="${prefix}.regions.bed.gz"
  [[ -s "$regions_gz" ]] || { echo -e "NA\tNA\tNA"; return 0; }

  # Compute mean and sample sd over depth column (4th col)
  # sd = sqrt( (sum(x^2) - n*mean^2) / (n-1) ) for n>1
  # cv% = sd/mean*100
  zcat -f "$regions_gz" | awk '
    {
      x=$4+0
      n++
      s+=x
      ss+=x*x
    }
    END{
      if(n==0){ print "NA\tNA\tNA"; exit }
      mean=s/n
      if(n>1){
        var=(ss - n*mean*mean)/(n-1)
        if(var<0) var=0
        sd=sqrt(var)
      } else {
        sd=0
      }
      if(mean>0){
        cv=(sd/mean)*100
        printf "%.6f\t%.6f\t%.6f\n", mean, sd, cv
      } else {
        printf "%.6f\t%.6f\tNA\n", mean, sd
      }
    }'
}

# Make a safe sample label from filename
sample_from_path() {
  local p="$1"
  local b
  b="$(basename "$p")"
  b="${b%.bam}"
  b="${b%.cram}"
  echo "$b"
}



run_qualimap() {
  local bam="$1"
  local sample="$2"
  local qdir="$3"

  mkdir -p "$qdir"
  # qualimap bamqc writes HTML/PDF/plots under outdir
  # Use --java-mem-size conservatively; user can adjust by env var if needed.
  qualimap bamqc \
    -bam "$bam" \
    -outdir "$qdir" \
    -nt "$QUALI_THREADS" \
    >/dev/null 2>&1 || {
      echo "[warn] qualimap failed for $sample ($bam); continuing." >&2
      return 0
    }
}

# ---------------------------- Write TSV header --------------------------------
# Keep columns stable + explicit.
{
  printf "sample\tpath\tbytes\tidx_present\t"
  printf "flagstat_total\tflagstat_secondary\tflagstat_supplementary\tflagstat_duplicates\tflagstat_mapped\tflagstat_paired\tflagstat_properly_paired\tflagstat_singletons\tflagstat_mate_mapped_diffchr\t"
  printf "stats_raw_total_sequences\tstats_total_length\tstats_avg_length\tstats_insert_size_avg\t"
  printf "depth_mean_bins\tdepth_sd_bins\tcv_depth_bins\t"
  printf "qualimap_outdir\n"
} > "$TSV_OUT"

# ---------------------------- Main loop ---------------------------------------
for bam in "${BAMS[@]}"; do
  [[ -f "$bam" ]] || { echo "[warn] missing file (skip): $bam" >&2; continue; }

  sample="$(sample_from_path "$bam")"
  bytes="$(wc -c < "$bam" | awk '{print $1}')"

  # Index
  idx_present=0
  if [[ "$bam" == *.cram || "$bam" == *.CRAM ]]; then
    [[ -s "${bam}.crai" ]] && idx_present=1
  else
    [[ -s "${bam}.bai" ]] && idx_present=1
    [[ -s "$(dirname "$bam")/$(basename "$bam" .bam).bai" ]] && idx_present=1
  fi
  if [[ "$idx_present" -eq 0 ]]; then
    ensure_index "$bam"
    # re-check
    if [[ "$bam" == *.cram || "$bam" == *.CRAM ]]; then
      [[ -s "${bam}.crai" ]] && idx_present=1
    else
      [[ -s "${bam}.bai" ]] && idx_present=1
      [[ -s "$(dirname "$bam")/$(basename "$bam" .bam).bai" ]] && idx_present=1
    fi
  fi

  # Metrics
  IFS=$'\t' read -r fs_total fs_secondary fs_supp fs_dups fs_mapped fs_paired fs_proper fs_single fs_diffchr < <(flagstat_metrics "$bam")
  IFS=$'\t' read -r st_reads st_bases st_avglen st_isize < <(stats_metrics "$bam")

  # CV
  depth_mean="NA"; depth_sd="NA"; cv="NA"
  if [[ "$DO_CV" -eq 1 ]]; then
    cv_prefix="$OUTDIR/cv/${sample}.bins${BINSIZE}"
    mkdir -p "$OUTDIR/cv"
    IFS=$'\t' read -r depth_mean depth_sd cv < <(cv_from_mosdepth_bins "$bam" "$cv_prefix" "$BINSIZE")
  fi

  # Qualimap
  qualimap_out="NA"
  if [[ "$DO_QUALIMAP" -eq 1 ]]; then
    mkdir -p "$OUTDIR/qualimap"
    qdir="$OUTDIR/qualimap/$sample"
    qualimap_out="$qdir"

    qbam="$bam"
    tmpqdir="$OUTDIR/.tmp_qualimap"
    mkdir -p "$tmpqdir"

   fi

    run_qualimap "$qbam" "$sample" "$qdir"
  fi

  # Row
  {
    printf "%s\t%s\t%s\t%d\t" "$sample" "$bam" "$bytes" "$idx_present"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" \
      "$fs_total" "$fs_secondary" "$fs_supp" "$fs_dups" "$fs_mapped" "$fs_paired" "$fs_proper" "$fs_single" "$fs_diffchr"
    printf "%s\t%s\t%s\t%s\t" "$st_reads" "$st_bases" "$st_avglen" "$st_isize"
    printf "%s\t%s\t%s\t" "$depth_mean" "$depth_sd" "$cv"
    printf "%s\n" "$qualimap_out"
  } >> "$TSV_OUT"

done

echo "[bam_summarize] wrote: $TSV_OUT"
echo "[bam_summarize] outdir: $OUTDIR"
