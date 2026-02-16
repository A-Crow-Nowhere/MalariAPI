#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# crows_to_svcrows.sh
# Convert segments.calls.tsv -> SVCROWS TSV
# Optional: prepend known features from TSV
# ============================================================

say() { echo "[crows_to_svcrows] $*" >&2; }
die() { say "ERROR: $*"; exit 1; }

usage() {
  cat >&2 <<'EOF'
Usage:
  crows_to_svcrows.sh -i <in.tsv> -o <out.tsv> [--known <known.tsv>] [--known-n 1] [--known-var3 <STRING>]
  crows_to_svcrows.sh -i <in.tsv> -O <outdir>  [--known <known.tsv>] [--known-n 1] [--known-var3 <STRING>]

Main input (segments.calls-style) required columns:
  chrom start end name ratio cn direction ratio_tier seg_conf_z

Main filter:
  - ratio_tier in {weak,strong} (case-insensitive)
  - seg_conf_z numeric and >= 0
  - handles blanks/NA/NaN by skipping those rows

Known features (optional):
  --known <FILE>      Can be EITHER:
                      (A) SVCROWS-format (same as output header):
                          Chr Start End Length Type ID Var1 Var2 Var3 IsKnown QScore NumReads
                          -> first N data rows are prepended as-is (but Var3/IsKnown can be overridden)
                      OR
                      (B) segments.calls-style:
                          chrom start end name ratio cn direction ...
                          -> converted and prepended (no filters applied)
  --known-n <INT>     Number of known rows to prepend (default 1; takes first N data lines)
  --known-var3 <STR>  Var3 value for known rows (default "KNOWN")

Output columns (SVCROWS):
  Chr Start End Length Type ID Var1 Var2 Var3 IsKnown QScore NumReads

Mapping for main rows:
  Chr=chrom Start=start End=end Length=End-Start
  Type=direction ID=1..N Var1=ratio Var2=name Var3=SAMPLE(from filename) IsKnown=FALSE
  QScore=cn*Length NumReads=Length

Weighted CN recovery after merges (if sums):
  cn_weighted = QScore_merged / NumReads_merged
EOF
}

in_tsv=""
out_tsv=""
outdir=""

known_tsv=""
known_n=1
known_var3="KNOWN"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in) in_tsv="${2:-}"; shift 2;;
    -o|--out) out_tsv="${2:-}"; shift 2;;
    -O|--outdir) outdir="${2:-}"; shift 2;;
    --known) known_tsv="${2:-}"; shift 2;;
    --known-n) known_n="${2:-}"; shift 2;;
    --known-var3) known_var3="${2:-}"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (try --help)";;
  esac
done

[[ -n "$in_tsv" ]] || die "Missing -i/--in"
[[ -f "$in_tsv" ]] || die "Input not found: $in_tsv"

if [[ -n "$outdir" && -n "$out_tsv" ]]; then
  die "Use either -o OR -O, not both"
fi

if [[ -n "$known_tsv" ]]; then
  [[ -f "$known_tsv" ]] || die "Known TSV not found: $known_tsv"
  [[ "$known_n" =~ ^[0-9]+$ ]] || die "--known-n must be an integer"
  [[ "$known_n" -ge 1 ]] || die "--known-n must be >= 1"
fi

# Derive SAMPLE from input filename prefix
bn="$(basename "$in_tsv")"
sample="${bn%%.segments.calls.tsv}"
if [[ "$sample" == "$bn" ]]; then
  die "Input filename must look like <SAMPLE>.segments.calls.tsv (got: $bn)"
fi

if [[ -n "$outdir" ]]; then
  mkdir -p "$outdir"
  stem="${bn%.tsv}"
  out_tsv="$outdir/${stem}.svcrows.tsv"
else
  [[ -n "$out_tsv" ]] || die "Missing -o/--out (or use -O/--outdir)"
  mkdir -p "$(dirname "$out_tsv")"
fi

say "Input : $in_tsv"
say "Sample: $sample"
say "Output: $out_tsv"
if [[ -n "$known_tsv" ]]; then
  say "Known : $known_tsv (prepend first $known_n row(s), Var3=$known_var3)"
fi

awk -v OFS='\t' -v SAMPLE="$sample" -v KNOWN_TSV="$known_tsv" -v KNOWN_N="$known_n" -v KNOWN_VAR3="$known_var3" '
  BEGIN { FS = "[ \t]+"; }

  function norm(s) { gsub(/\r/, "", s); return s; }
  function lower(s){ s=norm(s); return tolower(s); }
  function is_blank(s){ s=lower(s); return (s=="" || s=="na" || s=="nan" || s=="null"); }
  function is_num(s){ s=norm(s); return (s ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/); }

  # ------------------------------------------------------------
  # emit_known:
  # - autodetect header type:
  #   * SVCROWS header: has "Chr" and "QScore" and "NumReads"
  #   * segments header: has "chrom" and "cn" and "direction"
  # - prepend first KNOWN_N data rows
  # - known rows get new sequential ID in final output
  # - known rows force Var3=KNOWN_VAR3 and IsKnown=TRUE
  # ------------------------------------------------------------
  function emit_known(file, nmax,   line, nf, i, count,
                      has_svc, has_seg,
                      # svc cols
                      chr, st, en, len, typ, var1, var2, qscore, numreads,
                      # seg cols
                      nm, rat, cnr) {

    if (file == "") return;

    if ((getline line < file) <= 0) {
      print "[crows_to_svcrows] ERROR: could not read known TSV: " file > "/dev/stderr";
      exit 2;
    }
    gsub(/\r/, "", line);
    nf = split(line, H, FS);
    delete kidx;
    for (i=1; i<=nf; i++) kidx[H[i]] = i;

    has_svc = (("Chr" in kidx) && ("Start" in kidx) && ("End" in kidx) && ("QScore" in kidx) && ("NumReads" in kidx));
    has_seg = (("chrom" in kidx) && ("start" in kidx) && ("end" in kidx) && ("cn" in kidx) && ("direction" in kidx));

    if (!has_svc && !has_seg) {
      print "[crows_to_svcrows] ERROR: known TSV header not recognized. Need either SVCROWS header (Chr..QScore..NumReads) or segments header (chrom..cn..direction)." > "/dev/stderr";
      exit 2;
    }

    count = 0;
    while ((getline line < file) > 0) {
      gsub(/\r/, "", line);
      if (line ~ /^[ \t]*$/) continue;
      nf = split(line, F, FS);

      if (has_svc) {
        # already in svcrows format
        chr = F[kidx["Chr"]];
        st  = F[kidx["Start"]];
        en  = F[kidx["End"]];
        len = F[kidx["Length"]];
        typ = F[kidx["Type"]];
        var1 = F[kidx["Var1"]];
        var2 = F[kidx["Var2"]];
        qscore = F[kidx["QScore"]];
        numreads = F[kidx["NumReads"]];

        if (is_blank(chr) || !is_num(st) || !is_num(en) || !is_num(qscore) || !is_num(numreads)) continue;

        id++;
        # force Var3 + IsKnown, keep everything else
        print chr, (st+0), (en+0), (len+0), typ, id, var1, var2, KNOWN_VAR3, "TRUE", (qscore+0), (numreads+0);

      } else {
        # segments.calls style: convert minimal required cols (no filters applied)
        chr = F[kidx["chrom"]];
        st  = F[kidx["start"]];
        en  = F[kidx["end"]];
        nm  = ("name" in kidx ? F[kidx["name"]] : "");
        rat = ("ratio" in kidx ? F[kidx["ratio"]] : "");
        cnr = F[kidx["cn"]];
        typ = F[kidx["direction"]];

        if (is_blank(chr) || !is_num(st) || !is_num(en) || !is_num(cnr)) continue;

        st += 0; en += 0;
        len = en - st;
        if (len < 0) continue;

        qscore = (cnr+0) * len;
        numreads = len;

        id++;
        print chr, st, en, len, typ, id, rat, nm, KNOWN_VAR3, "TRUE", qscore, numreads;
      }

      count++;
      if (count >= nmax) break;
    }

    close(file);
  }

  # ------------------------------------------------------------
  # Main input stream
  # ------------------------------------------------------------
  NR==1 {
    for (i=1; i<=NF; i++) { $i = norm($i); idx[$i]=i; }

    reqm[1]="chrom";
    reqm[2]="start";
    reqm[3]="end";
    reqm[4]="name";
    reqm[5]="ratio";
    reqm[6]="cn";
    reqm[7]="direction";
    reqm[8]="ratio_tier";
    reqm[9]="seg_conf_z";

    for (i=1; i<=9; i++) {
      if (!(reqm[i] in idx)) {
        print "[crows_to_svcrows] ERROR: missing required column in main header: " reqm[i] > "/dev/stderr";
        exit 2;
      }
    }

    print "Chr","Start","End","Length","Type","ID","Var1","Var2","Var3","IsKnown","QScore","NumReads";

    # prepend known rows before main rows
    emit_known(KNOWN_TSV, KNOWN_N);

    next
  }

  {
    tier = lower($(idx["ratio_tier"]));
    if (!(tier=="weak" || tier=="strong")) next;

    zraw = norm($(idx["seg_conf_z"]));
    if (is_blank(zraw) || !is_num(zraw)) next;
    if ((zraw+0) < 0) next;

    chr = norm($(idx["chrom"]));
    st  = norm($(idx["start"]));
    en  = norm($(idx["end"]));
    nm  = norm($(idx["name"]));
    rat = norm($(idx["ratio"]));
    cnr = norm($(idx["cn"]));
    typ = norm($(idx["direction"]));

    if (is_blank(chr) || !is_num(st) || !is_num(en)) next;
    if (is_blank(cnr) || !is_num(cnr)) next;

    st += 0; en += 0;
    len = en - st;
    if (len < 0) next;

    id++;
    qscore = (cnr+0) * len;
    numreads = len;

    print chr, st, en, len, typ, id, rat, nm, SAMPLE, "FALSE", qscore, numreads;
  }
' "$in_tsv" > "$out_tsv"

# quick post-check
n_out="$(awk 'END{print NR}' "$out_tsv")"
if [[ "$n_out" -le 1 ]]; then
  say "WARNING: output has only a header (no rows passed main filter / known prepend)."
fi

say "Done."
