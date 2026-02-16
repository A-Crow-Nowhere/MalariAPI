#!/usr/bin/env bash
set -euo pipefail

# BED4 (chrom start end label) -> IGV BED9(+Program) with itemRgb
# Gain -> 76,185,68 ; Loss -> 240,101,67
# Usage:
#   ./bed4_gainloss_to_igv.sh [-t gain|loss|auto] [-T "Track"] [-D "Desc"] [-p "Program"] [in.bed|in.bed.gz] > out.bed
# Examples:
#   ./bed4_gainloss_to_igv.sh -t auto jasmine.bed > jasmine.igv.bed
#   zcat in.bed.gz | ./bed4_gainloss_to_igv.sh -t auto > out.bed

TRACK_NAME="SVCROWS Default"
DESCRIPTION="CNV Events"
PROGRAM="SVCROWS"
MODE="auto"                 # gain|loss|auto
RGB_GAIN="76,185,68"
RGB_LOSS="240,101,67"
INPUT_FILE=""

usage() {
  cat <<'EOF' >&2
BED4 (chrom start end label) -> IGV BED9 with itemRgb

Options:
  -t MODE          gain|loss|auto  (default: auto)
  -T "NAME"        Track name (default: SVCROWS Default)
  -D "DESC"        Track description (default: CNV Events)
  -p "PROGRAM"     Final column tag (default: SVCROWS)
  -h|--help        Show help

If no file is given, reads stdin; errors out if nothing is provided.
EOF
  exit 1
}

# parse args
while (( "$#" )); do
  case "$1" in
    -t) MODE="${2:-}"; shift 2;;
    -T) TRACK_NAME="${2:-}"; shift 2;;
    -D) DESCRIPTION="${2:-}"; shift 2;;
    -p) PROGRAM="${2:-}"; shift 2;;
    -h|--help) usage;;
    --) shift; break;;
    -*)
      echo "Unknown option: $1" >&2; usage;;
    *)
      if [[ -z "$INPUT_FILE" ]]; then INPUT_FILE="$1"; else echo "Too many inputs" >&2; usage; fi
      shift;;
  esac
done

# ensure there will be input (file or pipe)
if [[ -z "$INPUT_FILE" && -t 0 ]]; then
  echo "No input file and no data on stdin." >&2
  usage
fi
if [[ -n "$INPUT_FILE" && ! -e "$INPUT_FILE" ]]; then
  echo "Input not found: $INPUT_FILE" >&2
  exit 2
fi

# header
printf 'track name="%s" description="%s" itemRgb=on\n' "$TRACK_NAME" "$DESCRIPTION"
printf '#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tProgram\n'

# reader
reader() {
  if [[ -z "$INPUT_FILE" || "$INPUT_FILE" == "-" ]]; then
    cat
  else
    case "$INPUT_FILE" in
      *.gz) zcat -- "$INPUT_FILE" ;;
      *)    cat -- "$INPUT_FILE" ;;
    esac
  fi
}

# transform
reader | awk -v MODE="$MODE" -v RG_GAIN="$RGB_GAIN" -v RG_LOSS="$RGB_LOSS" -v PROG="$PROGRAM" '
BEGIN{
  FS = "[ \t]+"; OFS = "\t"; IGNORECASE = 1
}
# classify label -> gain/loss/unknown
function classify(lbl,   t){
  t = lbl
  sub(/\r$/,"",t)                             # strip CR from CRLF
  gsub(/[[:space:]]/,"",t)
  if (t ~ /^(loss|del|deletion|down|cnv:loss|copy_number_loss)$/) return "loss"
  if (t ~ /^(gain|dup|dup:tandem|ampl|amp|increase|up|cnv:gain|copy_number_gain)$/) return "gain"
  return "unknown"
}

# skip empties and comment/track lines
/^[[:space:]]*$/ { next }
$0 ~ /^[[:space:]]*#/ { next }
$0 ~ /^[[:space:]]*track[[:space:]]/ { next }

{
  # Expect at least 4 fields; take first 4 literally (no guessing)
  if (NF < 4) next
  chrom = $1; start = $2; end = $3; label = $4
  sub(/\r$/,"",label)

  # simple sanity for coords
  if (start !~ /^[0-9]+$/ || end !~ /^[0-9]+$/) next
  if ((start+0) >= (end+0)) next

  name_out = label
  rgb = RG_GAIN
  if (MODE == "auto") {
    cls = classify(label)
    if (cls=="loss") { rgb = RG_LOSS; name_out = "Loss" }
    else if (cls=="gain") { rgb = RG_GAIN; name_out = "Gain" }
    else { rgb = RG_GAIN; name_out = label }
  } else if (MODE == "loss") {
    rgb = RG_LOSS
  } else { # MODE == "gain"
    rgb = RG_GAIN
  }

  print chrom, start, end, name_out, 1, ".", start, end, rgb, PROG
}
'
