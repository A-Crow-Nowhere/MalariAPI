#!/usr/bin/env bash
# MAPI universal validator (modules & pipelines) — v2.2
# - Human-readable report -> STDERR
# - Machine JSON summary  -> STDOUT
# Exit codes: 0 = ok/warning, 1 = error

set -Eeuo pipefail

# --------- helpers ---------
err_count=0
warn_count=0
errors=() warnings=() errors_txt=() warnings_txt=()
kind="unknown"
schema_name="" schema_ver="" summary=""

json_escape() {
  local s=${1//\\/\\\\}; s=${s//\"/\\\"}; s=${s//$'\t'/\\t}; s=${s//$'\r'/}; s=${s//$'\n'/}
  printf '%s' "$s"
}
add_error()   { local c="$1" m="$2" p="${3:-}"; errors+=("{\"code\":\"$c\",\"msg\":\"$(json_escape "$m")\",\"path\":\"$(json_escape "$p")\"}"); errors_txt+=("[$c] $m  ($p)"); ((err_count++))||true; }
add_warning() { local c="$1" m="$2" p="${3:-}"; warnings+=("{\"code\":\"$c\",\"msg\":\"$(json_escape "$m")\",\"path\":\"$(json_escape "$p")\"}"); warnings_txt+=("[$c] $m  ($p)"); ((warn_count++))||true; }

print_human_summary() {
  {
    echo "── MAPI Validate ─────────────────────────────────────────────"
    echo "Kind:    $kind"
    echo "Schema:  ${schema_name}:${schema_ver}"
    if   (( err_count>0 )); then echo "Status:  ERROR    (${err_count} error(s), ${warn_count} warning(s))"
    elif (( warn_count>0)); then echo "Status:  WARNING  (${warn_count} warning(s))"
    else                        echo "Status:  OK"; fi
    [[ -n "$summary" ]] && echo "Summary: $summary"
    (( err_count>0 )) && { echo; echo "Errors:";   nl -ba <<<"$(printf '%s\n' "${errors_txt[@]}")" | sed 's/^/  /'; }
    (( warn_count>0 )) && { echo; echo "Warnings:"; nl -ba <<<"$(printf '%s\n' "${warnings_txt[@]}")" | sed 's/^/  /'; }
    echo "──────────────────────────────────────────────────────────────"
  } >&2
}
print_json_and_exit() {
  local status="ok"; (( err_count>0 )) && status="error" || { (( warn_count>0 )) && status="warning"; }
  print_human_summary
  printf '{"status":"%s","kind":"%s","schema":{"name":"%s","version":"%s"},"errors":[%s],"warnings":[%s],"summary":"%s"}\n' \
    "$status" "$kind" "$schema_name" "$schema_ver" \
    "$(IFS=,; echo "${errors[*]-}")" "$(IFS=,; echo "${warnings[*]-}")" "$(json_escape "$summary")"
  (( err_count>0 )) && exit 1 || exit 0
}
usage() {
  cat >&2 <<'EOF'
Usage:
  validate.sh -in1 <module.sh|pipeline.sh|pipeline.yml> [-in2 <metadata.yml|lock.yml>]

- Pass a module script (bin/modules/<name>.sh) or pipeline script (bin/pipelines/<name>.sh)
  Optionally pass its YAML as -in2 (bin/modules/yaml/<name>.yml or bin/pipelines/yaml/<name>.yml)
- Accepts flag variants:
    -in1 / --in1 / --in / -i      (primary input)
    -in2 / --in2 / --meta / --manifest   (secondary/metadata input)
- Positional fallback: first bare path becomes -in1
- Human report -> STDERR, single-line JSON -> STDOUT
EOF
}

# --------- arg parse (tolerant) ---------
IN1=""; IN2=""; OUTDIR=""; THREADS=""
EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in|--in1|-in1)
      if [[ $# -ge 2 ]]; then IN1="$2"; shift 2; else echo "Missing value for $1" >&2; usage; exit 2; fi
      ;;
    --in2|-in2|--meta|--manifest)
      if [[ $# -ge 2 ]]; then IN2="$2"; shift 2; else echo "Missing value for $1" >&2; usage; exit 2; fi
      ;;
    -o|--out|--outdir)
      if [[ $# -ge 2 ]]; then OUTDIR="$2"; shift 2; else echo "Missing value for $1" >&2; usage; exit 2; fi
      ;;
    -t|--threads)
      if [[ $# -ge 2 ]]; then THREADS="$2"; shift 2; else echo "Missing value for $1" >&2; usage; exit 2; fi
      ;;
    --) shift; EXTRA+=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    -*)
      EXTRA+=("$1"); shift
      if [[ $# -gt 0 && ! "$1" =~ ^- ]]; then EXTRA+=("$1"); shift; fi
      ;;
    *)
      if [[ -z "$IN1" ]]; then IN1="$1"; shift
      else EXTRA+=("$1"); shift
      fi
      ;;
  esac
done

# Guard: need a primary input
if [[ -z "${IN1}" ]]; then
  add_error "E001" "Missing required -in1/--in1 path" '$.in1'
  summary="No primary input provided"
  print_json_and_exit
fi

# --------- detect kind (script vs YAML) ---------
case "$(printf '%s' "$IN1" | tr '[:upper:]' '[:lower:]')" in
  *.sh)             kind="script" ;;
  *.yaml|*.yml)     kind="yaml" ;;
  *) add_error "E002" "Cannot infer kind from extension" '$.in1'; summary="Unknown input kind"; print_json_and_exit ;;
esac

# map to module/pipeline by location/name if possible
if [[ "$kind" == "script" ]]; then
  case "$IN1" in
    */bin/modules/*)     kind="module" ;;
    */bin/pipelines/*)   kind="pipeline" ;;
    */bin/templates/*)   add_warning "W011" "Template path treated as module" "$IN1"; kind="module" ;;
    *)                   add_warning "W010" "Script not under bin/modules or bin/pipelines; assuming module" "$IN1"; kind="module" ;;
  esac
else
  if [[ "$IN2" =~ /bin/modules/yaml/ || "$IN1" =~ /bin/modules/yaml/ ]]; then
    kind="module"
  else
    kind="pipeline"
  fi
fi

# --------- schema stub ---------
if [[ "$kind" == "module" ]]; then schema_name="mapi.module";   schema_ver="0.2"
else                              schema_name="mapi.pipeline";  schema_ver="0.2"; fi

# --------- common checks ---------
REPO_ROOT="$(cd "$(dirname "$IN1")/../.." 2>/dev/null && pwd || pwd)"
LOCK_HELPER="$REPO_ROOT/bin/_mapi_conda_lock.sh"

# Robust lock checker: two fixed-string probes (no regex footguns)
require_conda_lock() {
  local f="$1"
  local pat_rel='source "$(dirname "$0")/../_mapi_conda_lock.sh"'
  local pat_abs="source \"$REPO_ROOT/bin/_mapi_conda_lock.sh\""
  if ! grep -Fq -- "$pat_rel" "$f" && ! grep -Fq -- "$pat_abs" "$f"; then
    add_error "E105" "Script must source bin/_mapi_conda_lock.sh to enforce MAPI miniconda" "$f"
  fi
}

check_yaml_pair_exists() {
  local base="$1" type="$2" yml=""
  if [[ "$type" == "module" ]]; then
    yml="$REPO_ROOT/bin/modules/yaml/${base}.yml"
  else
    yml="$REPO_ROOT/bin/pipelines/yaml/${base}.yml"
  fi
  if [[ -f "$yml" ]]; then
    echo "$yml"
  else
    add_error "E120" "Missing metadata YAML: $yml" "$yml"
    echo ""
  fi
}

check_env_in_yaml() {
  local y="$1"
  local rel_env
  rel_env="$(awk -F: '/^env:/ {f=1} f&&/conda:/ {gsub(/[\"'\'' ]/,""); print $2; exit}' "$y" || true)"
  if [[ -z "$rel_env" ]]; then
    add_warning "W121" "No env.conda declared in YAML (will rely on caller)" "$y"
    return
  fi
  if [[ ! -f "$REPO_ROOT/$rel_env" ]]; then
    add_error "E122" "env.conda points to missing file: $rel_env" "$y:env.conda"
  fi
}

# --------- validators ---------
validate_module() {
  local script="$1"
  [[ -f "$script" ]] || add_error "E101" "Module script not a file" "$script"
  [[ -x "$script" ]] || add_warning "W102" "Module script is not executable" "$script"
  [[ -s "$script" ]] || add_error "E103" "Module script is empty" "$script"

  require_conda_lock "$script"

  local base; base="$(basename "${script%.sh}")"
  local yml; yml="$(check_yaml_pair_exists "$base" "module")"
  [[ -n "$yml" ]] && check_env_in_yaml "$yml"

  grep -Eq '^#!' "$script" || add_warning "W104" "Missing shebang on first line" "$script"

  # Try to capture help text; if we can't, don't raise E111/E112 (just W110)
  local tmp; tmp="$(mktemp)"
  if "$script" --help >"$tmp" 2>&1 || "$script" -h >"$tmp" 2>&1; then :; fi

  if ! grep -q 'Usage:' "$tmp"; then
    add_warning "W110" "No 'Usage:' in help (-h/--help)" "$script"
  else
    local usage_block
    usage_block="$(awk '/^Usage:/{f=1; print; next} f && NF==0{exit} f{print}' "$tmp" || true)"
    if ! grep -Eq '\b(--in|--in1|-i|-in1)\b|<[^>]*\.(fa|fasta|fq|fastq|bam|vcf|bed|txt|tsv)>' <<<"$usage_block"; then
      add_error "E111" "Usage must document one input file (e.g., --in, -i, or <reads.fq.gz>)" "$script:usage"
    fi
    if ! grep -Eq '\b(--out|--outdir|-o)\b|<outdir>' <<<"$usage_block"; then
      add_error "E112" "Usage must document an output directory (e.g., --out or <outdir>)" "$script:usage"
    fi
  fi

  rm -f "$tmp" || true

  if (( err_count == 0 )); then
    summary="Module valid$([[ $warn_count -gt 0 ]] && printf ' with %d warning(s)' "$warn_count")."
  else
    summary="Module has $err_count error(s) and $warn_count warning(s)."
  fi
}

validate_pipeline() {
  local script_or_yaml="$1"
  if [[ "$script_or_yaml" =~ \.ya?ml$ ]]; then
    local base; base="$(basename "${script_or_yaml%.*}")"
    local script_guess="$REPO_ROOT/bin/pipelines/${base}.sh"
    [[ -f "$script_guess" ]] || add_warning "W210" "No pipeline driver script found for YAML: $script_guess" "$script_or_yaml"
    check_env_in_yaml "$script_or_yaml"
  else
    local script="$script_or_yaml"
    [[ -f "$script" ]] || add_error "E201" "Pipeline script not a file" "$script"
    [[ -x "$script" ]] || add_warning "W202" "Pipeline script is not executable" "$script"
    require_conda_lock "$script"
    local base; base="$(basename "${script%.sh}")"
    local yml; yml="$(check_yaml_pair_exists "$base" "pipeline")"
    [[ -n "$yml" ]] && check_env_in_yaml "$yml"
  fi

  if (( err_count == 0 )); then
    summary="Pipeline valid$([[ $warn_count -gt 0 ]] && printf ' with %d warning(s)' "$warn_count")."
  else
    summary="Pipeline has $err_count error(s) and $warn_count warning(s)."
  fi
}

# --------- run ---------
if [[ "$kind" == "module" ]]; then
  validate_module "$IN1"
else
  validate_pipeline "$IN1"
fi

print_json_and_exit
