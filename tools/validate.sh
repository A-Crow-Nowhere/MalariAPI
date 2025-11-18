#!/usr/bin/env bash
set -Eeuo pipefail

# ---------- locate repo root ----------
THIS="$(readlink -f "$0" 2>/dev/null || python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$0")"
REPO_ROOT="$(cd "$(dirname "$THIS")/.." && pwd)"
MOD_DIR="$REPO_ROOT/bin/modules"
MOD_YAML_DIR="$MOD_DIR/yaml"
ENV_PREFIX_ROOT="$REPO_ROOT/envs"

# ---------- formatting ----------
RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; RST=$'\e[0m'
err(){ printf "${RED}[E%03d]${RST} %s\n" "$1" "$2"; }
warn(){ printf "${YEL}[W%03d]${RST} %s\n" "$1" "$2"; }
ok(){ printf "${GRN}[OK]${RST} %s\n" "$1"; }

# Codes:
# E100 missing YAML
# E101 name mismatch (file base vs. name:)
# E102 entry missing or not executable
# E103 env_name missing
# E104 env not installed (prefix missing)  -> suggest tools/mapi_install.sh
# E105 invalid outputs contract (tag rule or out_dir_template)
# W010 entry not under bin/modules/         (pathing)
# W020 inputs incomplete (no flags/expect ext)
# W030 outputs products missing tag/ext list
# W040 hidden files detected in modules directory
# W050 unused fields (informational)

usage(){ cat <<EOF
Usage:
  tools/validate.sh                 # validate all modules
  tools/validate.sh <module_name>   # validate just one (by YAML basename)
EOF
}

# ---------- YAML helpers (narrow subset) ----------
yaml_get_scalar(){ # key file
  local key="$1" file="$2"
  grep -E "^$key:" "$file" | head -n1 | sed "s/^$key:[[:space:]]*//"
}

yaml_has_key(){ local key="$1" file="$2"; grep -qE "^$key:" "$file"; }

yaml_get_array_items(){ # section_key file  -> prints array items for a top-level list under key:
  # Only supports flat lists like "exts: [a,b]" or
  # multi-line "- item" under "products:" handled separately below.
  local key="$1" file="$2"
  grep -E "^$key:[[:space:]]*\[" "$file" | sed -E "s/^$key:[[:space:]]*\[(.*)\].*/\1/" | tr ',' '\n' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//; s/^\"//; s/\"$//'
}

# Pull the env: block (for installer hashing parity; not needed for checks beyond presence)
has_env_block(){ awk '/^env:[[:space:]]*$/ {print "Y"; exit}' "$1" >/dev/null 2>&1; }

# Extract 'products' items (dash lists) and report per-item tag/exts
parse_products(){
  # Output lines: "TAG=<tag>|EXTS=a,b,c"
  local file="$1"
  awk '
    BEGIN{inprod=0; tag=""; exts=""}
    /^products:[[:space:]]*$/ {inprod=1; next}
    inprod==1 && /^[^[:space:]]/ {inprod=0}
    inprod==1 {
      if ($1=="-") { if (tag!=""||exts!="") {print "TAG="tag"|EXTS="exts; tag=""; exts=""} }
      if ($1=="tag:"||$1=="tag:"){ sub(/^tag:[[:space:]]*/,""); if ($0!="") tag=$0 }
      if ($1=="exts:"||$1=="exts:"){
        # exts: [a,b,c]
        match($0,/\[(.*)\]/,m); if(m[1]!=""){ exts=m[1] }
      }
    }
    END{ if (inprod==1 && (tag!=""||exts!="")) print "TAG="tag"|EXTS="exts }
  ' "$file" | sed 's/[[:space:]]//g'
}

# ---------- core checks ----------
validate_one(){
  local yml="$1"
  [[ -f "$yml" ]] || { err 100 "YAML not found: $yml"; return 1; }

  local base; base="$(basename "$yml")"
  local name_noext="${base%.*}"

  # Basic scalars
  local name entry env_name
  name="$(yaml_get_scalar "name" "$yml")"
  entry="$(yaml_get_scalar "entry" "$yml")"
  env_name="$(yaml_get_scalar "env_name" "$yml")"

  # E101: name mismatch
  if [[ -z "$name" || "$name" != "$name_noext" ]]; then
    err 101 "name: '$name' must match YAML filename '$name_noext' in $base"
  fi

  # E106: env_name must equal module name (one-to-one env per module)
  if [[ -n "$name" && -n "$env_name" && "$env_name" != "$name" ]]; then
    err 106 "env_name ('$env_name') must equal module name ('$name') for 1:1 mapping"
  fi

  # W010 / E102: entry resolution
  if [[ -z "$entry" ]]; then
    err 102 "entry missing in $base"
  else
    # must live under bin/modules/
    local entry_abs="$MOD_DIR/$entry"
    if [[ "${entry%%/*}" != "$entry" ]]; then
      # entry has a path; normalize to absolute under bin/modules
      entry_abs="$MOD_DIR/$entry"
    fi
    if [[ ! -e "$entry_abs" ]]; then
      err 102 "entry not found: $entry_abs"
    elif [[ ! -x "$entry_abs" ]]; then
      err 102 "entry not executable: $entry_abs (chmod +x)"
    else
      # Good; ensure it is actually under bin/modules (not elsewhere)
      case "$entry_abs" in
        "$MOD_DIR"/*) : ;; # ok
        *) warn 10 "entry not under bin/modules/: $entry_abs" ;;
      esac
    fi
  fi

  # E103: env_name present
  if [[ -z "$env_name" ]]; then
    err 103 "env_name missing in $base"
  else
    # E104: env prefix installed?
    local prefix="$ENV_PREFIX_ROOT/$env_name"
    [[ -d "$prefix" ]] || err 104 "env not installed: $prefix (run tools/mapi_install.sh --only $env_name)"
  fi

  # env: block present?
  has_env_block "$yml" || warn 50 "env: block missing (installer will skip building unless env exists already)"

  # Inputs sanity
  if ! yaml_has_key "inputs:" "$yml"; then
    warn 20 "no inputs: block in $base"
  else
    # quick checks: at least one flag and at least one expected ext in some item
    local flags_count exts_count
    flags_count="$(grep -nE '^[[:space:]]+flags:[[:space:]]*\[' "$yml" | wc -l | tr -d ' ')"
    exts_count="$(grep -nE '^[[:space:]]+(expects_ext|expects_exts):[[:space:]]*\[' "$yml" | wc -l | tr -d ' ')"
    [[ "$flags_count" -ge 1 ]] || warn 20 "inputs missing flags list in $base"
    [[ "$exts_count" -ge 1 ]] || warn 20 "inputs missing expects_ext list in $base"
  fi

  # Outputs contract
  if ! yaml_has_key "outputs:" "$yml"; then
    err 105 "outputs block missing in $base"
  else
    # out_dir_template must be exactly {sample_prefix}_mapi-out
    local out_dir_template
    out_dir_template="$(yaml_get_scalar "out_dir_template" "$yml")"
    if [[ "$out_dir_template" != "{sample_prefix}_mapi-out" ]]; then
      err 105 "outputs.out_dir_template must be '{sample_prefix}_mapi-out' (got: '$out_dir_template') in $base"
    fi

    # products with tag + exts
    local had_product=0 bad_product=0
    while IFS= read -r line; do
      had_product=1
      # line looks like TAG=tag|EXTS=a,b
      local TAG EXTS
      TAG="${line#TAG=}"; TAG="${TAG%%|*}"
      EXTS="${line#*|EXTS=}"; EXTS="${EXTS#,}"; EXTS="${EXTS%,}"
      if [[ -z "$TAG" || -z "$EXTS" ]]; then
        bad_product=1
      fi
    done < <(parse_products "$yml")

    if [[ $had_product -eq 0 ]]; then
      err 105 "outputs.products missing in $base"
    elif [[ $bad_product -eq 1 ]]; then
      err 105 "outputs.products entries must have non-empty tag and at least one ext in $base"
    fi
  fi

  # Hidden file sweep (warn if any hidden scripts under modules/)
  if find "$MOD_DIR" -maxdepth 1 -type f -name '.*' | grep -q .; then
    warn 40 "hidden files present under bin/modules/ (ignored by mapi)"
  fi

  ok "module '$name_noext' validated"
}

main(){
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi

  if [[ $# -eq 0 ]]; then
    # all modules
    shopt -s nullglob
    for y in "$MOD_YAML_DIR"/*.yml; do
      echo "==> $(basename "$y")"
      validate_one "$y" || true
    done
  else
    local target="$1"
    local y="$MOD_YAML_DIR/$target.yml"
    if [[ ! -f "$y" ]]; then err 100 "module YAML not found: $y"; exit 1; fi
    echo "==> $(basename "$y")"
    validate_one "$y" || true
  fi
}

main "$@"
