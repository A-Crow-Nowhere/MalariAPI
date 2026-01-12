#!/usr/bin/env bash
# Generic metadata + CLI parser for MAPI modules
# Reads the MAPI_META header in a module and builds:
#   - OPTION list (flag, short, env, default, required)
# Then parses "$@" accordingly.

# ---------------------- Internal helpers ----------------------

_mapi_meta_trim() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"  # leading
  s="${s%"${s##*[![:space:]]}"}"  # trailing
  printf '%s' "$s"
}

# OPTIONS
MAPI_META_OPT_FLAGS=()
MAPI_META_OPT_SHORTS=()
MAPI_META_OPT_ENVS=()
MAPI_META_OPT_DEFAULTS=()
MAPI_META_OPT_REQUIRED=()

MAPI_META_LOADED=0

# ------------------------- mapi_meta_load -------------------------
# Usage: mapi_meta_load "$0"
# Extracts and parses the MAPI_META header from a module script.
mapi_meta_load() {
  local script="$1"
  local header

  # Reset arrays
  MAPI_META_OPT_FLAGS=()
  MAPI_META_OPT_SHORTS=()
  MAPI_META_OPT_ENVS=()
  MAPI_META_OPT_DEFAULTS=()
  MAPI_META_OPT_REQUIRED=()

  header="$(
    awk '/^# MAPI_META_BEGIN/{flag=1;next} /^# MAPI_META_END/{flag=0} flag{print}' "$script"
  )"

  local section=""
  local line val last_idx

  while IFS= read -r line; do
    # Strip leading "#"/"# " and spaces
    line="${line#\# }"
    line="${line#\#}"
    line="$(_mapi_meta_trim "$line")"

    case "$line" in
      OPTIONS:*)  section="OPTIONS";  continue ;;
      ENV:*)      section="ENV";      continue ;;
      OUTPUTS:*)  section="OUTPUTS";  continue ;;
      RESOURCES:*)section="RESOURCES";continue ;;
      "" )        continue ;;
    esac

    # OPTIONS section
    if [[ "$section" == "OPTIONS" ]]; then
      case "$line" in
        -\ flag:*)
          # "- flag: --threads"
          val="${line#- flag:}"
          val="$(_mapi_meta_trim "$val")"
          MAPI_META_OPT_FLAGS+=("$val")
          MAPI_META_OPT_SHORTS+=("")
          MAPI_META_OPT_ENVS+=("")
          MAPI_META_OPT_DEFAULTS+=("")
          MAPI_META_OPT_REQUIRED+=("false")
          ;;
        short:*)
          val="${line#short:}"
          val="$(_mapi_meta_trim "$val")"
          last_idx=$((${#MAPI_META_OPT_SHORTS[@]} - 1))
          (( last_idx >= 0 )) && MAPI_META_OPT_SHORTS[$last_idx]="$val"
          ;;
        env:*)
          val="${line#env:}"
          val="$(_mapi_meta_trim "$val")"
          last_idx=$((${#MAPI_META_OPT_ENVS[@]} - 1))
          (( last_idx >= 0 )) && MAPI_META_OPT_ENVS[$last_idx]="$val"
          ;;
        default:*)
          val="${line#default:}"
          val="$(_mapi_meta_trim "$val")"
          last_idx=$((${#MAPI_META_OPT_DEFAULTS[@]} - 1))
          (( last_idx >= 0 )) && MAPI_META_OPT_DEFAULTS[$last_idx]="$val"
          ;;
        required:*)
          val="${line#required:}"
          val="$(_mapi_meta_trim "$val")"
          last_idx=$((${#MAPI_META_OPT_REQUIRED[@]} - 1))
          (( last_idx >= 0 )) && MAPI_META_OPT_REQUIRED[$last_idx]="$val"
          ;;
      esac
      continue
    fi

    # ENV/OUTPUTS/RESOURCES are currently ignored here.
  done <<< "$header"

  MAPI_META_LOADED=1
}

# ---------------------- mapi_meta_parse_args ----------------------
# Usage: mapi_meta_parse_args "$@"
# Populates variables named by OPTIONS.env and EXTRA[@].
# Enforces required: true for those env vars.
mapi_meta_parse_args() {
  if [[ "$MAPI_META_LOADED" -ne 1 ]]; then
    echo "[mapi_meta] ERROR: mapi_meta_load must be called before mapi_meta_parse_args" >&2
    return 1
  fi

  # Apply defaults
  local i env_name def
  for i in "${!MAPI_META_OPT_ENVS[@]}"; do
    env_name="${MAPI_META_OPT_ENVS[$i]}"
    def="${MAPI_META_OPT_DEFAULTS[$i]}"
    if [[ -n "$env_name" && -n "$def" ]]; then
      printf -v "$env_name" '%s' "$def"
    fi
  done

  EXTRA=()

  local arg matched j long short env
  while [[ $# > 0 ]]; do
    arg="$1"

    if [[ "$arg" == "--" ]]; then
      shift
      EXTRA+=("$@")
      break
    fi

    if [[ "$arg" == -* ]]; then
      matched=0
      for j in "${!MAPI_META_OPT_FLAGS[@]}"; do
        long="${MAPI_META_OPT_FLAGS[$j]}"
        short="${MAPI_META_OPT_SHORTS[$j]}"
        env="${MAPI_META_OPT_ENVS[$j]}"

        if [[ "$arg" == "$long" || ( -n "$short" && "$arg" == "$short" ) ]]; then
          matched=1
          if [[ $# -lt 2 ]]; then
            echo "[mapi_meta] ERROR: Option $arg requires a value" >&2
            return 2
          fi
          shift
          printf -v "$env" '%s' "$1"
          shift
          break
        fi
      done

      if [[ "$matched" -eq 0 ]]; then
        EXTRA+=("$arg")
        shift
      fi
    else
      EXTRA+=("$arg")
      shift
    fi
  done

  # Validate required options
  local req env_var val idx
  for idx in "${!MAPI_META_OPT_FLAGS[@]}"; do
    req="${MAPI_META_OPT_REQUIRED[$idx]}"
    env_var="${MAPI_META_OPT_ENVS[$idx]}"

    if [[ "$req" == "true" || "$req" == "yes" ]]; then
      # avoid set -u issues using ${!env_var-}
      val="${!env_var-}"
      if [[ -z "${val:-}" ]]; then
        echo "[mapi_meta] ERROR: required option for '$env_var' is not set. Check your flags." >&2
        return 3
      fi
    fi
  done

  return 0
}
