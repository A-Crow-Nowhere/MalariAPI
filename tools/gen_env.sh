#!/usr/bin/env b
# gen_env.sh — auto-generate minimal conda YAMLs from wrappers & source
# MalariAPI-standard (2025-10-29)
set -euo pipefail

# ===== Repo-aware config =====
ROOT="${ROOT:-$HOME/MalariAPI}"
YAML_OUT_DIR_DEFAULT="$ROOT/bin/modules/yaml"               # primary specs live here
PKG_YAML_SUBDIR="$ROOT/bin/packages/yaml"                               # inside bin/packages/<name>/yaml/
CHANNELS_DEFAULT=("conda-forge" "bioconda" "defaults")

# try to find a base YAML spec to subtract from (priority: envs/base.yml)
BASE_YAML_CANDIDATES=(
  "$ROOT/envs/base.yml"          # canonical runtime base (your setup)
  "$ROOT/tools/yaml/base.yml"    # fallbacks (if you keep a spec there)
  "$ROOT/tools/yaml/base.yaml"
  "$ROOT/pipeline/yaml/base.yml"
)

# ===== Tiny mapping helpers (extend these as you go) =====
declare -A PY_CONDA_MAP=(
  [numpy]="numpy" [pandas]="pandas" [scipy]="scipy" [pyyaml]="pyyaml" [ruamel.yaml]="ruamel.yaml"
  [networkx]="networkx" [matplotlib]="matplotlib" [seaborn]="seaborn" [pysam]="pysam"
  [cython]="cython" [biopython]="biopython" [pyarrow]="pyarrow" [polars]="polars"
  [tqdm]="tqdm" [requests]="requests" [jinja2]="jinja2"
)

declare -A R_CONDA_MAP=(
  [tidyverse]="r-tidyverse" [dplyr]="r-dplyr" [readr]="r-readr" [ggplot2]="r-ggplot2"
  [data.table]="r-data.table" [stringr]="r-stringr" [optparse]="r-optparse"
  [jsonlite]="r-jsonlite" [devtools]="r-devtools"
)

# lowercase keys for Bioc
declare -A BIOC_CONDA_MAP=(
  [biostrings]="bioconductor-biostrings"
  [genomicranges]="bioconductor-genomicranges"
  [iranges]="bioconductor-iranges"
  [s4vectors]="bioconductor-s4vectors"
  [summarizedexperiment]="bioconductor-summarizedexperiment"
  [delayedarray]="bioconductor-delayedarray"
  [edger]="bioconductor-edger"
  [limma]="bioconductor-limma"
)

# CLI tools to ignore (builtins/common)
IGNORE_CMDS=(
  awk bash cat cut date dirname echo env false find grep gunzip gzip head
  ln ls mkdir mv paste printf readlink realpath rm sed sort tail tar tee
  test touch tr true uniq wc xargs yes zcat jq parallel
  # extra builtins/keywords/launchers/oddities
  source exec exit break continue set shift trap ulimit umask export readonly
  local declare typeset hash times caller enable help alias unalias read getopts :
  then do done fi esac case function elif else time until select in if while for
  conda python python3
  print   # some systems have /usr/bin/print; not a dep we care about
)



# ===== CLI =====
usage() {
  cat <<EOF
Usage: $(basename "$0") --name <env-name> --kind <module|pipeline|package> [options]

Required:
  --name NAME               Logical name (module/pipeline/package name)
  --kind KIND               One of: module | pipeline | package

What to scan (repeatable):
  --file PATH               Add a file to scan (.sh, .py, .R, .r)
  --dir  DIR                Recursively scan a directory

Output control:
  --yaml-out DIR            Where to write YAML (default: $YAML_OUT_DIR_DEFAULT)
  --channels CSV            e.g. "conda-forge,bioconda,defaults"
  --no-base-subtract        Do NOT subtract base env dependencies
  --also-copy-package       If --kind=package, copy to bin/packages/NAME/$PKG_YAML_SUBDIR/NAME.yml

Force-add dependencies (all are version-aware):
  --add CMD[=VER]           Add a CLI tool (e.g., samtools=1.20)
  --python PKG[=VER]        Add a Python package
  --r CRANPKG[=VER]         Add a CRAN package (maps to r-<pkg>[=ver])
  --bioc PKG[=VER]          Add a Bioconductor package (bioconductor-<pkg>[=ver])
  --conda NAME[=VER]        Add an exact conda dep (raw passthrough)

Examples:
  $0 --name fastp --kind module --file "$ROOT/bin/modules/fastp.sh" --add fastp
  $0 --name bwa_align --kind module --dir "$ROOT/bin/modules/bwa_align" --conda "bwa=0.7.17"
  $0 --name circrowseq --kind pipeline --dir "$ROOT/bin/pipeline" --python "pandas=2.2.2"
  $0 --name mypkg --kind package --dir "$ROOT/bin/packages/mypkg" --bioc GenomicRanges=1.54.1 --also-copy-package
EOF
  exit 2
}

NAME="" KIND="" ;
FILES=() DIRS=()
EXTRA_CMDS=() EXTRA_PY=() EXTRA_R=() EXTRA_BIOC=() EXTRA_CONDA=()
YAML_OUT_DIR="$YAML_OUT_DIR_DEFAULT"
CHANNELS=("${CHANNELS_DEFAULT[@]}")
SUBTRACT_BASE=1
COPY_IN_PACKAGE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --name) NAME="$2"; shift 2;;
    --kind) KIND="$2"; shift 2;;
    --file) FILES+=("$2"); shift 2;;
    --dir)  DIRS+=("$2"); shift 2;;
    --yaml-out) YAML_OUT_DIR="$2"; shift 2;;
    --channels) IFS=',' read -r -a CHANNELS <<< "$2"; shift 2;;
    --no-base-subtract) SUBTRACT_BASE=0; shift;;
    --also-copy-package) COPY_IN_PACKAGE=1; shift;;
    --add) EXTRA_CMDS+=("$2"); shift 2;;
    --python) EXTRA_PY+=("$2"); shift 2;;
    --r) EXTRA_R+=("$2"); shift 2;;
    --bioc) EXTRA_BIOC+=("$2"); shift 2;;
    --conda) EXTRA_CONDA+=("$2"); shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -n "$NAME" && -n "$KIND" ]] || { echo "ERROR: --name and --kind are required"; usage; }
case "$KIND" in module|pipeline|package) :;; *) echo "ERROR: --kind must be module|pipeline|package"; usage;; esac

# ===== Gather files =====
for D in "${DIRS[@]}"; do
  while IFS= read -r -d '' f; do FILES+=("$f"); done < <(find "$D" -type f \( -name "*.sh" -o -name "*.py" -o -name "*.R" -o -name "*.r" \) -print0)
done
if [[ ${#FILES[@]} -eq 0 ]]; then
  case "$KIND" in
    module)
      [[ -f "$ROOT/bin/modules/$NAME.sh" ]] && FILES+=("$ROOT/bin/modules/$NAME.sh")
      [[ -d "$ROOT/bin/modules/$NAME" ]] && DIRS+=("$ROOT/bin/modules/$NAME")
      ;;
    pipeline)
      [[ -f "$ROOT/bin/pipeline.sh" ]] && FILES+=("$ROOT/bin/pipeline.sh")
      [[ -d "$ROOT/bin/pipeline" ]] && DIRS+=("$ROOT/bin/pipeline")
      ;;
    package)
      [[ -d "$ROOT/bin/packages/$NAME" ]] && DIRS+=("$ROOT/bin/packages/$NAME")
      ;;
  esac
  for D in "${DIRS[@]}"; do
    while IFS= read -r -d '' f; do FILES+=("$f"); done < <(find "$D" -type f \( -name "*.sh" -o -name "*.py" -o -name "*.R" -o -name "*.r" \) -print0)
  done
fi
[[ ${#FILES[@]} -gt 0 ]] || { echo "ERROR: No files to scan. Use --file/--dir or ensure defaults exist."; exit 1; }

# ===== Utils =====
normlist() { tr '[:upper:]' '[:lower:]' | tr -s ' ' '\n' | sed '/^$/d' | sort -u; }
join_by() { local IFS="$1"; shift; echo "$*"; }
contains() { local item="$1"; shift; for x in "$@"; do [[ "$x" == "$item" ]] && return 0; done; return 1; }

# Parse name=version (or plain name). Echo "<name>" and "<=version or empty>" via two vars.
split_nv() {
  local s="$1"
  if [[ "$s" == *"="* ]]; then
    NV_NAME="${s%%=*}"
    NV_VER="${s#*=}"
  else
    NV_NAME="$s"
    NV_VER=""
  fi
}

# Version-aware mapper: given input like "dplyr=1.1.4" and a mapped base "r-dplyr",
# it returns "r-dplyr=1.1.4". If no version: "r-dplyr".
apply_version_to_mapped() {
  local mapped="$1" ver="$2"
  if [[ -n "$ver" ]]; then echo "${mapped}=${ver}"; else echo "$mapped"; fi
}

# Bash-native: strip heredoc bodies from a shell script (prints modified file to stdout)
# Handles: <<EOF, <<'EOF', <<-"EOF" (delimiter is captured as bare word)
# Strip heredoc bodies from a shell script (prints modified file to stdout)
# Handles: <<DELIM, <<'DELIM', <<-"DELIM"
sh_strip_heredocs() {
  local f="$1"
  local in=0 end="" line rest delim first

  while IFS= read -r line || [[ -n "$line" ]]; do
    if (( in )); then
      if [[ "$line" == "$end" ]]; then
        in=0
      fi
      # NOTE: we do not print heredoc body lines
      continue
    fi

    if [[ "$line" == *"<<"* ]]; then
      # extract text after first << (supports <<- too)
      rest="${line#*<<}"
      # optional -
      if [[ "$rest" == -* ]]; then rest="${rest#-}"; fi
      # optional starting quote
      if [[ "$rest" == \"* || "$rest" == \'* ]]; then
        first="${rest:0:1}"
        rest="${rest:1}"
        delim="${rest%%"$first"*}"
      else
        delim="${rest%%[[:space:]]*}"
      fi
      if [[ -n "$delim" ]]; then
        end="$delim"; in=1
      fi
      echo "$line"
      continue
    fi

    echo "$line"
  done < "$f"
}

# Emit ALL heredoc bodies (concatenated) from a shell script to stdout.
# We don't assume what invoked the heredoc; we inspect bodies later.
sh_collect_all_heredoc_bodies() {
  local f="$1"
  local in=0 end="" line rest delim first
  while IFS= read -r line || [[ -n "$line" ]]; do
    if (( in )); then
      if [[ "$line" == "$end" ]]; then
        in=0
        continue
      fi
      echo "$line"
      continue
    fi
    if [[ "$line" == *"<<"* ]]; then
      rest="${line#*<<}"
      if [[ "$rest" == -* ]]; then rest="${rest#-}"; fi
      if [[ "$rest" == \"* || "$rest" == \'* ]]; then
        first="${rest:0:1}"
        rest="${rest:1}"
        delim="${rest%%"$first"*}"
      else
        delim="${rest%%[[:space:]]*}"
      fi
      if [[ -n "$delim" ]]; then
        end="$delim"; in=1
        continue
      fi
    fi
  done < "$f"
}

# Emit only the bodies of heredocs that start from a python launcher line
# Matches lines that mention "python" (python/python3/conda run ... python) AND start a heredoc.
sh_collect_python_imports_from_heredocs() {
  local f="$1"
  local in=0 end="" want=0 line rest delim first

  while IFS= read -r line || [[ -n "$line" ]]; do
    if (( in )); then
      if [[ "$line" == "$end" ]]; then
        in=0; want=0
        continue
      fi
      (( want )) && echo "$line"
      continue
    fi

    # Must contain "python" and "<<"
    if [[ "$line" == *python* && "$line" == *"<<"* ]]; then
      rest="${line#*<<}"
      # optional -
      if [[ "$rest" == -* ]]; then rest="${rest#-}"; fi
      # extract delimiter (quoted or bare)
      if [[ "$rest" == \"* || "$rest" == \'* ]]; then
        first="${rest:0:1}"
        rest="${rest:1}"
        delim="${rest%%"$first"*}"
      else
        delim="${rest%%[[:space:]]*}"
      fi
      if [[ -n "$delim" ]]; then
        end="$delim"; in=1; want=1
        continue
      fi
    fi
  done < "$f"
}

# Find .py files referenced on the command lines of a shell script (after heredocs stripped).
sh_find_referenced_pyfiles() {
  local f="$1"
  sh_strip_heredocs "$f" | awk '
    {
      for (i=1;i<=NF;i++) {
        if ($i ~ /\.py$/) print $i
      }
    }
  ' | sed 's/[;"'\'']//g' | sort -u
}

# ===== Extractors =====
PY_PKGS=() R_PKGS=() BIOC_PKGS=() CMDS=()

declare -A BASH_FUNS
collect_bash_functions() {
  local file="$1"
  grep -E '^[a-zA-Z_][a-zA-Z0-9_]*\(\)[[:space:]]*\{' -n "$file" 2>/dev/null | \
  while IFS= read -r line; do
    local fname="${line%%(*}"; fname="${fname##*:}"
    [[ -n "$fname" ]] && BASH_FUNS["$fname"]=1
  done || true
}

for f in "${FILES[@]}"; do
  [[ "$f" == *.sh ]] && collect_bash_functions "$f"
done

extract_python_imports() { sed -E 's/#.*$//' "$1" | grep -E '^\s*(from|import)\s+' || true; }
extract_r_libs()       { sed -E 's/#.*$//' "$1" | grep -E 'library\(|require\(' || true; }

extract_shell_cmds() {
  local f="$1"
  sh_strip_heredocs "$f" | \
  awk '
    /^[[:space:]]*#/ {next}
    {
      n=split($0, segs, /(\&\&|;|\|\||\|)/)
      for (i=1; i<=n; i++) {
        s=segs[i]
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
        if (s ~ /^[[:alnum:]_.+-]+([[:space:]]|$)/) {
          cmd=s
          sub(/[[:space:]].*$/, "", cmd)
          print cmd
        }
      }
    }
  ' 2>/dev/null | \
  grep -Ev '^(then|do|done|fi|esac|case|function|elif|else|time|until|select|in|if|while|for)$' || true
}



for f in "${FILES[@]}"; do
  case "$f" in
    *.py)
      while read -r line; do
        [[ -z "$line" ]] && continue
        if [[ "$line" =~ ^[[:space:]]*from[[:space:]]+([a-zA-Z0-9_\.]+)[[:space:]]+import ]]; then
          mod="${BASH_REMATCH[1]%%.*}"
          PY_PKGS+=("$mod")
        elif [[ "$line" =~ ^[[:space:]]*import[[:space:]]+(.+) ]]; then
          mods="${BASH_REMATCH[1]}"; IFS=',' read -r -a arr <<< "$mods"
          for m in "${arr[@]}"; do m="${m%% *}"; m="${m%%.*}"; [[ -n "$m" ]] && PY_PKGS+=("$m"); done
        fi
      done < <(extract_python_imports "$f")
      ;;
    *.R|*.r)
      while read -r line; do
        [[ -z "$line" ]] && continue
        # Match: library(pkg), library("pkg"), require('pkg'), etc.
        if [[ "$line" =~ (library|require)\([\"']?([A-Za-z0-9_.]+)[\"']?\) ]]; then
          pkg="${BASH_REMATCH[2]}"
          # Collect raw; we’ll classify as Bioc or CRAN later
          R_PKGS+=("$pkg")
        fi
      done < <(extract_r_libs "$f")
      ;;
    *.sh)
      # 1) Track functions so we don't count them as external tools
      collect_bash_functions "$f"

      # 2) Extract real CLI executables only (exclude builtins/keywords/functions/aliases)
      while read -r c; do
        [[ -n "${c:-}" ]] || continue
        # ignore common builtins/keywords and stuff we never want as deps
        if contains "$c" "${IGNORE_CMDS[@]}"; then continue; fi
        case "$(type -t -- "$c" 2>/dev/null || true)" in
          file)  ;;     # keep only real executables found in PATH
          *)     continue;;
        esac
        # ignore conda/python launchers as "deps"
        [[ "$c" == "conda" || "$c" == "python" || "$c" == "python3" ]] && continue
        # ignore functions defined in this script as a final guard
        [[ -n "${BASH_FUNS[$c]:-}" ]] && continue
        CMDS+=("$c")
      done < <(extract_shell_cmds "$f")

      # 3) Detect Python imports from ANY heredoc body (even if invoked via ${RUN_PY[@]})
      if all_blocks="$(sh_collect_all_heredoc_bodies "$f")"; then
        # Heuristic: treat a heredoc as Python if it has a python shebang
        # or contains python-style import statements.
        # We parse imports regardless and just collect matches.
        while read -r line; do
          [[ -z "$line" ]] && continue
          if [[ "$line" =~ ^[[:space:]]*from[[:space:]]+([A-Za-z0-9_\.]+)[[:space:]]+import ]]; then
            mod="${BASH_REMATCH[1]%%.*}"
            PY_PKGS+=("$mod")
          elif [[ "$line" =~ ^[[:space:]]*import[[:space:]]+(.+) ]]; then
            mods="${BASH_REMATCH[1]}"
            IFS=',' read -r -a arr <<< "$mods"
            for m in "${arr[@]}"; do
              m="${m%% *}"; m="${m%%.*}"
              [[ -n "$m" ]] && PY_PKGS+=("$m")
            done
          fi
        done <<< "$all_blocks"
      fi


      # 4) Also scan any referenced .py files invoked by the shell script
      while read -r pyref; do
        [[ -z "${pyref:-}" ]] && continue
        # strip simple quotes/trailing punctuation
        pyref="${pyref%\"}"; pyref="${pyref%\'}"; pyref="${pyref%;}"
        # resolve relative path against the shell file's directory
        if [[ "$pyref" != /* ]]; then
          pyref="$(cd "$(dirname "$f")" && printf '%s' "$(pwd)/$pyref")"
        fi
        [[ -f "$pyref" ]] || continue
        while read -r line; do
          [[ -z "$line" ]] && continue
          if [[ "$line" =~ ^[[:space:]]*from[[:space:]]+([A-Za-z0-9_\.]+)[[:space:]]+import ]]; then
            mod="${BASH_REMATCH[1]%%.*}"
            PY_PKGS+=("$mod")
          elif [[ "$line" =~ ^[[:space:]]*import[[:space:]]+(.+) ]]; then
            mods="${BASH_REMATCH[1]}"
            IFS=',' read -r -a arr <<< "$mods"
            for m in "${arr[@]}"; do
              m="${m%% *}"; m="${m%%.*}"
              [[ -n "$m" ]] && PY_PKGS+=("$m")
            done
          fi
        done < <(sed -E 's/#.*$//' "$pyref" | grep -E '^\s*(from|import)\s+' || true)
      done < <(sh_find_referenced_pyfiles "$f")
      ;;
  esac
done

# ===== Merge user extras (version-aware) =====
# CLI tools (as-is; may include =ver)
CMDS+=("${EXTRA_CMDS[@]}")

# Python (keep version if provided)
PY_EXTRAS=()
for e in "${EXTRA_PY[@]:-}"; do PY_EXTRAS+=("$e"); done
# R CRAN (raw values like dplyr or dplyr=1.1.4 — mapping happens later)
R_EXTRAS=()
for e in "${EXTRA_R[@]:-}"; do R_EXTRAS+=("$e"); done
# Bioc (raw; mapping later)
BIOC_EXTRAS=()
for e in "${EXTRA_BIOC[@]:-}"; do BIOC_EXTRAS+=("$e"); done
# Raw conda deps
CONDA_RAW=("${EXTRA_CONDA[@]:-}")

# Normalize & unique (for detected lists)
mapfile -t CMDS       < <(printf "%s\n" "${CMDS[@]:-}" | normlist)
mapfile -t PY_PKGS    < <(printf "%s\n" "${PY_PKGS[@]:-}" | sed 's/[[:space:]]//g' | sed '/^$/d' | sort -u)

# R/BIoc detected need classification; we keep raw for extras so users can pass versions
# Split detected R into CRAN vs Bioc by name map
R_DET=() BIOC_DET=()
for r in "${R_PKGS[@]:-}"; do
  [[ -z "$r" ]] && continue
  low="$(echo "$r" | tr '[:upper:]' '[:lower:]')"
  if [[ -n "$low" && -n "${BIOC_CONDA_MAP[$low]+x}" ]]; then
    BIOC_DET+=("$r")
  else
    R_DET+=("$r")
  fi
done

# Compose final pools (detected first, then extras, keeping potential duplicates for now)
R_ALL=( "${R_DET[@]:-}" "${R_EXTRAS[@]:-}" )
BIOC_ALL=( "${BIOC_DET[@]:-}" "${BIOC_EXTRAS[@]:-}" )

# ===== Map to conda deps (version-aware for extras) =====
CONDA_DEPS=()
PIP_DEPS=()  # kept for future use

# ==== Build CONDA_DEPS safely (all guards for set -u) =========================

# 1) CLI tools (as-is, allow optional =version pins)
for c in "${CMDS[@]:-}"; do
  [[ -n "${c:-}" ]] || continue
  CONDA_DEPS+=("$c")
done

# 2) Python skip set / helper
PY_SKIP_SET="sys os pathlib argparse typing importlib pkgutil logging json re subprocess shutil glob math random time datetime itertools collections"
py_skip_has() { [[ " $PY_SKIP_SET " == *" ${1:-} "* ]]; }

# 3) Python (detected) — guarded and map-checked
for p in "${PY_PKGS[@]:-}"; do
  [[ -n "${p:-}" ]] || continue
  [[ "$p" =~ ^_ ]] && continue
  [[ "$p" == "__future__" ]] && continue
  low="$(echo "$p" | tr '[:upper:]' '[:lower:]')"
  py_skip_has "$low" && continue

  if [[ -n "${PY_CONDA_MAP[$p]+x}" ]]; then
    mapped="${PY_CONDA_MAP[$p]}"
  elif [[ -n "${PY_CONDA_MAP[$low]+x}" ]]; then
    mapped="${PY_CONDA_MAP[$low]}"
  else
    mapped="$low"
  fi
  [[ -n "${mapped:-}" ]] && CONDA_DEPS+=("$mapped")
done

# 4) Python extras (may include versions) — guarded map lookup
for pe in "${PY_EXTRAS[@]:-}"; do
  [[ -n "${pe:-}" ]] || continue
  split_nv "$pe"  # sets NV_NAME NV_VER
  [[ -n "${NV_NAME:-}" ]] || continue
  low="$(echo "$NV_NAME" | tr '[:upper:]' '[:lower:]')"

  if [[ -n "${PY_CONDA_MAP[$NV_NAME]+x}" ]]; then
    base="${PY_CONDA_MAP[$NV_NAME]}"
  elif [[ -n "${PY_CONDA_MAP[$low]+x}" ]]; then
    base="${PY_CONDA_MAP[$low]}"
  else
    base="$low"
  fi
  [[ -n "${base:-}" ]] && CONDA_DEPS+=("$(apply_version_to_mapped "$base" "${NV_VER:-}")")
done

# 5) R (CRAN; detected+extras merged into R_ALL earlier) — guarded map lookup
for r in "${R_ALL[@]:-}"; do
  [[ -n "${r:-}" ]] || continue
  split_nv "$r"; [[ -n "${NV_NAME:-}" ]] || continue
  low="$(echo "$NV_NAME" | tr '[:upper:]' '[:lower:]')"

  if [[ -n "${R_CONDA_MAP[$NV_NAME]+x}" ]]; then
    base="${R_CONDA_MAP[$NV_NAME]}"
  elif [[ -n "${R_CONDA_MAP[$low]+x}" ]]; then
    base="${R_CONDA_MAP[$low]}"
  else
    base="r-$low"
  fi
  [[ -n "${base:-}" ]] && CONDA_DEPS+=("$(apply_version_to_mapped "$base" "${NV_VER:-}")")
done

# 6) Bioconductor (detected+extras in BIOC_ALL) — guarded map lookup
for b in "${BIOC_ALL[@]:-}"; do
  [[ -n "${b:-}" ]] || continue
  split_nv "$b"; [[ -n "${NV_NAME:-}" ]] || continue
  low="$(echo "$NV_NAME" | tr '[:upper:]' '[:lower:]')"

  if [[ -n "${BIOC_CONDA_MAP[$low]+x}" ]]; then
    base="${BIOC_CONDA_MAP[$low]}"
  else
    base="bioconductor-$low"
  fi
  [[ -n "${base:-}" ]] && CONDA_DEPS+=("$(apply_version_to_mapped "$base" "${NV_VER:-}")")
done

# 7) Raw conda deps (passthrough, e.g., "bwa=0.7.17")
for raw in "${CONDA_RAW[@]:-}"; do
  [[ -n "${raw:-}" ]] || continue
  CONDA_DEPS+=("$raw")
done

# 8) Unique conda deps (case-insensitive by name; keep first occurrence to preserve pins)
declare -A SEEN
UNIQED=()
for d in "${CONDA_DEPS[@]:-}"; do
  [[ -n "${d:-}" ]] || continue
  base="${d%%=*}"
  key="$(echo "$base" | tr '[:upper:]' '[:lower:]')"
  [[ -n "${key:-}" ]] || continue
  if [[ -z "${SEEN[$key]+x}" ]]; then
    SEEN[$key]=1
    UNIQED+=("$d")
  fi
done
CONDA_DEPS=("${UNIQED[@]:-}")

# ===== Subtract base env if requested (robust under set -u) ===================
BASE_DEPS=()  # ensure declared
if (( ${SUBTRACT_BASE:-0} )); then
  tmpfile="$(mktemp)"
  for cand in "${BASE_YAML_CANDIDATES[@]}"; do
    [[ -n "$cand" && -f "$cand" ]] || continue
    awk '
      BEGIN{in_dep=0}
      /^[[:space:]]*dependencies:[[:space:]]*$/ {in_dep=1; next}
      /^[[:alpha:]]/ && !/^[[:space:]]/ {in_dep=0}
      in_dep==1 {
        if ($0 ~ /^[-][[:space:]]+pip$/) next
        if ($0 ~ /^[-][[:space:]]+/) {
          gsub(/^[-][[:space:]]+/, "", $0)
          gsub(/=.*/, "", $0)
          print tolower($0)
        }
      }
    ' "$cand" | sed '/^$/d' | sort -u > "$tmpfile"
    mapfile -t BASE_DEPS < "$tmpfile"
    rm -f "$tmpfile"
    break
  done
fi

if (( ${#BASE_DEPS[@]} > 0 )); then
  FILTERED=()
  for d in "${CONDA_DEPS[@]:-}"; do
    [[ -n "${d:-}" ]] || continue
    lowbase="$(echo "${d%%=*}" | tr '[:upper:]' '[:lower:]')"
    if contains "$lowbase" "${BASE_DEPS[@]}"; then
      : # present in base -> drop
    else
      FILTERED+=("$d")
    fi
  done
  CONDA_DEPS=("${FILTERED[@]:-}")
fi

# ===== Produce YAML =====
mkdir -p "$YAML_OUT_DIR"
YAML_PATH="$YAML_OUT_DIR/$NAME.yml"

{
  echo "name: $NAME"
  echo "channels:"
  for ch in "${CHANNELS[@]}"; do echo "  - $ch"; done
  echo "dependencies:"
  # include r-base if any R/Bioc deps exist
  have_r=0
  for d in "${CONDA_DEPS[@]}"; do
    [[ "$d" == r-* || "$d" == bioconductor-* ]] && have_r=1
  done
  [[ $have_r -eq 1 ]] && echo "  - r-base"
  for d in "${CONDA_DEPS[@]}"; do echo "  - $d"; done
  if [[ ${#PIP_DEPS[@]} -gt 0 ]]; then
    echo "  - pip"
    echo "  - pip:"
    for p in "${PIP_DEPS[@]}"; do echo "      - $p"; done
  fi
} > "$YAML_PATH"

echo "[gen_env] Wrote: $YAML_PATH"

# Optionally copy into package subdir for publication
if [[ "$KIND" == "package" && $COPY_IN_PACKAGE -eq 1 ]]; then
  DEST="$ROOT/bin/packages/$NAME/$PKG_YAML_SUBDIR"
  mkdir -p "$DEST"
  cp -f "$YAML_PATH" "$DEST/$NAME.yml"
  echo "[gen_env] Copied package YAML to: $DEST/$NAME.yml"
fi

# ===== Summary (show cleaned sets) ============================================
# Build a clean, de-duplicated view of detected Python imports (non-stdlib)
PY_PKGS_DET_VIEW=()
for p in "${PY_PKGS[@]:-}"; do
  [[ -n "${p:-}" ]] || continue
  [[ "$p" =~ ^_ ]] && continue
  [[ "$p" == "__future__" ]] && continue
  low="$(echo "$p" | tr '[:upper:]' '[:lower:]')"
  # same stdlib skip list used during mapping
  [[ " $PY_SKIP_SET " == *" $low "* ]] && continue
  PY_PKGS_DET_VIEW+=("$low")
done
mapfile -t PY_PKGS_DET_VIEW < <(printf "%s\n" "${PY_PKGS_DET_VIEW[@]:-}" | sort -u)

# Unique views (already filtered earlier)
mapfile -t CMDS_VIEW < <(printf "%s\n" "${CMDS[@]:-}" | sort -u)

join_csv() { local IFS=,; echo "$*"; }

echo "[gen_env] Summary for $NAME ($KIND):"
echo "  Scanned files: ${#FILES[@]}"
echo "  CLI tools:    $( ((${#CMDS_VIEW[@]})) && join_csv "${CMDS_VIEW[@]}" )"
echo "  Python pkgs:  $( ((${#PY_PKGS_DET_VIEW[@]})) && join_csv "${PY_PKGS_DET_VIEW[@]}" )"
echo "  R pkgs:       $( ((${#R_DET[@]})) && join_csv "${R_DET[@]}" )"
echo "  Bioc pkgs:    $( ((${#BIOC_DET[@]})) && join_csv "${BIOC_DET[@]}" )"
if (( ${#BASE_DEPS[@]} > 0 )); then
  echo "  Base subtract: ${#BASE_DEPS[@]} deps removed"
else
  echo "  Base subtract: none"
fi
