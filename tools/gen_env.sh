#!/usr/bin/env bash
# Generate a minimal conda env YAML for a MAPI module wrapper.
# Scans: inline Python (heredoc PYCODE), R libraries, and common CLI tools.
# Writes: ~/MalariAPI/env/yaml/<module>.yml by default.
# Requires: standard POSIX tools (awk/sed/grep/sort/comm); mamba (optional for --test).
set -euo pipefail

# ---------- Defaults (MAPI-rooted) ----------
MAPI_ROOT="${MAPI_ROOT:-$HOME/MalariAPI}"
DEFAULT_OUT_DIR="$MAPI_ROOT/env/yaml"
DEFAULT_BASE_CONDA="$MAPI_ROOT/env/base/conda.txt"
DEFAULT_BASE_PIP="$MAPI_ROOT/env/base/pip.txt"
CHANNELS="conda-forge,bioconda,defaults"
ENV_PY_VERSION=">=3.10"   # adjust if you prefer pinning exact versions
ASSUME_PIP=()             # force these to pip even if conda exists
FORCE_CONDA=()            # force these to conda even if we might put in pip
DRY_RUN=0
OUT_YML=""
ENV_NAME=""               # default: <module>-env
QUIET=0

usage() {
  cat <<EOF
Usage: $(basename "$0") --from <wrapper.sh> [options]

Required:
  --from <path>           Wrapper to scan (e.g., $MAPI_ROOT/bin/modules/<name>.sh)

Outputs:
  --out-yml <path>        Output YAML file (default: $DEFAULT_OUT_DIR/<module>.yml)
  --env-name <name>       Conda env name in YAML (default: <module>-env)
  --channels <csv>        Channels (default: $CHANNELS)

Base subtraction (avoid double-installs):
  --base-conda <file>     List of conda packages already in base (default: $DEFAULT_BASE_CONDA)
  --base-pip   <file>     List of pip packages already in base (default: $DEFAULT_BASE_PIP)

Control:
  --assume-pip <pkg>      Force listed pkg(s) to pip (repeatable)
  --force-conda <pkg>     Force listed pkg(s) to conda (repeatable)
  --py-version <spec>     Python version spec (default: $ENV_PY_VERSION)
  --test                  Run: mamba create --dry-run -n <env> -f <yaml>
  -q|--quiet              Less chatter
  -h|--help               This help

Notes:
  * Full MAPI-rooted paths are always used in output.
  * We map Python imports & R libraries to conda names via curated rules + heuristics.
  * Unmappable Python imports go to pip: section (unless forced to conda).
EOF
}

log() { [[ $QUIET -eq 1 ]] || printf '%s\n' "$*" >&2; }

# ---------- Parse args ----------
WRAPPER=""
BASE_CONDA="$DEFAULT_BASE_CONDA"
BASE_PIP="$DEFAULT_BASE_PIP"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from)           WRAPPER="$2"; shift 2;;
    --out-yml)        OUT_YML="$2"; shift 2;;
    --env-name)       ENV_NAME="$2"; shift 2;;
    --channels)       CHANNELS="$2"; shift 2;;
    --base-conda)     BASE_CONDA="$2"; shift 2;;
    --base-pip)       BASE_PIP="$2"; shift 2;;
    --assume-pip)     ASSUME_PIP+=("$2"); shift 2;;
    --force-conda)    FORCE_CONDA+=("$2"); shift 2;;
    --py-version)     ENV_PY_VERSION="$2"; shift 2;;
    --test)           DRY_RUN=1; shift;;
    -q|--quiet)       QUIET=1; shift;;
    -h|--help)        usage; exit 0;;
    *) log "Unknown arg: $1"; usage; exit 2;;
  esac
done

[[ -n "$WRAPPER" ]] || { log "[ERROR] --from <wrapper.sh> is required"; usage; exit 2; }
[[ -f "$WRAPPER" ]] || { log "[ERROR] wrapper not found: $WRAPPER"; exit 2; }

# ---------- Derive names / paths ----------
MODULE="$(basename "${WRAPPER%.sh}")"
[[ -n "$ENV_NAME" ]] || ENV_NAME="${MODULE}-env"
[[ -n "$OUT_YML"  ]] || { mkdir -p "$DEFAULT_OUT_DIR"; OUT_YML="$DEFAULT_OUT_DIR/${MODULE}.yml"; }

# Ensure base files exist (empty OK)
mkdir -p "$(dirname "$BASE_CONDA")" "$(dirname "$BASE_PIP")"
touch "$BASE_CONDA" "$BASE_PIP"

log "[+] Module:     $MODULE"
log "[+] Wrapper:    $WRAPPER"
log "[+] Out YAML:   $OUT_YML"
log "[+] Env name:   $ENV_NAME"
log "[+] Channels:   $CHANNELS"
log "[+] Base conda: $BASE_CONDA"
log "[+] Base pip:   $BASE_PIP"

# ---------- Helpers: unique, set minus ----------
uniq_sorted() { awk 'NF{print tolower($0)}' | sed 's/[[:space:]]\+$//' | sort -u; }
set_minus()    { comm -23 <(printf '%s\n' "$1" | uniq_sorted) <(printf '%s\n' "$2" | uniq_sorted); }

# ---------- Curated mappings ----------
# Map Python import -> conda package name (extend as needed)
declare -A PY_TO_CONDA=(
  [pysam]=pysam
  [numpy]=numpy
  [pandas]=pandas
  [scipy]=scipy
  [pyyaml]=pyyaml
  [biopython]=biopython
  [matplotlib]=matplotlib
  [pyspark]=pyspark
  [requests]=requests
  [tqdm]=tqdm
  [pyarrow]=pyarrow
  [pydantic]=pydantic
  [click]=click
)

# Map CLI tool -> conda package (extend as needed)
declare -A BIN_TO_CONDA=(
  [bwa]=bwa
  [samtools]=samtools
  [bcftools]=bcftools
  [htslib]=htslib
  [minimap2]=minimap2
  [fastp]=fastp
  [bedtools]=bedtools
  [pigz]=pigz
  [seqtk]=seqtk
  [tabix]=htslib
)

# ---------- Extract: inline Python imports ----------
extract_python_imports() {
  # Capture lines inside any <<'PYCODE' ... PYCODE heredocs and pull first-level imports
  awk '
    $0 ~ /<<\047PYCODE\047/ {inpy=1; next}
    inpy && $0 ~ /^PYCODE$/ {inpy=0; next}
    inpy {print}
  ' "$WRAPPER" \
  | awk '
      $1=="import" { 
        # import a, b as c  -> a b
        gsub(/[,]/," "); 
        for (i=2;i<=NF;i++){ 
          gsub(/^as$/,"", $i); 
          if ($i!="" && $i!="as") { split($i,parts,"."); print parts[1]; } 
        } 
      }
      $1=="from" && NF>=3 { 
        split($2,parts,"."); print parts[1]; 
      }
    ' \
  | uniq_sorted
}

# ---------- Extract: R libraries ----------
extract_r_libs() {
  grep -Eo 'library\([^)]+\)|require\([^)]+\)' "$WRAPPER" 2>/dev/null \
    | sed -E 's/.*\(["'\'' ]*([^"'\'' )]+).*/\1/i' \
    | uniq_sorted
}

# ---------- Extract: CLI tools ----------
extract_bins() {
  # Heuristic: any bare command tokens that match BIN_TO_CONDA keys
  # Look at executable lines (ignore comments, functions, case labels)
  grep -Ehv '^[[:space:]]*#|^[[:space:]]*$|^[[:space:]]*(function|case|esac|;;)' "$WRAPPER" \
    | grep -Eo '\b[a-zA-Z0-9._+-]+\b' \
    | while read -r tok; do
        if [[ -n "${BIN_TO_CONDA[$tok]:-}" ]]; then echo "$tok"; fi
      done \
    | uniq_sorted
}

# ---------- Build candidate sets ----------
PY_IMPORTS="$(extract_python_imports || true)"
R_LIBS="$(extract_r_libs || true)"
BINS="$(extract_bins || true)"

log "[+] Python imports: $(printf '%s' "$PY_IMPORTS" | tr '\n' ' ')"
log "[+] R libraries:    $(printf '%s' "$R_LIBS" | tr '\n' ' ')"
log "[+] CLI tools:      $(printf '%s' "$BINS" | tr '\n' ' ')"

# ---------- Map Python imports -> conda or pip ----------
map_py_to_conda_or_pip() {
  local pkgs_conda=() pkgs_pip=()
  while read -r imp || [[ -n "${imp:-}" ]]; do
    [[ -z "${imp:-}" ]] && continue
    # forced routes
    for f in "${ASSUME_PIP[@]:-}"; do [[ "$imp" == "$f" ]] && { pkgs_pip+=("$imp"); continue 2; }; done
    for f in "${FORCE_CONDA[@]:-}"; do [[ "$imp" == "$f" ]] && { pkgs_conda+=("$f"); continue 2; }; done
    # curated mapping or heuristic
    if [[ -n "${PY_TO_CONDA[$imp]:-}" ]]; then
      pkgs_conda+=("${PY_TO_CONDA[$imp]}")
    else
      # heuristic: many python libs are conda packages with same/lower-dashed name
      cand="${imp//_/-}"
      pkgs_conda+=("$cand")
    fi
  done <<< "$PY_IMPORTS"

  printf '%s\n' "__CONDALIST__"
  printf '%s\n' "${pkgs_conda[@]:-}" | uniq_sorted
  printf '%s\n' "__PIPLIST__"
  printf '%s\n' "${pkgs_pip[@]:-}" | uniq_sorted
}

# Run mapping
MAP_OUT="$(map_py_to_conda_or_pip)"
PY_CONDA="$(printf '%s\n' "$MAP_OUT" | awk '/^__CONDALIST__/{f=1;next} /^__PIPLIST__/{exit} f' | uniq_sorted)"
PY_PIP="$(printf '%s\n' "$MAP_OUT" | awk '/^__PIPLIST__/{f=1;next} f' | uniq_sorted)"

# R -> conda package prefix r-
R_CONDA="$(printf '%s\n' "$R_LIBS" | awk 'NF{print "r-"$0}' | uniq_sorted)"
# Binaries -> conda via mapping
BIN_CONDA="$(printf '%s\n' "$BINS" | while read -r b; do [[ -n "$b" ]] && echo "${BIN_TO_CONDA[$b]}"; done | uniq_sorted)"

# ---------- Subtract base (conda & pip) ----------
BASE_CONDA_SET="$(awk 'NF{print tolower($0)}' "$BASE_CONDA" 2>/dev/null || true)"
BASE_PIP_SET="$(awk 'NF{print tolower($0)}' "$BASE_PIP" 2>/dev/null || true)"

ALL_CONDA="$(printf '%s\n%s\n%s\n' "$PY_CONDA" "$R_CONDA" "$BIN_CONDA" | uniq_sorted)"
KEEP_CONDA="$(set_minus "$ALL_CONDA" "$BASE_CONDA_SET")"

KEEP_PIP="$(set_minus "$PY_PIP" "$BASE_PIP_SET")"

# Never duplicate python itself if base has it; else include a python line
PYTHON_LINE="python$ENV_PY_VERSION"
if printf '%s\n' "$BASE_CONDA_SET" | grep -q '^python\>'; then
  PYTHON_LINE=""  # base already has python
fi

log "[+] Conda deps (post-base): $(printf '%s' "$KEEP_CONDA" | tr '\n' ' ')"
log "[+] Pip deps   (post-base): $(printf '%s' "$KEEP_PIP" | tr '\n' ' ')"

# ---------- Write YAML ----------
mkdir -p "$(dirname "$OUT_YML")"
{
  echo "name: $ENV_NAME"
  IFS=, read -r -a ch <<<"$CHANNELS"
  printf 'channels: ['
  local first=1
  for c in "${ch[@]}"; do
    [[ $first -eq 1 ]] && printf '%s' "$c" || printf ', %s' "$c"
    first=0
  done
  echo "]"
  echo "dependencies:"
  if [[ -n "$PYTHON_LINE" ]]; then
    echo "  - $PYTHON_LINE"
  fi
  if [[ -n "$KEEP_CONDA" ]]; then
    printf '%s\n' "$KEEP_CONDA" | sed 's/^/  - /'
  fi
  if [[ -n "$KEEP_PIP" ]]; then
    echo "  - pip"
    echo "  - pip:"
    printf '%s\n' "$KEEP_PIP" | sed 's/^/    - /'
  fi
} > "$OUT_YML"

log "[+] Wrote $OUT_YML"

# ---------- Optional dry-run test ----------
if [[ $DRY_RUN -eq 1 ]]; then
  if command -v mamba >/dev/null 2>&1; then
    log "[+] Testing with mamba (dry-run)..."
    mamba create -n "$ENV_NAME" -f "$OUT_YML" --dry-run || {
      log "[!] mamba dry-run reported problems"; exit 1;
    }
    log "[+] Dry-run OK."
  else
    log "[!] mamba not found; skipping --test."
  fi
fi
