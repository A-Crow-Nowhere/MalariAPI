# Writing MAPI Modules (Next‑Gen Format)

This guide explains how to write a MAPI **module** using the new header‑driven format.
It is written for beginners who may have only used Unix/Linux a few times.

You can copy and paste examples directly from this document.

---

## What is a MAPI module?

A **module** is a single, executable shell script that performs one task
(e.g. align reads, summarize a BAM, filter variants).

Modules live here:

```
~/MalariAPI/bin/modules/<name>.sh
```

Each module:
- Has a metadata **header** that MAPI reads
- Is run inside a **conda environment**
- Can be executed:
  - via `mapi modules <name> ...`
  - or directly as a standalone script

---

## Directory structure (important)

```
MalariAPI/
├── bin/
│   ├── mapi
│   ├── _mapi_meta.sh
│   └── modules/
│       ├── example_module.sh
│       └── yaml/
│           └── example_module.yml
├── envs/
│   ├── default/
│   └── bwa/
└── genomes/
```

Required files:
- `bin/mapi` – main CLI
- `bin/_mapi_meta.sh` – header‑driven argument parser

If `_mapi_meta.sh` is missing, **modules will not run**.

---

## Make the module executable

After creating a module file:

```bash
chmod +x ~/MalariAPI/bin/modules/example_module.sh
```

If you forget this step, MAPI will say:
```
module not found or not executable
```

---

## The four parts of a module

Every module has four main sections:

1. **Shebang**
2. **MAPI_META header**
3. **Parser boilerplate**
4. **Your tool code**

---

## 1. Shebang (required)

Always start with:

```bash
#!/usr/bin/env bash
```

This tells the system to run the script with Bash.

---

## 2. The MAPI_META header

The header is where you describe:
- Inputs
- Outputs
- Options
- Environment
- Required resources

It starts and ends with:

```bash
# MAPI_META_BEGIN
...
# MAPI_META_END
```

Everything inside is **metadata**, not executable code.

---

## Full example header

```bash
# ======================================================================
# MAPI_META_BEGIN
# NAME: example_module
# SUMMARY: Demonstration module showing the MAPI format.
#
# ENV:
#   primary: default
#   fallback: default
#
# OUTPUTS:
#   - key: example_out
#     ext: txt
#     description: Main output file
#
# OPTIONS:
#   # ===================== INPUTS =====================
#   - flag: --in
#     short: -1
#     env: IN1
#     required: true
#     description: Primary input file
#
#   - flag: --in2
#     short: -2
#     env: IN2
#     required: false
#     description: Optional secondary input
#
#   # ===================== OUTPUT CONTROL =====================
#   - flag: --out
#     short: -o
#     env: OUTDIR
#     required: false
#     description: Output directory
#
#   # ===================== PERFORMANCE =====================
#   - flag: --threads
#     short: -t
#     env: THREADS
#     default: 8
#     required: false
#     description: Number of threads
#
# RESOURCES:
#   - name: example_reference
#     path: genomes/example/ref.fa
#     type: file
#     fetch: |
#       Download from:
#         https://example.org/ref.fa
# MAPI_META_END
# ======================================================================
```

---

## Header concepts explained

### ENV

```bash
# ENV:
#   primary: bwa
#   fallback: default
```

When run via MAPI:
1. Try `~/MalariAPI/envs/bwa`
2. If missing, fall back to `~/MalariAPI/envs/default`

If you run the module directly, you must already be in a working environment.

---

### OPTIONS (this replaces “inputs vs options”)

**Everything is an option.**

Required inputs are just options with:

```bash
required: true
```

Example:

```bash
- flag: --in
  env: IN1
  required: true
```

This means:
- User must supply `--in <value>`
- Inside the script, `$IN1` will contain that value

---

### env: VARIABLE_NAME

This line controls the variable name in your script:

```bash
env: IN1
```

Later in your code:

```bash
echo "$IN1"
```

If the names do not match, your script will break.

---

### Defaults

```bash
default: 8
```

If the user does not provide the flag, the variable is set to the default.

---

### Visual separators are allowed

These lines are **safe and ignored** by the parser:

```bash
# ===================== INPUTS =====================
```

They exist only for readability.

---

## 3. Parser boilerplate (required)

This code:
- Reads the header
- Builds an argument parser
- Fills variables automatically

```bash
set -euo pipefail

REPO_ROOT="${MAPI_REPO_ROOT:-"$(cd "$(dirname "$0")/../.." && pwd)"}"

# Load parser
# shellcheck disable=SC1090
source "$REPO_ROOT/bin/_mapi_meta.sh"

# Parse args
mapi_meta_load "$0"
mapi_meta_parse_args "$@"
```

After this runs:
- Variables like `IN1`, `IN2`, `THREADS`, `OUTDIR` exist
- Unknown flags go into `EXTRA[@]`

---

## 4. Output directory and naming

Standalone default behavior:

```bash
bn="$(basename "$IN1")"
sample_prefix="${bn%%.*}"
OUTDIR="${OUTDIR:-${sample_prefix}_mapi-out}"
mkdir -p "$OUTDIR"
```

When run inside a **MAPI Sample Folder** (`mapi-sampleName/`):
- `mapi` automatically sets `OUTDIR` to `sampleName-output/`

---

## Output filename helper

```bash
DEFAULT_TAG="example_out"

out_path(){
  printf "%s/%s.%s.%s\n" "$OUTDIR" "$sample_prefix" "$DEFAULT_TAG" "$1"
}
```

Example usage:

```bash
out_path txt
# -> sample-output/sample.example_out.txt
```

---

## Multiple outputs

You may list many outputs in the header:

```bash
# OUTPUTS:
#   - key: bam
#     ext: bam
#   - key: bai
#     ext: bai
```

This is documentation for now; your script decides what to write.

---

## 5. Resources (genomes, indexes, models)

Resources should live under:

```
~/MalariAPI/genomes/
```

Example check:

```bash
ensure_resources(){
  local ref="$REPO_ROOT/genomes/example/ref.fa"
  if [[ ! -f "$ref" ]]; then
    echo "[error] Missing reference: $ref" >&2
    exit 3
  fi
}
```

Call this before running your tool.

---

## 6. Tool execution example

```bash
tool_binary \
  --threads "$THREADS" \
  -i "$IN1" \
  ${IN2:+-I "$IN2"} \
  -o "$(out_path txt)" \
  "${EXTRA[@]}"
```

Explanation of the odd syntax:

```bash
${IN2:+-I "$IN2"}
```

- If `IN2` is empty → nothing is added
- If `IN2` is set → `-I <value>` is added

This avoids writing an explicit `if` statement.

---

## Common mistakes and fixes

### Module not executable
```
chmod +x bin/modules/<name>.sh
```

### Required option missing
- Header has `required: true`
- User forgot the flag
- Flag name mismatch

### Wrong variable name
- `env:` name does not match what you use in code

### Tool writes no output
- Forgot `mkdir -p "$OUTDIR"`
- Tool writes to current directory instead of `OUTDIR`

---

## Testing

### Run via MAPI (recommended)

```bash
mapi modules example_module --in input.txt
```

### Run standalone

```bash
bin/modules/example_module.sh --in input.txt
```

---

## Minimal working example (copy & test)

```bash
#!/usr/bin/env bash
# MAPI_META_BEGIN
# NAME: copy_file
# SUMMARY: Copies a file
# ENV:
#   primary: default
#   fallback: default
# OUTPUTS:
#   - key: copy
#     ext: txt
# OPTIONS:
#   - flag: --in
#     env: IN1
#     required: true
#   - flag: --out
#     env: OUTDIR
#     required: false
# MAPI_META_END

set -euo pipefail

REPO_ROOT="${MAPI_REPO_ROOT:-"$(cd "$(dirname "$0")/../.." && pwd)"}"
source "$REPO_ROOT/bin/_mapi_meta.sh"

mapi_meta_load "$0"
mapi_meta_parse_args "$@"

bn="$(basename "$IN1")"
sample_prefix="${bn%%.*}"
OUTDIR="${OUTDIR:-${sample_prefix}_mapi-out}"
mkdir -p "$OUTDIR"

DEFAULT_TAG="copy"
out_path(){ printf "%s/%s.%s.%s\n" "$OUTDIR" "$sample_prefix" "$DEFAULT_TAG" "$1"; }

cp "$IN1" "$(out_path txt)"
```
## TO SEE THE FULL CAPABILITIES OF MAPI MODULES [SEE THIS TECHNICAL EXAMPLE](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/Full_Example_module)
---

## Final advice

- Start simple
- Use `required: true` for real inputs
- Trust the header parser
- Test standalone and via `mapi`

This format is designed so **the header is the single source of truth**.
