# Executable Wrappers for Your Tools (`~/bin`)

Make your scripts runnable from anywhere with clean commands. This doc shows two patterns:

1) **Single one-off script** → a tiny wrapper in `~/bin/<name>`  
2) **Namespace tool (multiple scripts)** → a dispatcher in `~/bin/<tool>` that runs subcommands in `<tool>/bin/*`

Each includes a blank, copy-pasteable wrapper and an optional **environment activation** (conda/mamba/venv).

---

## Prerequisites

Ensure `~/bin` is on your `PATH`:

```bash
mkdir -p "$HOME/bin"
case ":$PATH:" in *":$HOME/bin:"*) ;; *) export PATH="$HOME/bin:$PATH" ;; esac
```

After adding new executables to `~/bin`, refresh the command cache:

```bash
hash -r   # (bash)    # or: rehash (zsh)
```

---

## A) Single One-Off Script Wrapper

**Use when:** you have one script you want to call by a short command (e.g., run `align` to execute `~/projects/aligner/run.sh`).

### 1) Blank wrapper (no environment)

Create `~/bin/<command>`, make it executable.

```bash
# ✂️ copy from here
cat > "$HOME/bin/<command>" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# Absolute path to the real script:
REAL_SCRIPT="$HOME/path/to/your/script.sh"

# If the script expects to run from its own folder, uncomment:
# SCRIPT_DIR="$(cd -- "$(dirname -- "$REAL_SCRIPT")" && pwd)"
# cd "$SCRIPT_DIR"

# Replace wrapper with the real script, forwarding all args
exec "$REAL_SCRIPT" "$@"
EOF
chmod +x "$HOME/bin/<command>"
hash -r
# ✂️ copy to here
```

**Usage**
```bash
<command> [args...]
```

### 2) Same wrapper with **environment activation** (choose one)

**Conda / Mamba**
```bash
# ✂️ copy
cat > "$HOME/bin/<command>" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
REAL_SCRIPT="$HOME/path/to/your/script.sh"

# Run inside a named environment
exec conda run -n <env-name> --no-capture-output "$REAL_SCRIPT" "$@"
# micromamba alternative:
# exec micromamba run -n <env-name> "$REAL_SCRIPT" "$@"
EOF
chmod +x "$HOME/bin/<command>"
hash -r
# ✂️
```

**Python venv**
```bash
# ✂️ copy
cat > "$HOME/bin/<command>" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
REAL_SCRIPT="$HOME/path/to/your/script.sh"

# Activate a Python virtualenv, then exec
. "$HOME/.venvs/<venv-name>/bin/activate"
exec "$REAL_SCRIPT" "$@"
EOF
chmod +x "$HOME/bin/<command>"
hash -r
# ✂️
```

---

## B) Multi-Script Tool Dispatcher (Git-style subcommands)

**Use when:** you have a tool folder with many executables, e.g.:

```
~/bin/<tool>/bin/        # if you keep tool folders under ~/bin
```

You’ll run subcommands like:
```bash
<tool> <subcommand> [args...]
# e.g.
breakpoints intergenic.sh --flag
blast make_db.sh --db mydb
```

### 1) Blank dispatcher (no environment)

Create `~/bin/<tool>`. Set `TOOL_ROOT` to match your layout.

```bash
# ✂️ copy
cat > "$HOME/bin/<tool>" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# Where your tool’s scripts live (pick ONE and delete the other):
# TOOL_ROOT="$HOME/tools/<tool>/bin"
TOOL_ROOT="$HOME/bin/<tool>/bin"

usage() {
  echo "Usage: $(basename "$0") <subcommand> [args...]"
  echo
  echo "Available subcommands:"
  # Show executables (+ also list .py helpers)
  if command -v find >/dev/null 2>&1; then
    find "$TOOL_ROOT" -maxdepth 1 -type f -perm -u+x -printf "  %f\n" 2>/dev/null | sort
    find "$TOOL_ROOT" -maxdepth 1 -type f -name "*.py" -printf "  %f\n" 2>/dev/null | sort
  else
    ls -1 "$TOOL_ROOT" 2>/dev/null | sort
  fi
  exit "${1:-0}"
}

[[ $# -lt 1 ]] && usage 1

sub="$1"; shift
cmd="$TOOL_ROOT/$sub"

if [[ -x "$cmd" ]]; then
  exec "$cmd" "$@"
elif [[ -f "$cmd" && "$cmd" == *.py ]]; then
  exec python3 "$cmd" "$@"
else
  echo "$(basename "$0"): unknown subcommand: $sub" >&2
  usage 1
fi
EOF
chmod +x "$HOME/bin/<tool>"
hash -r
# ✂️
```

**Usage**
```bash
<tool> <subcommand> [args...]
<tool>                 # lists available subcommands
```

### 2) Dispatcher **with environment activation** (choose one)

**Conda / Mamba**
```bash
# ✂️ copy
cat > "$HOME/bin/<tool>" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# TOOL_ROOT="$HOME/tools/<tool>/bin"
TOOL_ROOT="$HOME/bin/<tool>/bin"
ENV_NAME="<env-name>"

usage() {
  echo "Usage: $(basename "$0") <subcommand> [args...]"
  echo
  echo "Available subcommands:"
  ls -1 "$TOOL_ROOT" 2>/dev/null | sort
  exit "${1:-0}"
}

[[ $# -lt 1 ]] && usage 1
sub="$1"; shift
cmd="$TOOL_ROOT/$sub"

if [[ -x "$cmd" ]]; then
  exec conda run -n "$ENV_NAME" --no-capture-output "$cmd" "$@"
elif [[ -f "$cmd" && "$cmd" == *.py ]]; then
  exec conda run -n "$ENV_NAME" --no-capture-output python "$cmd" "$@"
else
  echo "$(basename "$0"): unknown subcommand: $sub" >&2
  usage 1
fi
EOF
chmod +x "$HOME/bin/<tool>"
hash -r
# ✂️
```

**Python venv**
```bash
# ✂️ copy
cat > "$HOME/bin/<tool>" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# TOOL_ROOT="$HOME/tools/<tool>/bin"
TOOL_ROOT="$HOME/bin/<tool>/bin"
VENV="$HOME/.venvs/<venv-name>"

[[ -d "$VENV" ]] || { echo "Missing venv: $VENV" >&2; exit 1; }
. "$VENV/bin/activate"

usage() {
  echo "Usage: $(basename "$0") <subcommand> [args...]"
  echo
  echo "Available subcommands:"
  ls -1 "$TOOL_ROOT" 2>/dev/null | sort
  exit "${1:-0}"
}

[[ $# -lt 1 ]] && usage 1
sub="$1"; shift
cmd="$TOOL_ROOT/$sub"

if [[ -x "$cmd" ]]; then
  exec "$cmd" "$@"
elif [[ -f "$cmd" && "$cmd" == *.py ]]; then
  exec python "$cmd" "$@"
else
  echo "$(basename "$0"): unknown subcommand: $sub" >&2
  usage 1
fi
EOF
chmod +x "$HOME/bin/<tool>"
hash -r
# ✂️
```

---

## Optional: Bash completion for a dispatcher

Add lightweight tab-completion for subcommands:

```bash
# ✂️ copy
cat > "$HOME/.bash_completion_<tool>" <<'EOF'
_<tool>_complete() {
  local cur="${COMP_WORDS[COMP_CWORD]}"
  local subdir="$HOME/bin/<tool>/bin"   # or "$HOME/tools/<tool>/bin"
  COMPREPLY=( $(compgen -W "$(command ls -1 "$subdir" 2>/dev/null)" -- "$cur") )
}
complete -F _<tool>_complete <tool>
EOF

# Load now and persist (add `source` to ~/.bashrc)
source "$HOME/.bash_completion_<tool>" 2>/dev/null || true
# ✂️
```

---

## Tips & Gotchas

- **Name collisions:** If two tools export the same subcommand name, adjust `PATH` order or prefix names in your tool’s `bin/` (e.g., `toolx-…`).
- **Use absolute paths** in wrappers for reliability. If you move the tools, update `REAL_SCRIPT`/`TOOL_ROOT`.
- **Working directory:** If your scripts expect their own folder, `cd` there before `exec` (see comments above).
- **WSL & Windows Terminal:** To open directly into your Linux home, set the profile’s starting directory to:
  ```
  \\wsl.localhost\<YourDistro>\home\<your-user>
  ```

---

## Quick Examples

**One-off**
```bash
# Run your aligner anywhere
align reads.fq --ref ref.fa
```

**Dispatcher**
```bash
# Namespaced tool with subcommands
breakpoints intergenic.sh input.tsv > out.tsv
blast make_db.sh --db mydb
```

You now have clean, portable commands for both small scripts and larger tool suites—with optional environment activation baked in.
