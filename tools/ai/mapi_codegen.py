#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

def log(msg: str) -> None:
    print(f"[codegen] {msg}", file=sys.stderr)

def die(msg: str, rc: int = 2) -> None:
    print(f"[codegen] ERROR: {msg}", file=sys.stderr)
    raise SystemExit(rc)

def read_text(path: Optional[str], text: Optional[str]) -> str:
    if path and text:
        die("Provide only one of --in or --text (not both).")
    if path:
        p = Path(path)
        if not p.exists():
            die(f"Input file not found: {p}")
        return p.read_text(encoding="utf-8", errors="replace")
    if text is not None:
        return text
    die("Provide one of --in or --text.")
    raise AssertionError

PY_IMPORT_RE = re.compile(r"^\s*(?:from\s+([a-zA-Z0-9_\.]+)\s+import|import\s+([a-zA-Z0-9_\.]+))", re.M)
R_LIB_RE = re.compile(r"^\s*(?:library|require)\(\s*['\"]?([A-Za-z0-9_.]+)['\"]?\s*\)", re.M)
SHELL_CMD_RE = re.compile(r"^\s*([a-zA-Z0-9_.+-]+)\b")
ARGPARSE_FLAG_RE = re.compile(r"\.add_argument\(\s*(['\"]--?[A-Za-z0-9][A-Za-z0-9_-]*['\"])", re.M)
GETOPTS_RE = re.compile(r"getopts\s+['\"]([A-Za-z0-9:]+)['\"]")

def detect_language(path: Optional[str], code: str) -> str:
    if path:
        ext = Path(path).suffix.lower()
        if ext == ".py":
            return "python"
        if ext in (".r", ".rscript"):
            return "r"
        if ext in (".sh", ".bash"):
            return "bash"
    first = code.splitlines()[0] if code.splitlines() else ""
    if first.startswith("#!"):
        if "python" in first:
            return "python"
        if "Rscript" in first or re.search(r"\bR\b", first):
            return "r"
        if "bash" in first or "sh" in first:
            return "bash"
    if re.search(r"^\s*import\s+\w+|^\s*from\s+\w+\s+import", code, re.M):
        return "python"
    if re.search(r"^\s*(library|require)\(", code, re.M):
        return "r"
    if re.search(r"^\s*getopts\b|^\s*set\s+-e", code, re.M):
        return "bash"
    return "other"

def python_imports(code: str) -> List[str]:
    mods = set()
    for m in PY_IMPORT_RE.finditer(code):
        mod = m.group(1) or m.group(2) or ""
        top = mod.split(".")[0]
        if top and top not in {"__future__"}:
            mods.add(top)
    return sorted(mods)

def r_packages(code: str) -> List[str]:
    return sorted(set(R_LIB_RE.findall(code)))

def bash_commands(code: str) -> List[str]:
    cmds = set()
    for line in code.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", line):
            continue
        m = SHELL_CMD_RE.match(line)
        if not m:
            continue
        cmd = m.group(1)
        if cmd in {"if","then","fi","for","do","done","case","esac","while","until","function","{","}","export","local","set","echo","printf","test","["}:
            continue
        if cmd.startswith((".", "source")):
            continue
        cmds.add(cmd)
    return sorted(cmds)

def guess_inputs_from_argparse(code: str) -> List[Dict[str, Any]]:
    flags = sorted(set([s.strip("'\"") for s in ARGPARSE_FLAG_RE.findall(code)]))
    out: List[Dict[str, Any]] = []
    for fl in flags[:120]:
        out.append({
            "flag": fl,
            "type": "string",
            "required": False,
            "default": None,
            "description": "Auto-detected flag (argparse). Please edit description/types."
        })
    return out

def guess_inputs_from_getopts(code: str) -> List[Dict[str, Any]]:
    m = GETOPTS_RE.search(code)
    if not m:
        return []
    optstr = m.group(1)
    out = []
    for ch in optstr:
        if ch == ":":
            continue
        out.append({
            "flag": f"-{ch}",
            "type": "string",
            "required": False,
            "default": None,
            "description": "Auto-detected flag (getopts). Please edit."
        })
    return out

def minimal_spec(name: str, description: str, language: str) -> Dict[str, Any]:
    return {
        "name": name,
        "description": description,
        "language": language,
        "entrypoint": {"command": "", "notes": ""},
        "inputs": [],
        "outputs": [],
        "dependencies": {"cli_tools": [], "conda_packages": [], "python_pip": [], "r_packages": [], "notes": ""},
        "resources": {"cpus": 1, "mem_gb": 2, "time": "01:00:00", "gpu": False},
        "example": "",
        "notes": ""
    }

def cheap_validate(spec: Dict[str, Any]) -> None:
    req_top = ["name","description","language","entrypoint","inputs","outputs","dependencies"]
    for k in req_top:
        if k not in spec:
            raise ValueError(f"Spec missing required key: {k}")
    if not re.match(r"^[a-zA-Z0-9][a-zA-Z0-9._-]{1,63}$", spec["name"]):
        raise ValueError("Spec.name must be 2-64 chars: letters/numbers/._- and start with alnum")
    if spec["language"] not in {"python","r","bash","other"}:
        raise ValueError("Spec.language invalid")
    if "command" not in spec["entrypoint"]:
        raise ValueError("Spec.entrypoint.command missing")
    dep = spec["dependencies"]
    for k in ["cli_tools","conda_packages","python_pip","r_packages"]:
        if k not in dep or not isinstance(dep[k], list):
            raise ValueError(f"Spec.dependencies.{k} must be a list")

PROMPT = """You are generating a MAPI module SPEC as strict JSON only.
Rules:
- Output JSON only (no markdown, no commentary).
- Be conservative: only list dependencies you see strongly implied.
- If unsure, leave fields empty and explain uncertainty in spec.notes.
- Prefer flags you can confirm from argparse/add_argument or getopts.

MAPI layout (important):
- Command wrappers live in: bin/modules/<name>.sh
- Templates live in: tools/templates/
- Env YAMLs live in: modules/yaml/<name>.yml
- Runtime conda envs live in: envs/<name>
- Underlying tool code can be copied to: tools/generated/<name>/<name>.(sh|py|R)

Task:
Given CODE and STATIC_HINTS, fill the JSON spec.

STATIC_HINTS:
{static_hints}

CODE:
{code}
"""

def ollama_generate_spec(model: str, static_hints: str, code: str, name: str, language: str) -> Dict[str, Any]:
    if shutil.which("ollama") is None:
        raise RuntimeError("ollama not found in PATH. Install Ollama or use --backend offline.")
    prompt = PROMPT.format(static_hints=static_hints, code=code[:120000])
    log(f"LLM backend: ollama (model={model})")
    proc = subprocess.run(["ollama", "run", model], input=prompt.encode("utf-8"),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"ollama failed (rc={proc.returncode}): {proc.stderr.decode('utf-8', errors='replace')[:2000]}")
    out = proc.stdout.decode("utf-8", errors="replace").strip()
    out = re.sub(r"^\s*```(?:json)?\s*", "", out)
    out = re.sub(r"\s*```\s*$", "", out).strip()
    m = re.search(r"\{.*\}\s*$", out, re.S)
    if not m:
        raise RuntimeError("Model did not return a JSON object.")
    obj = json.loads(m.group(0))
    obj.setdefault("name", name)
    obj.setdefault("language", language)
    return obj

def shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"

def build_usage_block(inputs: List[Dict[str, Any]]) -> str:
    if not inputs:
        return "  (no auto-detected options; pass-through to underlying command)\n"
    lines = []
    for it in inputs[:200]:
        flag = it.get("flag","")
        desc = it.get("description","")
        req = " (required)" if it.get("required") else ""
        default = it.get("default", None)
        dflt = f" [default: {default}]" if default not in (None, "", "null") else ""
        lines.append(f"  {flag:<18} {desc}{req}{dflt}")
    return "\n".join(lines) + "\n"

def load_template(mapi_root: Path, language: str) -> str:
    tmpl_map = {
        "python": "module.python.tmpl",
        "r": "module.r.tmpl",
        "bash": "module.bash.tmpl",
        "other": "module.other.tmpl"
    }
    tmpl_name = tmpl_map.get(language, "module.other.tmpl")
    p = mapi_root / "tools" / "templates" / tmpl_name
    if not p.exists():
        die(f"Template not found: {p}")
    return p.read_text(encoding="utf-8")

def deps_block(spec: Dict[str, Any]) -> str:
    lang = spec.get("language")
    deps: List[str] = []
    if lang == "python":
        deps.append("  - python>=3.10")
        deps.append("  - pip")
        deps.extend([f"  - {p}" for p in spec["dependencies"].get("conda_packages", [])])
        pip_pkgs = spec["dependencies"].get("python_pip", [])
        if pip_pkgs:
            deps.append("  - pip:")
            deps.extend([f"      - {p}" for p in pip_pkgs])
    elif lang == "r":
        deps.append("  - r-base")
        deps.extend([f"  - {p}" for p in spec["dependencies"].get("conda_packages", [])])
        r_pkgs = spec["dependencies"].get("r_packages", [])
        deps.extend([f"  - r-{p.lower().replace('.', '-')}" for p in r_pkgs])
    elif lang == "bash":
        conda_pkgs = spec["dependencies"].get("conda_packages", [])
        deps.extend([f"  - {p}" for p in conda_pkgs] or ["  - bash"])
    else:
        deps.append("  - bash")
    cli = spec["dependencies"].get("cli_tools", [])
    if cli:
        deps.append("  # Suggested CLI tools (verify names):")
        deps.extend([f"  - {c}" for c in cli])
    return "\n".join(deps) + ("\n" if deps else "")

ENV_YAML_TEMPLATE = """name: {name}
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
{deps_block}
"""

def ensure_dirs(mapi_root: Path) -> None:
    needed = [
        "bin/modules",
        "bin/modules/yaml",
        "tools/templates",
        "tools",
        "tools/ai/spec_out",
        "tools/generated",
    ]
    for rel in needed:
        p = mapi_root / rel
        if not p.exists():
            log(f"Creating missing directory: {p}")
            p.mkdir(parents=True, exist_ok=True)

def copy_input_script_into_mapi(mapi_root: Path, name: str, in_path: str) -> Path:
    src = Path(in_path).resolve()
    if not src.exists():
        die(f"Cannot copy script; input does not exist: {src}")
    outdir = mapi_root / "tools" / "generated" / name
    outdir.mkdir(parents=True, exist_ok=True)

    ext = src.suffix
    if ext.lower() in (".bash",):
        ext = ".sh"
    if ext.lower() not in (".sh", ".py", ".r", ".R"):
        # keep original extension
        pass

    dest = outdir / f"{name}{ext}"
    shutil.copy2(src, dest)
    try:
        os.chmod(dest, 0o755)
    except Exception:
        pass

    log(f"Copied input script -> {dest}")
    return dest

def set_entrypoint_for_copied_script(spec: Dict[str, Any], copied: Path) -> None:
    name = spec["name"]
    rel = f"$MAPI_ROOT/tools/generated/{name}/{copied.name}"

    lang = spec.get("language", "other")
    if lang == "bash":
        cmd = f"bash \"{rel}\""
    elif lang == "python":
        cmd = f"python3 \"{rel}\""
    elif lang == "r":
        cmd = f"Rscript \"{rel}\""
    else:
        cmd = f"bash \"{rel}\""

    spec.setdefault("entrypoint", {})
    spec["entrypoint"]["command"] = cmd
    log(f"Entrypoint set -> {cmd}")


def render_wrapper(spec: Dict[str, Any], mapi_root: Path) -> Path:
    name = spec["name"]
    desc = spec.get("description","").strip().replace("\n"," ")
    usage_block = build_usage_block(spec.get("inputs", [])).rstrip("\n")

    entrypoint = spec.get("entrypoint", {}).get("command", "").strip()
    if not entrypoint:
        log("Decision: entrypoint.command empty -> defaulting to bash (no tool path known)")
        entrypoint = "bash"

    # We expect that for copied scripts, spec contains:
    # spec["_rel_tool"] = "tools/generated/<name>/<file>"
    rel_tool = spec.get("_rel_tool", "")
    if not rel_tool:
        log("WARNING: spec missing _rel_tool; wrapper may not be runnable until you set it.")
        rel_tool = "tools/generated/REPLACE_ME/REPLACE_ME.sh"

    lang = spec.get("language", "other")
    cmd_exe = {"bash": "bash", "python": "python3", "r": "Rscript"}.get(lang, "bash")

    def bash_q(s: str) -> str:
        return "'" + s.replace("'", "'\"'\"'") + "'"

    tmpl = load_template(mapi_root, lang)
    rendered = (
        tmpl.replace("{{NAME}}", name)
            .replace("{{DESCRIPTION}}", desc or "(no description)")
            .replace("{{USAGE_BLOCK}}", usage_block)
            .replace("{{ENTRYPOINT}}", entrypoint.replace("\n", " ").strip())
            .replace("{{REL_TOOL}}", rel_tool)
            .replace("{{CMD_EXE}}", bash_q(cmd_exe))
    )

    wrapper_path = mapi_root / "bin" / "modules" / f"{name}.sh"
    wrapper_path.parent.mkdir(parents=True, exist_ok=True)
    wrapper_path.write_text(rendered, encoding="utf-8")
    try:
        os.chmod(wrapper_path, 0o755)
    except Exception:
        pass
    return wrapper_path

def render_env_yaml(spec: Dict[str, Any], mapi_root: Path) -> Path:
    name = spec["name"]
    yml = ENV_YAML_TEMPLATE.format(name=name, deps_block=deps_block(spec).rstrip("\n"))
    yaml_path = mapi_root / "bin" / "modules" / "yaml" / f"{name}.yml"
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    yaml_path.write_text(yml, encoding="utf-8")
    return yaml_path

def write_spec_dump(spec: Dict[str, Any], mapi_root: Path) -> Path:
    outdir = mapi_root / "tools" / "ai" / "spec_out"
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / f"{spec['name']}.json"
    p.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    return p

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="mapi_codegen",
        description="Generate MAPI commands (bin/modules/*.sh) + env YAMLs (modules/yaml/*.yml) from code."
    )
    ap.add_argument("--in", dest="in_path", help="Input code file path (.py/.R/.sh)")
    ap.add_argument("--text", help="Raw code as a string (no file copy will occur)")
    ap.add_argument("--name", required=True, help="Command name (also env name). Example: count_files")
    ap.add_argument("--description", default="", help="Short human description (optional)")
    ap.add_argument("--out-root", default=".", help="MAPI root")
    ap.add_argument("--backend", default=os.environ.get("MAPI_CODEGEN_BACKEND", "offline"),
                    choices=["offline", "ollama"], help="Spec backend")
    ap.add_argument("--spec-in", help="Path to JSON spec (offline backend)")
    ap.add_argument("--ollama-model", default=os.environ.get("MAPI_CODEGEN_OLLAMA_MODEL", "llama3.1:8b-instruct"),
                    help="Ollama model name")
    ap.add_argument("--print-spec", action="store_true", help="Print computed JSON spec to stdout")
    args = ap.parse_args(argv)

    mapi_root = Path(args.out_root).resolve()
    log(f"out_root={mapi_root}")

    ensure_dirs(mapi_root)

    code = read_text(args.in_path, args.text)
    language = detect_language(args.in_path, code)
    desc = args.description or f"Auto-generated MAPI command for {args.name}"

    log(f"Input: {'file=' + args.in_path if args.in_path else 'text=<inline>'}")
    log(f"Detected language: {language}")
    log(f"Command name: {args.name}")
    log(f"Description: {desc}")
    log(f"Backend: {args.backend}")

    # Build spec
    if args.backend == "offline":
        if args.spec_in:
            log(f"offline mode: reading spec from --spec-in {args.spec_in}")
            spec = json.loads(Path(args.spec_in).read_text(encoding="utf-8"))
        else:
            log("offline mode: no spec provided -> building conservative static spec")
            spec = minimal_spec(args.name, desc, language)
    else:
        hints: List[str] = [f"detected_language={language}"]
        if language == "python":
            hints.append("python_imports=" + ", ".join(python_imports(code)))
            hints.append("argparse_flags=" + ", ".join([i["flag"] for i in guess_inputs_from_argparse(code)]))
        elif language == "r":
            hints.append("r_packages=" + ", ".join(r_packages(code)))
        elif language == "bash":
            hints.append("shell_commands=" + ", ".join(bash_commands(code)[:80]))
            hints.append("getopts_flags=" + ", ".join([i["flag"] for i in guess_inputs_from_getopts(code)]))
        static_hints = "\n".join(hints) + "\n"
        log("Prepared static hints for LLM.")
        spec = ollama_generate_spec(args.ollama_model, static_hints, code, args.name, language)

    # Normalize
    spec.setdefault("name", args.name)
    spec.setdefault("description", desc)
    spec.setdefault("language", language)
    spec.setdefault("entrypoint", {"command": "", "notes": ""})
    spec.setdefault("inputs", [])
    spec.setdefault("outputs", [])
    spec.setdefault("dependencies", {"cli_tools": [], "conda_packages": [], "python_pip": [], "r_packages": [], "notes": ""})
    spec.setdefault("resources", {"cpus": 1, "mem_gb": 2, "time": "01:00:00", "gpu": False})
    spec.setdefault("example", "")
    spec.setdefault("notes", "")

    # Backfill deps/inputs
    if spec["language"] == "python":
        if not spec["dependencies"].get("python_pip"):
            stdlib = {
                "argparse","sys","os","re","json","math","pathlib","typing","subprocess","shutil","datetime","time",
                "itertools","functools","collections","statistics","random","gzip","bz2","lzma"
            }
            imps = [m for m in python_imports(code) if m not in stdlib]
            if imps:
                spec["dependencies"]["python_pip"] = sorted(set(spec["dependencies"].get("python_pip", []) + imps))
                log(f"Backfill: python_pip inferred -> {spec['dependencies']['python_pip']}")
        if not spec["inputs"]:
            spec["inputs"] = guess_inputs_from_argparse(code)
            if spec["inputs"]:
                log(f"Backfill: inputs inferred from argparse -> {len(spec['inputs'])} flags")
    elif spec["language"] == "r":
        if not spec["dependencies"].get("r_packages"):
            pkgs = r_packages(code)
            if pkgs:
                spec["dependencies"]["r_packages"] = pkgs
                log(f"Backfill: r_packages inferred -> {pkgs}")
    elif spec["language"] == "bash":
        if not spec["dependencies"].get("cli_tools"):
            cmds = bash_commands(code)
            common = [c for c in cmds if c not in {"bash","sh","cat","echo","printf","grep","awk","sed","cut","sort","head","tail","tr","xargs","find"}]
            if common:
                spec["dependencies"]["cli_tools"] = common[:25]
                log(f"Backfill: cli_tools inferred -> {spec['dependencies']['cli_tools']}")

    # THIS IS THE CRITICAL NEW BEHAVIOR:
    # If --in was provided, copy the file into tools/generated/<name>/ and wire entrypoint to it.
    if args.in_path:
        copied = copy_input_script_into_mapi(mapi_root, spec["name"], args.in_path)
        spec["_rel_tool"] = f"tools/generated/{spec['name']}/{copied.name}"
        log(f"Recorded rel tool path -> {spec['_rel_tool']}")
        set_entrypoint_for_copied_script(spec, copied)
        
    else:
        log("No --in file provided, so no script was copied into tools/generated/. "
            "If entrypoint.command is empty, wrapper may contain a placeholder.")

    
    # Validate
    try:
        cheap_validate(spec)
        log("Spec validation: OK")
    except Exception as ex:
        die(f"Spec validation failed: {ex}")

    if args.print_spec:
        print(json.dumps(spec, indent=2))

    # Render outputs
    log("Rendering outputs...")
    wrapper_path = render_wrapper(spec, mapi_root)
    log(f"Wrote wrapper: {wrapper_path}")

    env_yaml_path = render_env_yaml(spec, mapi_root)
    log(f"Wrote env YAML: {env_yaml_path}")

    spec_path = write_spec_dump(spec, mapi_root)
    log(f"Wrote spec dump: {spec_path}")

    log("Done.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
