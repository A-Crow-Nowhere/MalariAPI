# MAPI Codegen MVP v2 (modules/ + templates/)

This bundle updates the previous MVP to match your current layout:
- NO `mapi.d/`
- Modules live under `modules/<name>/`
- Templates live under `templates/`
- AI implementation lives under `tools/ai/` (hidden from users)
- User runs everything via `mapi <command> ...`

## Files included

- `modules/codegen/run` : user-facing entrypoint (`mapi codegen ...`)
- `tools/ai/mapi_codegen.py` : generator (static analysis + optional LLM via Ollama)
- `tools/ai/spec_schema.json` : schema documentation for the spec
- `templates/module.{python,r,bash,other}.tmpl` : wrapper templates
- `examples/filter_tsv.py` : tiny example script
- `examples/make_codegen_demo.sh` : demo runner

## Install

Unzip into your `~/MalariAPI/` so it merges these folders:
- `modules/`
- `tools/ai/`
- `templates/`
- `examples/`

## Quick start (offline, no LLM)

From your MalariAPI root:

```bash
bash examples/make_codegen_demo.sh
```

## Using Ollama (optional)

```bash
mapi codegen   --in path/to/script.py   --name my_module   --out-root ~/MalariAPI   --backend ollama   --ollama-model llama3.1:8b-instruct
```

## Outputs generated

For module `<name>`:

- `modules/<name>/run`        (wrapper; pass-through args; loud logging)
- `modules/<name>/module.yml` (spec dump; JSON formatted)
- `tools/yaml/<name>.yml`     (conda env spec)

## Logging

Both `mapi codegen` and the generated wrappers print explicit messages about:
- inferred MAPI_ROOT
- inferred language and dependencies (when applicable)
- output paths written
- any fallback decisions (e.g., placeholder entrypoint)

