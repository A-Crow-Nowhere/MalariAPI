
---

## `MalariAPI/docs/env-yaml-contract.md`

```markdown
# MAPI Environment YAML Contract

Each tool/module owns a **small, separate** conda/mamba environment. YAMLs live in `~/envs/yaml/` and are auto-created the first time the module runs.

## Location & Naming

- File: `~/envs/yaml/<tool>.yaml`
- The YAML’s `name:` **must equal** the module’s `ENV_NAME`.
- The filename’s basename **must equal** the module’s `YAML_BASE`.

Example mapping:
- Module: `~/bin/mapi.d/fastp.sh`
  - `ENV_NAME="fastp-env"`
  - `YAML_BASE="fastp"`
- YAML: `~/envs/yaml/fastp.yaml` with `name: fastp-env`

## Minimal Schema

```yaml
name: yourtool-env                  # MUST match module's ENV_NAME
channels: [conda-forge, bioconda, defaults]
dependencies:
  - yourtool=1.*                    # pin at least major.minor
  # Add explicit CLIs the module calls:
  # - samtools=1.20
  # - pigz
  # If needed:
  # - python=3.11
  # - pip
  # - pip:
  #   - something==x.y.z
