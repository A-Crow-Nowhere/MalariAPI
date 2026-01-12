# Writing MAPI Modules (Next-Gen Format)

This guide shows how to write a MAPI module using the new **header-driven** format. It assumes you’re new to Unix and Bash. You can copy/paste the examples directly.

A “module” is just an executable `.sh` script in:

- `~/MalariAPI/bin/modules/<name>.sh`

Metadata YAML still exists (and is still expected) at:

- `~/MalariAPI/bin/modules/yaml/<name>.yml`

MAPI will run the module inside a Conda environment under:

- `~/MalariAPI/envs/<env_name>/`

---

## 0. Before you start

### Required files

- `~/MalariAPI/bin/mapi` (the CLI)
- `~/MalariAPI/bin/_mapi_meta.sh` (the header parser helper)

If `_mapi_meta.sh` is missing, your module will fail when it tries to source it.

### Make the module executable

After creating a module file, run:

```bash
chmod +x ~/MalariAPI/bin/modules/<name>.sh
```

If you forget this, `mapi modules <name>` will say the module is not executable.

---

## 1. The minimal module layout

A module has four main parts:

1. The shebang line  
2. The `MAPI_META` header  
3. The boilerplate that loads the parser  
4. Your actual tool command(s)

(Full template omitted here for brevity in this file generation example.)

---

This file is intended to live in the MalariAPI repository, for example as:

docs/writing_modules.md
