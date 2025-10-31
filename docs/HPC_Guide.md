# 🚀 MalariAPI HPC Integration

This document explains how to set up, install, and use the **MalariAPI HPC module** — a fully integrated extension that allows users to **submit, monitor, and retrieve jobs from remote computing clusters (HPCs)** such as UVA Rivanna or other institutional clusters.

---

## 📘 Overview

MalariAPI can now interact directly with high-performance computing systems (HPCs) through a single command-line interface:

```bash
mapi <hpc-name> <command> [options]
```

This means once your HPC connection is initialized, you can run jobs remotely, view logs, sync data, and manage results **without leaving your local environment**.

All HPC configuration, mirroring, and communication are handled through:

```
~/MalariAPI/tools/
├── hpc_install            ← first-time installer
├── hpc_lib                ← shared helper library
├── hpc_templates/         ← command templates copied per HPC
├── <hpc-name>/            ← generated command set for that HPC
└── hpc_config.json        ← your cluster configuration file
```

---

## ⚙️ System Requirements

Before installing, ensure you have the following tools on your **local system** (Ubuntu or WSL is fine):

| Tool | Purpose | Install Command |
|------|----------|-----------------|
| `ssh` | Secure connection to HPC | included in most systems |
| `rsync` | File synchronization | `sudo apt install rsync` |
| `jq` | JSON parsing | `sudo apt install jq` or `mamba install -c conda-forge jq` |
| `sshfs` | Mount remote directories locally | `sudo apt install sshfs` |

> ⚠️ If you don’t have sudo access, all of these can be installed in your `base` Conda environment (`mamba install -c conda-forge jq rsync openssh`).

---

## 🧰 What the Installer Does

The installer (`tools/hpc_install`) automates **everything** about setting up your HPC connection.

When you run:

```bash
~/MalariAPI/tools/hpc_install <hpc-name> <user@cluster.address>
```

it will:

1. **Create or update** your configuration file `tools/hpc_config.json`.
2. **Generate an SSH key** if you don’t already have one.
3. **Add your key to the remote cluster** (`ssh-copy-id`).
4. **Mirror your local MalariAPI repository** to:
   - Remote **HOME** → `~/MalariAPI`
   - Remote **SCRATCH** → `/scratch/<user>/MalariAPI/scratch`
5. **Mount both locations locally** (if `sshfs` is available):
   ```
   ~/MalariAPI/_hpc_home/
   ~/MalariAPI/_hpc_scratch/
   ```
6. **Create a command set** under `~/MalariAPI/tools/<hpc-name>/` that you can call through `mapi`.

---

## 🧩 Folder Architecture

```text
~/MalariAPI/
├── bin/mapi
├── tools/
│   ├── hpc_install
│   ├── hpc_lib
│   ├── hpc_config.json
│   ├── hpc_templates/
│   ├── <hpc-name>/
│   ├── _hpc_home/
│   └── _hpc_scratch/
└── scratch/
```

---

## 🧭 Understanding File Locations

| Space | Example Path | Description |
|--------|---------------|-------------|
| **Local** | `~/MalariAPI/...` | Your personal copy of the code and small test data. |
| **Remote HOME** | `~/MalariAPI/...` on HPC | A mirrored copy for job scripts, configs, and logs. |
| **Remote SCRATCH** | `/scratch/$USER/MalariAPI/...` on HPC | High-performance storage for large inputs and job outputs. |

**Rules:**
- Code and small data live in **Local** and **Remote HOME**.
- Large data and outputs live in **Remote SCRATCH**.
- Syncing happens **one way**: Local → Remote.
- Pull results manually if needed.

---

## 🧮 Commands (via `mapi <hpc-name> <cmd>`)

| Command | Description |
|----------|-------------|
| `submit` | Submit a job remotely. |
| `status` | Check the status of a job. |
| `cancel` | Cancel jobs. |
| `look` | Tail the `.out` or `.err` log file. |
| `peak` | `ls` files on the remote scratch. |
| `pull` | Copy results from remote scratch → local scratch. |
| `push` | Manually copy data from local scratch → remote scratch. |
| `sync` | Mirror everything per the specs (code + scratch inputs). |

### `submit` Syntax

```bash
mapi rivanna submit   --name testJob   --env-name fastp   --cmd 'fastp -i [local]/scratch/data/R1.fq.gz -I [local]/scratch/data/R2.fq.gz -o [scratch]/scratch/out/clean_R1.fq.gz'   --force-local
```

**Path tokens:**
| Token | Meaning |
|--------|----------|
| `[local]` | Expands locally before submission. |
| `[home]` | Refers to `~/MalariAPI` on the HPC. |
| `[scratch]` | Refers to `/scratch/$USER/MalariAPI` on the HPC. |
| *(no tag)* | Treated as a local path. |

**`--force-local`**: After the job completes, automatically copy back `/scratch/$USER/MalariAPI/scratch/out/` → `~/MalariAPI/scratch/out/`.

---

## 🧱 Example Workflow

```bash
~/MalariAPI/tools/hpc_install rivanna njb8sg@login.hpc.virginia.edu
mapi rivanna sync
mapi rivanna submit   --name demo --env-name fastp   --cmd 'echo "Running on $(hostname)" > [scratch]/scratch/out/test.txt'
mapi rivanna peak scratch/out
mapi rivanna look --name demo --jobid 123456 --out
mapi rivanna pull --from /scratch/$USER/MalariAPI/scratch/out --to ~/MalariAPI/scratch/out
```

---

## 🧩 Troubleshooting

| Issue | Explanation | Fix |
|--------|--------------|-----|
| `Missing: jq` | `jq` not installed locally | `sudo apt install jq` |
| `rsync: mkdir ... failed` | Scratch directory missing | Re-run installer |
| Password prompts after setup | Key not installed | Run `ssh-copy-id <user@cluster>` manually |
| `sshfs` mount fails | Missing package | `sudo apt install sshfs` |
| Job output missing locally | Output remains on HPC | Run `mapi <hpc> pull` |
| Reinstall or rename HPC | Safe to rerun installer | No duplication issues |

---

## 💡 Why Three Places?

| Location | Purpose | Sync Direction |
|-----------|----------|----------------|
| **Local (~)** | Working copy | source of truth |
| **Remote HOME (~)** | Mirror of code | Local → Remote |
| **Remote SCRATCH (/scratch/$USER)** | Large files & results | Pull manually |

---

## 🧭 Navigation Cheat Sheet

| Purpose | Command |
|----------|----------|
| View remote HPC home | `ls ~/MalariAPI/_hpc_home` |
| View remote HPC scratch | `ls ~/MalariAPI/_hpc_scratch` |
| Pull results | `mapi rivanna pull --from /scratch/$USER/MalariAPI/scratch/out` |
| Sync code | `mapi rivanna sync` |
| Submit job | `mapi rivanna submit --cmd "<...>"` |
| Cancel job | `mapi rivanna cancel <jobid>` |
| Watch logs | `mapi rivanna look --name <job> --jobid <id> --out` |

---

## 🧠 Philosophy

MalariAPI’s HPC extension is designed to be:
- **Simple** – single setup step, single command interface.
- **Secure** – key-based SSH only, no passwords stored.
- **Reproducible** – mirrored code, fixed environments.
- **Flexible** – works with any SSH-accessible Slurm cluster.
- **Non-destructive** – never deletes scratch data automatically.

---

## ✅ Summary

| Feature | Description |
|----------|-------------|
| One-line setup | `~/MalariAPI/tools/hpc_install <name> <user@host>` |
| Secure login | Installs and tests SSH key |
| Auto mirroring | Rsyncs code and inputs |
| Local mounts | View HPC dirs locally |
| Unified interface | `mapi <hpc> <cmd>` |
| Tag-aware paths | `[local]`, `[home]`, `[scratch]` |
| Safe sync | Local → Remote only |
| Works anywhere | Any SSH + Slurm HPC |

---

© 2025 Noah Brown (A-Crow-Nowhere)
