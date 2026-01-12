# HPC Workflows: Traditional vs MAPI

This document compares a **typical academic HPC workflow** with the same workflow using **MAPI’s streamlined HPC passthrough**.

The goal is not to claim one approach is wrong, but to show how MAPI reduces friction, repetition, and required background knowledge.

---

## Scenario

You want to:

1. Run a command on an HPC cluster
2. Write output to scratch
3. Check job status
4. Retrieve results back to your local machine

Example task:

```bash
echo "hello world" > out.txt
```

---

## Traditional HPC Workflow

### 1. Remember where you are connecting

```bash
ssh username@login.hpc.university.edu
```

You must remember:
- the login node
- your username
- which cluster you are on

---

### 2. Remember scratch vs home paths

```bash
cd /scratch/username/project/
```

Paths are cluster-specific and easy to mistype.

---

### 3. Create a job script (or inline sbatch)

```bash
nano job.sh
```

```bash
#!/usr/bin/env bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

echo "hello world" > /scratch/username/project/out.txt
```

Submit it:

```bash
sbatch job.sh
```

---

### 4. Check job status

```bash
squeue -u username
```

You must remember scheduler commands and flags.

---

### 5. Retrieve results manually

Exit the cluster:

```bash
exit
```

Copy files back:

```bash
scp username@login.hpc.university.edu:/scratch/username/project/out.txt .
```

If paths or usernames change, this breaks.

---

## MAPI Workflow

MAPI treats the HPC as an extension of your local environment.

### 1. Run the command from your local machine

```bash
mapi rivanna submit --cmd "echo hello world > [scratch]/out.txt"
```

What MAPI handles automatically:
- SSH connection
- Scheduler invocation
- Scratch path resolution

You never type:
- `/scratch/username`
- `ssh`
- `sbatch`

---

### 2. Check job status

```bash
mapi rivanna status
```

No scheduler flags required.

---

### 3. Retrieve results

```bash
mapi rivanna pull
```

Results appear locally under:

```
~/MalariAPI/scratch/
```

No `scp`. No guessing paths.

---

## Side-by-Side Comparison

| Step | Traditional HPC | MAPI |
|----|------------------|------|
Connect | `ssh login.node` | implicit |
Scratch paths | hard-coded | `[scratch]` token |
Job submission | `sbatch job.sh` | `mapi <hpc> submit` |
Job status | `squeue`, `sacct` | `mapi <hpc> status` |
File transfer | `scp` | `mapi <hpc> pull` |
Error surface | high | reduced |
Required HPC knowledge | high | low |

---

## What MAPI is not hiding

MAPI does **not** remove:
- scheduler limits
- queue times
- resource constraints

It removes:
- repetitive boilerplate
- memorization of paths
- manual file transfer

---

## Why this matters

In real research environments:
- HPC details differ between institutions
- students rotate frequently
- documentation goes out of date

MAPI centralizes HPC knowledge **once**, so users can focus on the science.

---

## When to use traditional workflows

Traditional workflows are still appropriate when:
- you are developing scheduler-specific optimizations
- you need fine-grained control over submission scripts
- you are debugging low-level performance issues

MAPI does not prevent this — you can always drop down to raw HPC commands.

---

## Summary

MAPI streamlines common HPC interactions while preserving transparency.

It reduces friction without removing control.
