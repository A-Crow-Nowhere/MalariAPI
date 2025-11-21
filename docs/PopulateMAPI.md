# 🚀 Populate MAPI `bin`, `envs`, and `genomes` (Cross‑Platform)

This guide pulls the latest **MalariAPI** content from GitHub into your local checkout, and downloads the current **Pf3D7 reference genome** with indexes.

> Assumes you completed the base setup and you are inside your working tree at `~/MalariAPI` (or another path of your choice).

---

## 0) Verify you’re in the right folder

```bash
pwd
# Expect: /home/<you>/MalariAPI   (or /Users/<you>/MalariAPI on macOS)
```

---

## 1) Get the latest repo content (bin, templates, env YAMLs)

If you **haven’t cloned** the repo yet:

```bash
cd ~
git clone https://github.com/A-Crow-Nowhere/MalariAPI.git
cd MalariAPI
```

If you already cloned it earlier, just update:

```bash
cd ~/MalariAPI
git pull --rebase --autostash
```

This gives you the latest `bin/`, `envs/yaml/`, and templates.  
(If you prefer not to clone, you can fetch single files with `curl -L -o`, but cloning is simpler and keeps everything in sync.)

> Optional: symlink the launcher so `mapi` is on your PATH without copying files:
>
> ```bash
> mkdir -p "$HOME/bin"
> ln -sfn "$PWD/bin/mapi" "$HOME/bin/mapi"
> # or: sudo ln -sfn "$PWD/bin/mapi" /usr/local/bin/mapi
> ```

---

## 2) Create/refresh Conda environments from repo YAMLs

List the provided YAMLs:

```bash
ls -1 envs/yaml
```

Create an environment (example: `base.yml` + a tool env like `bwa.yml`):

```bash
# Loop through provided enbviornments
for y in envs/yaml/*.yml; do n=$(basename "${y%.yml}"); conda env list|awk '{print $1}'|grep -qx "$n" && conda env update -n "$n" -f "$y" --prune || conda env create -n "$n" -f "$y"; done

```

---

## 3) Download the **latest Pf3D7 reference** into `genomes/`

We’ll fetch the current NCBI RefSeq assembly **GCF_000002765.6** FASTA (Pf3D7) from the UCSC GenArk mirror and index it. This is convenient and consistent across platforms.

```bash
cd ~/MalariAPI/genomes
mkdir -p pf3d7 && cd pf3d7

# Download the primary FASTA (about ~25 MB compressed)
wget -c https://hgdownload.soe.ucsc.edu/hubs/GCF/000/002/765/GCF_000002765.6/GCF_000002765.6.fa.gz

# Uncompress
gunzip -f GCF_000002765.6.fa.gz

# Standardize name for tools
mv GCF_000002765.6.fa pf3d7_GCF000002765v6.fa
```

> Alternative sources (if you prefer ENSEMBL or NCBI Datasets CLI) can be used; just keep file naming consistent with your pipelines.

---

## 4) Build common indexes (samtools / bwa / picard)

```bash
# samtools FASTA index
samtools faidx pf3d7_GCF000002765v6.fa

# bwa index (activate your bwa env first if needed)
# conda activate bwa
bwa index pf3d7_GCF000002765v6.fa

# (Optional) Picard sequence dictionary for GATK-compatible tools
# conda activate samtools  # or any env containing picard
picard CreateSequenceDictionary   R=pf3d7_GCF000002765v6.fa   O=pf3d7_GCF000002765v6.dict
```

Results you should see in `~/MalariAPI/genomes/pf3d7/`:

```
pf3d7_GCF000002765v6.fa
pf3d7_GCF000002765v6.fa.fai
pf3d7_GCF000002765v6.fa.amb
pf3d7_GCF000002765v6.fa.ann
pf3d7_GCF000002765v6.fa.bwt
pf3d7_GCF000002765v6.fa.pac
pf3d7_GCF000002765v6.fa.sa
pf3d7_GCF000002765v6.dict   # (optional, Picard/GATK)
```

---

## 5) Quick sanity check

```bash
# Where is MAPI rooted?
mapi where

# Can the launcher list modules?
mapi list modules

# Verify genome index loads (no output means OK)
bwa mem pf3d7_GCF000002765v6.fa <(echo -e ">t\nACGT") /dev/null 2>/dev/null | head -n1
```

---

## 6) Notes & tips

- Keep `genomes/` **out of Git** — it’s large and machine-specific. Use `.gitignore`.
- If you keep multiple references, use subfolders like `genomes/pf3d7/`, `genomes/hg38/`, etc.
- When environments change, re‑export with `conda env export -n <env> --from-history > envs/yaml/<env>.yml` for clean, portable specs.
- Use `mapi help <module>` to view per‑module metadata from `bin/modules/<name>/meta/module.yml`.
- On macOS (Apple Silicon), prefer **mambaforge** or **conda-forge** packages to avoid cross‑arch issues.
