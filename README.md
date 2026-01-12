# Making 

**(MAPI)** is a lightweight, human-centered workflow framework for bioinformatics, designed to make complex analyses readable, reproducible, and accessible — locally, on HPC, and in low-resource environments.

MAPI is not just about running pipelines. It is about reducing friction at every stage of scientific computation.

---

   ### Set up base MAPI
   1. [Set up a functional environment](docs/setup_MalariAPI.md) 
   2. [Install base MAPI]
   3. [Beginner user guide for MAPI]
   ### Set up advanced MAPI
   4. [Contributer user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/contributer.md)
   5. [Set up job subisison and sycing to a high performance computing (HPC) cluster](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/HPC_install.md)   
   6. [Creating modules and pipelines](docs/wrapping_tools.md)
   ### Helpful scripts
   7. [A guide on how to backup, and cleanup local distros of Ubuntu](docs/distro_backup.md)


## Why MAPI?

### Human-readable, predictable structure

MAPI enforces a clear and discoverable project layout:

- Modules live in one place  
- Pipelines live in one place  
- Environments live in one place  
- Genomes and references live in one place  

You never have to guess:
- where outputs are written  
- which environment was used  
- how a result was produced  

This makes projects easier to understand, debug, and share — even months later.

---

### MAPI Sample Folders (MSFs): zero-path workflows

MAPI introduces **MAPI Sample Folders (MSFs)**:

```
mapi-sampleName/
├── sampleName-output/
├── input_files.fastq.gz
└── summary.txt
```

When working inside an MSF:

- You do not need to type file paths
- Outputs are automatically written to the correct backend folder
- File naming follows a consistent, machine-readable convention

Example:

```bash
mapi modules bwa --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz
```

No long paths. No guessing where results went.

---

### Header-driven modules (less code, fewer mistakes)

Traditional bioinformatics scripts require manual argument parsing and duplicated documentation.

MAPI replaces this with a single metadata header that declares:

- inputs and options  
- defaults and required arguments  
- outputs  
- environments  
- required resources  

From this header, MAPI automatically builds the argument parser, validates required inputs, documents usage, and standardizes outputs.

---

### Conda environments without the pain

Each module or pipeline declares the environment it needs.

MAPI automatically runs code inside the correct environment and never requires users to manually activate environments.

---

### Seamless HPC passthrough

MAPI treats HPC clusters as first-class citizens.

Local and remote workflows use the same commands, and path tokens like `[scratch]` and `[home]` expand automatically.

---

### Git integration without Git expertise

MAPI includes safe Git wrappers designed for scientists.

You can collaborate without becoming a Git expert.

---

### Resource awareness (genomes, indexes, models)

MAPI makes external dependencies explicit and reproducible by design.

---

### Backwards compatibility by design

Any executable script can be run as a module. Migration is gradual and opt-in.

---

### Designed for global and low-resource contexts

MAPI is suitable for teaching, collaboration, and institutions without dedicated DevOps support.

---

### A foundation for AI-assisted science

MAPI is designed to grow with AI, not be replaced by it.

---

## In short

MAPI prioritizes clarity, reproducibility, and human usability.
