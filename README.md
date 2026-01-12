# Making 

**(MAPI)** is a lightweight, human-centered workflow framework for bioinformatics, designed to make complex analyses readable, reproducible, and accessible — locally, on HPC, and in low-resource environments.

MAPI is not just about running pipelines. It is about reducing friction at every stage of scientific computation.

---

   ### Set up base MAPI
   1. [Get your computer ready](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/setupyourcomputer.md)
   2. [Install base MAPI](docs/setup_MalariAPI.md)
   3. [Get familiar with MAPI's structure and workflow](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/GetFamiliarWithMapi.md)
   ### Set up advanced MAPI
   4. [Contributer user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/contributer.md)
   5. [Set up job subisison and sycing to a high performance computing (HPC) cluster](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/HPC_install.md)   
   6. [Creating modules and pipelines](docs/wrapping_tools.md)
   ### Helpful scripts
   7. [A guide on how to backup, and cleanup local distros of Ubuntu](docs/distro_backup.md)


## Why MAPI?

### Human-readable, predictable structure

MAPI enforces a clear and discoverable project layout,
This makes projects easier to understand, debug, and share — even months later.
To increase speed and de-clutter workflows, MAPI includes sever quality of life features.
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

No long paths. No guessing where results went.

## [Look at a minimal **MSF** workflow example here](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/MSF.example.md)

---

### Idempotent modules (less code, fewer mistakes)

MAPI replaces this with a single metadata header that declares:

- inputs and options  
- defaults and required arguments  
- outputs  
- environments  
- required resources  

From this header, MAPI automatically builds the argument parser, validates required inputs, documents usage, and standardizes outputs.

---

### Seamless HPC (High performance computer/cluster) passthrough

MAPI treats HPC clusters as first-class citizens.

Local and remote workflows use the same commands, and path tokens like `[scratch]` and `[home]` expand automatically.

[Look at a minimal **HPC** workflow example here](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/HPC.example.md)
---

### Git integration without Git expertise

MAPI includes safe Git wrappers designed for lab-based scientists (and computer scientists).

You can collaborate without becoming a Git expert.

---

### Backwards compatibility by design

Any executable script can be run as a module. Migration is gradual and opt-in.
All other organizational bennefits (besides MSFs) will still apply. 

---

### Designed for global and low-resource contexts

MAPI is suitable for teaching, collaboration, and institutions without dedicated DevOps support.

---


## In short

MAPI prioritizes clarity, reproducibility, and human usability.
