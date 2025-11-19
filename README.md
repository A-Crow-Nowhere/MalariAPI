# Welcome to MalariAPI (MAPI)
_An all in one, beginner friendly, and easy to use system for compuational investigations of malaria._ 
Because the organisms that cause malaria (and releated apicomplexan parasites) are so complex, tools that are designed for use in broader contexts are not neccisarily optimzed for malaria. Here, every module and pipeline has been pre-optimized for specific parameters, allowing for a precise and starndardized workflow. MAPI includes steps for setting up a ready-to-use bioninformatics workspace; easy enough for coding beginners, but robust and flexible enough to be implemented by verteran bioinformaticians. 
## Quick Links:
   ### Set up base MAPI
   1. [Set up a functional environment](docs/setup_MalariAPI.md) 
   2. [Install base MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/PopulateMAPI.md)
   3. [Beginner user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/user.md)
   ### Set up advanced MAPI
   4. [Contributer user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/contributer.md)
   5. [Set up job subisison and sycing to a high performance computing (HPC) cluster](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/HPC_Guide.md)   
   6. [Creating modules and pipelines](docs/wrapping_tools.md)
   ### Helpful scripts
   7. [A guide on how to backup, and cleanup local distros of Ubuntu](docs/distro_backup.md)

## Quick Install MalariAPI (MAPI)
## 🚀 Quick install (local + HPC)

These steps assume:

- You have **git**, **SSH**, and either **curl** or **wget** installed.   # See [initial setup](docs/setup_MalariAPI.md) 
- You have access to an HPC and a valid username there.                   # Optional, see [HPC setup](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/HPC_Guide.md)

> 🧠 If you're new to git: just follow the commands as written, replacing the
> placeholders (`<your-github-username>`, `<your-hpc-username>`, etc.).

---

### 1. Fork (optional but recommended)

If you plan to contribute code:

1. Go to the canonical repo:  
   https://github.com/A-Crow-Nowhere/MalariAPI
2. Click **Fork** and create your own copy under your GitHub account.

You’ll use your fork as `origin`, and the canonical repo as `upstream`.

---

### 2. Clone MalariAPI to your local machine

On your laptop / WSL / workstation:

```bash
# If you forked:
git clone git@github.com:<your-github-username>/MalariAPI.git ~/MalariAPI

# Or, if you just want a read-only clone of the canonical repo:
# git clone https://github.com/A-Crow-Nowhere/MalariAPI.git ~/MalariAPI

cd ~/MalariAPI


```MalariAPI/
├── bin
│   ├── mapi
│   ├── modules
│   │   ├── easteregg.sh
│   │   ├── summarize_bam.sh
│   │   └── yaml
│   ├── packages
│   ├── pipelines
│   │   └── yaml
│   ├── scripts
│   │   ├── bed_to_igv.sh
│   │   ├── bed_to_vcf.sh
│   │   └── vcf_to_bed.sh
│   └── templates
│       ├── module_template.sh
│       ├── module_template.yml
│       ├── package_wrapper_template
│       ├── pipeline_template.sh
│       └── pipeline_template.yml
├── envs
│   └── base.yml
├── genomes
├── package_depot
└── tools
    ├── gen_env.sh
    ├── install_mapi.sh
    └── validate.sh

```
