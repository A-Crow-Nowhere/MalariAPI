# Welcome to MalariAPI (MAPI)
_An all in one, beginner friendly, and easy to use system for compuational investigations of malaria._ 
Because the organisms that cause malaria (and releated apicomplexan parasites) are so complex, tools that are designed for use in broader contexts are not neccisarily optimzed for malaria. Here, every module and pipeline has been pre-optimized for specific parameters, allowing for a precise and starndardized workflow. MAPI includes steps for setting up a ready-to-use bioninformatics workspace; easy enough for coding beginners, but robust and flexible enough to be implemented by verteran bioinformaticians. 
## Quick Links:
   ### Set up base MAPI
   1. [Setup a functional environment](docs/setup_MalariAPI.md) 
   2. [Install base MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/PopulateMAPI.md)
   3. [Beginner user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/user.md)
   ### Set up advanced MAPI
   4. [Contributer user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/contributer.md)
      
   6. [Creating modules and pipelines](docs/wrapping_tools.md)
   5. [A guide on how to backup, and cleanup local distros of Ubuntu](docs/distro_backup.md)
   6. A list and description of the tools avalible in this Github repo

## Quick Install MalariAPI (MAPI)
### MAPI (MalariAPI) — Bash-first Modular Runner
```
# Quick install
curl -fsSL https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/tools/install_mapi.sh -o install_mapi.sh
chmod +x install_mapi.sh
MINICONDA_HOME="$HOME/tools/miniconda3" ./install_mapi.sh
```

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
