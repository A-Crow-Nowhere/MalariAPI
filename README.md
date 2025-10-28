# Welcome to MalariAPI (MAPI)
_An all in one, beginner friendly, and easy to use system for compuational investigations of malaria._ 
Because the organisms that cause malaria (and releated apicomplexan parasites) are so complex, tools that are designed for use in broader contexts are not neccisarily optimzed for malaria. Here, every module and pipeline has been pre-optimized for specific parameters, allowing for a precise and starndardized workflow. MAPI includes steps for setting up a ready-to-use bioninformatics workspace; easy enough for coding beginners, but robust and flexible enough to be implemented by verteran bioinformaticians. 
## Quick Links:
   1. [Setup a functional environment](docs/setup_MalariAPI.md) 
   2. [Install base MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/PopulateMAPI.md)
   3. [Beginners guide to navigating MAPI]()
   4. [Creating modules and pipelines](docs/wrapping_tools.md)
   5. [A guide on how to backup, and cleanup local distros of Ubuntu](docs/distro_backup.md)
   6. A list and description of the tools avalible in this Github repo

## Quick Install MalariAPI (MAPI)
### MAPI (MalariAPI) — Bash-first Modular Runner
```
## Quick install
bash <(curl -fsSL https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/setup/mapi/install.sh)
# or:
# wget -qO- https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/setup/mapi/install.sh | bash
```

```MalariAPI/
MalariAPI/
├── bin/					   # all executables (entrypoints live here)
│ ├── mapi					# main launcher (single file)
│ ├── modules/				# (NEW) flat module scripts + central YAMLs
│ │ ├── fastp.sh			# module executables live directly here
│ │ ├── bwa.sh				# e.g., "bwa.sh" runs the BWA module
│ │ ├── markdup.sh
│ │ └── yaml/				# all module metadata & parameter files
│ │ ├── fastp.yml
│ │ ├── bwa.yml
│ │ ├── markdup.yml
│ │ └── <name>.params.yml 
│ ├── packages/				# 3rd-party tools with MAPI shims
│ │ └── <tool_name>/
│ │ ├── <tool_name>			# required wrapper name (entrypoint)
│ │ ├── bin/				   # (optional) vendor binaries, if any
│ │ ├── package.yml			# metadata + env resolution
│ │ └── README.md
│ ├── pipelines/			   # mapi-conformant pipelines
│ │ └── <pipeline_name>/
│ │ ├── run.sh				   # required pipeline entrypoint
│ │ ├── pipeline.yml		   # required metadata & DAG hints
│ │ └── README.md
│ ├── scripts/				   # small helpers (non-conformant OK)
│ │ └── vcf_to_bed.sh (ex)
│ └── templates/			   # scaffolds for new modules/pipelines
│ ├── module_skel/
│ └── pipeline_skel/
├── envs/
│ └── yaml/		    		   # conda/mamba environments (full envs)
├── genomes/				   # (optional) references (gitignored)
├── scratch/				   # empty workspace (gitignored)
└── docs/		   			# user docs including this guide		
```
