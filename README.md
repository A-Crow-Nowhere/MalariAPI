# In this directory is a easy-to-use list of databses, tools, and code that are optimized for the investigation of the the malaria genome.
## It will include:
   1. [Setup and install of Ubuntu and Mobaxterm for running scripts](docs/setup_moba_client.md)
   2. [Setup the structure of MalariAPI](docs/setup_MalariAPI.md)
   3. [Install basic tools needed to run in MalariAPI]()
   4. [A description on how to wrap code to make it accessable from anywhere](docs/wrapping_tools.md)
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
├── README.md
│
├── setup/
│   └── mapi/
│       ├── install.sh
│       │
│       ├── bin/
│       │   ├── mapi                        # dispatcher (executable)
│       │   └── mapi.d/
│       │       ├── fastqc.sh
│       │       ├── fastp.sh
│       │       ├── bwa.sh
│       │       ├── lumpy.sh
│       │       ├── run.sh
│       │       ├── pipeline.sh -> run.sh   # symlink alias
│       │       └── doctor.sh               # optional diagnostic subcommand
│       │
│       ├── envs/
│       │   └── yaml/
│       │       ├── fastqc.yaml
│       │       ├── fastp.yaml
│       │       ├── bwa-mem2.yaml
│       │       ├── lumpy.yaml
│       │       └── yourtool.yaml           # template placeholder
│       │
│       ├── tools/
│       │   └── mapi/
│       │       └── lib.sh                  # shared helper functions
│       │
│       ├── templates/
│       │   ├── yourtool.sh                 # blank module template
│       │   └── yourtool.yaml               # blank env template
│       │
│       └── README_MAPI.md                  # usage + quick install guide
│
└── docs/
    ├── module-contract.md                  # required structure for modules
    ├── env-yaml-contract.md                # required structure for env YAMLs
    └── templates.md                        # reference + templates for new tools
```
