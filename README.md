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
# Quick Install Guide (Local and HPC)

This document outlines a streamlined setup process for MalariAPI, covering both local installation and HPC configuration.

## 1. (Optional) Fork the Repository

If you plan to contribute:

1. Visit: https://github.com/A-Crow-Nowhere/MalariAPI
2. Select "Fork" to create your personal copy.

This lets you push to your own fork while still pulling updates from the main repository.

## 2. Clone MalariAPI to Your Local Machine

### If you forked:

```bash
git clone git@github.com:<your-github-username>/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI
```

### If you did not fork:

```bash
git clone https://github.com/A-Crow-Nowhere/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI
```

Replace `<your-github-username>` with your actual GitHub username.

## 3. Run the Local Installer

This installer will:

- Install or reuse a Miniconda distribution.
- Update the base environment using `envs/base.yml`.
- Add MalariAPI executables to your PATH.
- Configure git identity for this repository.
- Add an `upstream` remote if you cloned a fork.

### Recommended command:

```bash
./tools/install_mapi.sh   --repo-owner "<your-github-username>"   --upstream-owner "A-Crow-Nowhere"   --git-name "Your Name"   --git-email "you@example.edu"
```

If you cloned the canonical repo directly:

```bash
./tools/install_mapi.sh   --git-name "Your Name"   --git-email "you@example.edu"
```

After installation:

```bash
source ~/.bashrc    # or: source ~/.zshrc
mapi --help
```

## 4. Install MalariAPI on Your HPC

From your local clone:

```bash
cd ~/MalariAPI

./tools/hpc_install rivanna   --remote-user "<your-hpc-username>"   --remote-host "login.hpc.virginia.edu"
```

Replace `<your-hpc-username>` with your HPC account name.

Test the connection:

```bash
mapi rivanna status
```

A valid response confirms successful HPC setup.

## 5. Git Helper Commands in MAPI

MalariAPI provides simple wrappers for common git operations.

### Update your local branch from upstream:

```bash
mapi git update
```

### Upload your current branch to your fork:

```bash
mapi git upload
```

### Switch or create a branch:

```bash
mapi git switch feature/new_module
```

## 6. Next Steps

List available modules and pipelines:

```bash
mapi list
```

Create a development branch:

```bash
git checkout -b feature/<description>
```

Build any additional Conda environments as needed with MAPI’s environment tools.

## Installation Complete

You now have a fully functional local and HPC installation of MalariAPI.
Refer to the Contributing section for details on module and pipeline development.


