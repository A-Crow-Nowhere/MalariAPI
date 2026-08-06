# MalariAPI Setup — Cross-Platform Environment Guide

This document walks you through setting up a fully functional environment for MalariAPI (MAPI) on Windows (via WSL + MobaXterm), macOS (via UTM Ubuntu VM), or Linux. By the end, you will have:
- a working Ubuntu Linux shell,
- Miniconda installed in ~/MalariAPI/tools/miniconda3,
- base dependencies installed, and
- directories ready for the MAPI framework.

## Requirements

- Windows 10 (v2004+) / Windows 11 OR macOS 11+ OR Linux
- Administrator access
- Internet connection
- Approximately 20–40 GB of free disk space

# Option A — Windows (WSL + MobaXterm)

## 1. Enable WSL and Install Ubuntu

Open PowerShell as an administrator and run:

```
wsl --install
```

This installs WSL2 and the latest Ubuntu LTS. If you prefer to specify Ubuntu manually:

```
wsl --install -d Ubuntu
```

Reboot if prompted.

Launch Ubuntu:

```
wsl
```

Set up your Linux username and password.

## 2. Install MobaXterm (optional but recommended)

Download the installer from:
https://mobaxterm.mobatek.net/download.html

Inside MobaXterm:
- Session → WSL → Ubuntu
- Choose a convenient startup directory such as /home/<user>/MalariAPI

This provides improved terminal features, a GUI file browser, and easier copy/paste.

## 3. Enable SSH (optional)

```
sudo service ssh start
```

# Option C — Linux (Native)

If you are already using a Linux distribution, ensure the required tools are present:

```
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential git curl wget unzip htop openssh-server micro samtools bedtools
```

Adjust to your distribution’s package manager as necessary.

# Step 2 — Configure the MalariAPI Environment

Run the following inside your Linux shell (WSL Ubuntu, UTM Ubuntu, or native Linux).

## 1. Install core tools

```
sudo apt update && sudo apt upgrade -y
sudo apt install -y     build-essential git curl wget unzip htop openssh-server micro     samtools bedtools
```

At this point, your system is ready for the MAPI installer. The next steps include:
- cloning the MalariAPI repository (your fork),
- running the one-shot installer at tools/install_mapi.sh,
- enabling the mapi command in your PATH,
- configuring the base environment,
- optionally configuring HPC support,
- and beginning module or pipeline development.

