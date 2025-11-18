# 🧬 MalariAPI Setup — Cross-Platform Environment Guide

This document walks you through setting up a fully functional environment for **MalariAPI (MAPI)** on **Windows (via WSL + MobaXterm)** or **macOS/Linux (native shell)**.  
By the end, you’ll have:
- a working Ubuntu or macOS terminal,
- Miniconda installed in `~/tools/miniconda3`,
- base dependencies installed, and  
- folders ready for the MAPI framework.

---

## 🧰 Requirements

- **Windows 10 (v2004+) / Windows 11** or **macOS 11+**
- Admin access
- Internet connection

---

## 🪟 Option A — Windows (WSL + MobaXterm)

### 1️⃣ Enable WSL and Install Ubuntu

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This enables WSL 2 and installs the latest Ubuntu LTS (e.g. 22.04).  
If you already have WSL, you can manually install Ubuntu with:

```powershell
wsl --install -d Ubuntu
```

Reboot if prompted.

Launch Ubuntu for the first time:

```powershell
wsl
```

Follow the prompt to create a username and password.

---

### 2️⃣ Install MobaXterm (optional but recommended)

Download from [mobaxterm.mobatek.net/download.html](https://mobaxterm.mobatek.net/download.html).  
Use the **Installer Edition** → **Session → WSL → Ubuntu**.  
Pick a startup directory that’s easy to reach from File Explorer (e.g. `/home/<user>/MalariAPI`).

---

### 3️⃣ Enable SSH (only if you want to remote into WSL)

```bash
sudo service ssh start
```

---

## 🍏 Option B — macOS or Linux (Native Terminal)

### 1️⃣ Install Command Line Tools and Homebrew (macOS only)

```bash
xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install git wget curl
```

*(Linux users can skip this; `apt`/`dnf`/`yum` will suffice.)*

---

## ⚙️ Step 2 — Configure MalariAPI Environment

All following commands run **inside your Linux or macOS shell** (either Ubuntu WSL or native).

### 🔧 1. Update System and Install Core Tools

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential git curl unzip htop openssh-server micro samtools bedtools
```

*(macOS users with Homebrew can instead run  
`brew install git curl htop samtools bedtools micro`)*


✅ **Basic framework of MalariAPI is now set up!**  
Proceed to the next guide to install the **MAPI repo**, link the `bin/mapi` launcher, and start adding tools and pipelines.
