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

---

### 🗂 2. Create Base Folders

```bash
mkdir -p ~/MalariAPI/{tools,genomes,bin,envs,scratch}
mkdir -p ~/MalariAPI/bin/{modules,pipelines,packages,scripts,templates}
mkdir -p ~/MalariAPI/bin/modules/yaml
mkdir -p ~/MalariAPI/bin/pipelines/yaml
```

This mirrors the MAPI repo layout.

---

### 🐍 3. Install Miniconda to `~/tools/miniconda3`

```bash
cd ~/MalariAPI/tools
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ~/MalariAPI/tools/miniconda3
rm miniconda.sh
source ~/MalariAPI/tools/miniconda3/bin/activate
conda init --all
```

> 🧩 On macOS, replace the `-Linux-x86_64.sh` URL with  
> `Miniconda3-latest-MacOSX-x86_64.sh`.

Restart your terminal after initialization.

---

### 🧬 4. Configure Conda Channels

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

---

### 📦 5. (Optional) Download the Base MAPI Environment YAML

```bash
cd ~/MalariAPI/envs
curl -L -o base.yml https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/envs/base.yml
```

This YAML defines the minimal Conda environment required for MAPI.

---

### 🧰 6. Verify Core Tools

```bash
conda --version
git --version
samtools --version
bedtools --version
```

All should return valid versions.

---

### 🪶 7. Notes

- `samtools` and `bedtools` are installed system-wide to avoid conflicts.
- `micro` is used as a lightweight terminal text editor (you can use `vim` or `nano` if preferred).
- Your full MAPI directory should now exist at:

```
~/MalariAPI/
├── tools/
│   └── miniconda3/
├── genomes/
├── bin/
│   └── modules/
│   │   └── yaml/
│   └── pipelines/
│   │   └── yaml/
│   └── packages/ 
│   └── scripts/
│   └── templates/
├── envs/
├── scratch/
```

---

✅ **Basic framework of MalariAPI is now set up!**  
Proceed to the next guide to install the **MAPI repo**, link the `bin/mapi` launcher, and start adding tools and pipelines.
