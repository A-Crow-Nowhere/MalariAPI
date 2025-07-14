# In this directory is a easy-to-use list of databses, tools, and code that are optimized for the investigation of the the malaria genome.
## It will include:
   1. Setup and install of Ubuntu and Mobaxterm for running scripts
   2. A general install of the basic packages included in MalariAPI
   3. A description on how to run the code
   4. A description on how to write the code to match the MalriAPI format
   5. A list and description of the tools avalible in this Github repo
   6. A graphical representation of where to find files in the repo

# Ubuntu WSL + MobaXterm Setup Guide for Windows

This guide walks you through setting up a fully functional Ubuntu WSL environment on a Windows machine and connecting to it using MobaXterm.
---

## 🧰 Requirements

* Windows 10 (v2004+) or Windows 11
* Admin access on the system
* Internet connection

---

## 1. Install WSL and Ubuntu

### ✅ Enable WSL (Windows Subsystem for Linux)

Open **PowerShell as Administrator** and run in powershell:

```powershell
wsl --install
```

This will:

* Enable WSL
* Enable Virtual Machine Platform
* Download and install the latest Ubuntu LTS release (e.g., Ubuntu 22.04)

⚠️ If you already have WSL installed, you can manually install Ubuntu with in powershell:

```powershell
wsl --install -d Ubuntu
```

After installation, reboot if prompted.

### 🔄 First Launch

After installation, run in powershell:

```powershell
wsl
```
*note it may do this automatically

This will launch Ubuntu and prompt:

```
Installing, this may take a few minutes...
Please create a default UNIX user account...
```

Enter a username and password (this will be your Linux user).

---

## 2. Install MobaXterm

Download the free edition of MobaXterm from:

📦 [https://mobaxterm.mobatek.net/download.html](https://mobaxterm.mobatek.net/download.html)

* Use the Installer Edition for automatic installation
* Launch MobaXterm

---


## 3. Enable SSH Access for MobaXterm

Inside WSL:

```bash
sudo service ssh start
```

## 4. Connect with MobaXterm

1. Open MobaXterm
2. Click **Session → WSL → Ubuntu**

This will launch the WSL shell inside MobaXterm.
Pick a start up directory that makes it easy to access from your File Explorere




# Configure MalariAPI Environment in mobaxterm

### 🔧 Update and Install Useful Tools

Inside WSL:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential git curl unzip htop openssh-server micro samtools bedtools
mkdir tools
mkdir bin
mkdir -p envs/yaml
mkdir /tools/miniconda3
cd tools/miniconda3

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O .
bash ../miniconda.sh -b -u -p .
rm ../miniconda.sh
source ~/tools/miniconda3/bin/activate
conda init --all


#For a variety of dependency issues samtools and bedtools will be installed in the highest level
#executable bin. 
#I will write any text editing commands with 'micro' but you can use vim or nano if you prefer. 
```
You're ready to move on to recreating the MalariAPI codebase.
