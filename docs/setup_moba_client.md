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
