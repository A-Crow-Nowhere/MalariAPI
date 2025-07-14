


# Ubuntu WSL + MobaXterm Setup Guide for Windows

This guide walks you through setting up a fully functional Ubuntu WSL environment on a Windows machine and connecting to it using MobaXterm. It is suitable for new users and includes extra guidance for edge cases (like root login issues, SSH access, and future-proofing for Ubuntu updates).

---

## 🧰 Requirements

* Windows 10 (v2004+) or Windows 11
* Admin access on the system
* Internet connection

---

## 1. Install WSL and Ubuntu

### ✅ Enable WSL (Windows Subsystem for Linux)

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This will:

* Enable WSL
* Enable Virtual Machine Platform
* Download and install the latest Ubuntu LTS release (e.g., Ubuntu 22.04)

⚠️ If you already have WSL installed, you can manually install Ubuntu with:

```powershell
wsl --install -d Ubuntu
```

After installation, reboot if prompted.

### 🔄 First Launch

After installation, run:

```powershell
wsl
```

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

## 3. Configure Ubuntu WSL Environment

### 🔧 Update and Install Useful Tools

Inside WSL:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential git curl unzip htop openssh-server
```

### 👤 Ensure You’re Logged in as a Non-Root User

If WSL dropped you into a `root@host` shell without creating a user:

```bash
adduser yourname
usermod -aG sudo yourname
```

Then set WSL to auto-login as that user:

```bash
sudo nano /etc/wsl.conf
```

Add:

```ini
[user]
default=yourname
```

Restart WSL:

```powershell
wsl --shutdown
```

Then relaunch `wsl`.

---

## 4. Enable SSH Access for MobaXterm

Inside WSL:

```bash
sudo service ssh start
```

## 5. Connect with MobaXterm

1. Open MobaXterm
2. Click **Session → WSL → Ubuntu**

This will launch the WSL shell inside MobaXterm.
Pick a start up directory that makes it easy to access from your File Explorere

You're ready to move on to recreating the MalariAPI codebase.
