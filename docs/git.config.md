# GitHub Setup on Ubuntu (WSL)

This document describes how to set up Git and GitHub for the first time on an Ubuntu WSL distribution using SSH (recommended).

---

## Assumptions
- Ubuntu running under WSL
- A GitHub account already exists
- SSH authentication will be used (not HTTPS)

---

## 1. Install Git

```bash
sudo apt update
sudo apt install -y git
```

Verify installation:

```bash
git --version
```

---

## 2. Configure Global Git Identity

These values are attached to every commit you create.

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

Verify configuration:

```bash
git config --global --list
```

---

## 3. Generate an SSH Key (inside WSL)

IMPORTANT: Run this inside the WSL Ubuntu terminal, not Windows PowerShell.

```bash
ssh-keygen -t ed25519 -C "you@example.com"
```

When prompted:
- File location: press Enter to accept the default
- Passphrase: optional, but recommended

This creates:
- ~/.ssh/id_ed25519        (private key)
- ~/.ssh/id_ed25519.pub    (public key)

---

## 4. Start SSH Agent and Add Key

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

If you see:
```
Could not open a connection to your authentication agent
```
re-run the eval command and try again.

---

## 5. Add SSH Key to GitHub

Copy your public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Then on GitHub:
- Go to Settings → SSH and GPG keys
- Click “New SSH key”
- Title: something descriptive (e.g. WSL Ubuntu)
- Paste the key and save

---

## 6. Test GitHub SSH Connection

```bash
ssh -T git@github.com
```

Expected output:
```
Hi <username>! You've successfully authenticated, but GitHub does not provide shell access.
```

---

## 7. Ensure Repository Uses SSH (Not HTTPS)

Check remotes:

```bash
git remote -v
```

If the URL starts with https://github.com, switch it to SSH:

```bash
git remote set-url origin git@github.com:USERNAME/REPO.git
```

Example:

```bash
git remote set-url origin git@github.com:A-Crow-Nowhere/MalariAPI.git
```

---

## 8. First Push and Upstream Setup

For an existing GitHub repository:

```bash
git branch -M main
git push -u origin main
```

The -u flag sets the upstream so future git push and git pull commands work without extra arguments.

---

## 9. Recommended Global Defaults

Set default branch name:

```bash
git config --global init.defaultBranch main
```

Enable colored output:

```bash
git config --global color.ui auto
```

Optional short status alias:

```bash
git config --global alias.st "status -sb"
```

Usage:

```bash
git st
```

---

## 10. WSL-Specific Notes

Do not mix Git between Windows and WSL.
- Use Git inside WSL only
- Do not run Git on the same repository from Windows tools

Line endings (recommended for Linux-first repos):

```bash
git config --global core.autocrlf false
```

---

## 11. Sanity Check

Run these commands to confirm setup:

```bash
git status -sb
git branch -vv
git remote -v
ssh -T git@github.com
```

---

## Common Commands Reference

```bash
git status -sb
git branch -vv
git fetch origin
git reset --hard origin/main
git push --force-with-lease
```
