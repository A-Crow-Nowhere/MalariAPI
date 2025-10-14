# WSL Ubuntu Backup, Verify, and Restore Guide (PowerShell + WSL)

> Safely back up your entire WSL Ubuntu distro to a `.tar`, verify it, and (optionally) **rebuild/restore** from the tar to clean up virtual disk bloat — without needing Windows Hyper‑V tools. Everything here is copy‑pasteable.

---

## 0) Terminology & Notes

- **Export (`wsl --export`)** creates a `.tar` of *files inside* your Linux root. It does **not** include free space or Windows swap/pagefile.
- **Import (`wsl --import`)** creates a new distro with a fresh, right‑sized `ext4.vhdx` from your `.tar`.
- **Compact VHDX**: If you previously filled free space (e.g., with `/zero`) and the Windows disk grew, the cleanest fix is **export → unregister → import** (see §7). This avoids Windows‑side tools.
- Replace `Ubuntu` with your exact distro name from `wsl -l -v` where indicated.

---

## 1) Find your exact distro name

```powershell
wsl -l -v
# Note the NAME exactly (e.g., Ubuntu, Ubuntu-22.04, Debian, etc.)
$distro = 'Ubuntu'   # <-- replace with what you see
```

---

## 2) Choose a backup folder and create it

```powershell
$dest = "$env:USERPROFILE\Downloads\wsl_backups"   # change if you want
New-Item -ItemType Directory -Path $dest -Force | Out-Null
```

---

## 3) Cleanly shut down WSL and export

```powershell
wsl --shutdown
$stamp = Get-Date -Format yyyy-MM-dd_HH-mm
$tar   = Join-Path $dest "WSL-$($distro)-$stamp.tar"
wsl --export "$distro" "$tar"
```

**Verify the file exists and has size:**
```powershell
Get-Item "$tar" | Format-List Name,Length,LastWriteTime,FullName
```

---

## 4) Validate the archive (quick checks)

**List entries (structural read):**
```powershell
tar -tf "$tar" > $null
$LASTEXITCODE   # 0 means OK
```

**Peek a known file:**
```powershell
tar -xOf "$tar" etc/os-release | Select-Object -First 12
```

**Save a checksum for future re-checks:**
```powershell
Get-FileHash "$tar" -Algorithm SHA256 | Tee-Object -FilePath "$tar.sha256.txt"
```

---

## 5) (Optional) Create a smaller, compressed backup

**A. Compress after export (simple):**
```powershell
Compress-Archive -Path "$tar" -DestinationPath "$tar.zip"
```

**B. Stream-compress directly (saves space/time on slower disks):**
```powershell
wsl --shutdown
$stamp = Get-Date -Format yyyy-MM-dd_HH-mm
$gz    = Join-Path $dest "WSL-$($distro)-$stamp.tar.gz"
wsl --export "$distro" - | wsl -e gzip -c > "$gz"
# Validate gzip layer
wsl -e gzip -t "$gz"; $LASTEXITCODE
```

---

## 6) Strongest validation — test import safely

Imports as a **new** distro (won’t touch your current one).

```powershell
$inst = "$env:USERPROFILE\wsl_temp_import"
New-Item -ItemType Directory -Path $inst -Force | Out-Null

wsl --shutdown
wsl --import "Test-$distro" "$inst" "$tar" --version 2

wsl -d "Test-$distro" -- cat /etc/os-release
wsl -d "Test-$distro" -- ls -la /

wsl --unregister "Test-$distro"
Remove-Item -Recurse -Force "$inst"
```

---

## 7) **Recommended cleanup path**: Rebuild VHDX using export → import

> Use this when Windows reports low space after zero‑fill or heavy writes. This **does not** require Hyper‑V or Windows compaction tools.

```powershell
# 1) Export current distro (already done above, but do it again if needed)
$distro = 'Ubuntu'   # exact name
$dest   = "$env:USERPROFILE\Downloads\wsl_backups"
New-Item -ItemType Directory -Path $dest -Force | Out-Null
$stamp  = Get-Date -Format yyyy-MM-dd_HH-mm
$tar    = Join-Path $dest "WSL-$($distro)-$stamp.tar"

wsl --shutdown
wsl --export "$distro" "$tar"

# 2) Remove current distro (frees the old, bloated VHDX)
wsl --unregister "$distro"

# 3) Import from the tar to a fresh location (creates compact VHDX)
$newHome = "$env:USERPROFILE\WSL\$distro"
New-Item -ItemType Directory -Path $newHome -Force | Out-Null
wsl --import "$distro" "$newHome" "$tar" --version 2

# 4) First run sanity check
wsl -d "$distro" -- uname -a
```

> After import, the VHDX is right‑sized to your actual data. You may now move/zip the tar off C: to reclaim space.

---

## 8) Set default user and HOME after import (generalizable)

### A) Inside WSL — set your user’s HOME and make it default
```bash
# Replace 'noah' with your desired Linux username
USER=noah

# Ensure the home exists at /home/$USER and move content if needed
sudo mkdir -p /home/$USER
sudo usermod -d /home/$USER -m "$USER" 2>/dev/null || true

# Make WSL start as this user by default
printf "[user]\ndefault=%s\n" "$USER" | sudo tee /etc/wsl.conf

# (Optional) enable systemd & interop in /etc/wsl.conf
# printf "\n[boot]\nsystemd=true\n" | sudo tee -a /etc/wsl.conf
# printf "\n[interop]\nenabled=true\nappendWindowsPath=true\n" | sudo tee -a /etc/wsl.conf

exit
```

Back in PowerShell:
```powershell
wsl --shutdown
wsl -d "$distro"   # should open as your user in /home/<user>
```

### B) Windows Terminal — start in Linux HOME (UNC path)
- Settings → Ubuntu profile → **Starting directory**:
  ```text
  \\wsl.localhost\Ubuntu\home\<your-user>
  ```
- Or in `settings.json` for that profile:
  ```json
  "startingDirectory": "\\\\wsl.localhost\\Ubuntu\\home\\<your-user>"
  ```

> The UNC uses your distro name (`Ubuntu`, `Ubuntu-22.04`, etc.). Adjust accordingly.

---

## 9) Make your tar **generic** so others don’t inherit your username

### Option A (simple): default to **root** before export
```bash
# Inside WSL, as root or with sudo
printf "[user]\ndefault=root\n" | sudo tee /etc/wsl.conf
```
Export the distro. On import, recipients will land as root and can create their own user:
```powershell
# Recipient side (PowerShell)
wsl --import Ubuntu C:\WSL\Ubuntu C:\path\to\WSL-Ubuntu.tar --version 2
wsl -d Ubuntu -u root -- adduser newuser
wsl -d Ubuntu -u root -- usermod -aG sudo newuser
wsl -d Ubuntu -u root -- bash -lc 'printf "[user]\ndefault=newuser\n" > /etc/wsl.conf'
wsl --shutdown
```

### Option B (neutral image): remove your user before export
```bash
# Inside WSL, as root (WARNING: deletes that user's home!)
printf "[user]\ndefault=root\n" | sudo tee /etc/wsl.conf
sudo deluser --remove-home <your-user> || true
```
Export. Recipients import, then create their own user as above.

### Option C (guided): include a first‑run script
```bash
# /root/init-user.sh inside WSL
sudo tee /root/init-user.sh >/dev/null <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
read -rp "New Linux username: " U
adduser "$U"
usermod -aG sudo "$U"
printf "[user]\ndefault=%s\n" "$U" > /etc/wsl.conf
echo "User $U created and set as default. Run: wsl --shutdown, then reopen."
EOF
sudo chmod +x /root/init-user.sh
printf "[user]\ndefault=root\n" | sudo tee /etc/wsl.conf
```
Export. Recipient runs:
```powershell
wsl -d Ubuntu -u root -- /root/init-user.sh
wsl --shutdown
```

---

## 10) Optional: zero‑fill free space **inside Linux** (why/when)

**Why:** helps *compaction* or downstream imaging by turning free blocks into zeros.  
**Caution:** it fills to 100% temporarily. Prefer the export→import method in §7 for cleanup.

```bash
# Fills free space with zeros; shows progress; exits "successfully" thanks to || true
sudo dd if=/dev/zero of=/zero bs=1M status=progress || true
sync && sudo rm -f /zero
# Optional: also trim
sudo fstrim -av || true
```

> Note: On many systems, §7 (export→import) is simpler than using Windows tools to shrink VHDX.

---

## Troubleshooting

- **`Wsl/ERROR_PATH_NOT_FOUND`** when exporting:
  - The destination directory doesn’t exist. Create it:
    ```powershell
    New-Item -ItemType Directory -Path $dest -Force | Out-Null
    Test-Path (Split-Path -Parent $tar)   # should be True
    ```
- **Wrong distro name**: Use the exact NAME from `wsl -l -v` (quotes help).
- **What gets exported?** Everything inside the Linux filesystem tree. Windows pagefile and WSL’s external `swap.vhdx` are **not** included. A swap file **inside** Linux (e.g., `/swapfile`) **is** included.
- **Windows still low on space after cleanup?** Look for large `.tar/.zip` backups, OneDrive/DriveFS caches, Windows Update/Component Store, and Shadow Copies. Export→unregister→import (§7) ensures your WSL VHDX is right‑sized; the rest is Windows data.

---

## Handy one‑liner (minimal backup)

```powershell
$distro='Ubuntu'; $dest="$env:USERPROFILE\Downloads\wsl_backups"; New-Item -ItemType Directory -Path $dest -Force|Out-Null; wsl --shutdown; $stamp=Get-Date -Format yyyy-MM-dd_HH-mm; $tar=Join-Path $dest "WSL-$($distro)-$stamp.tar"; wsl --export "$distro" "$tar"; tar -tf "$tar" > $null; Get-FileHash "$tar" -Algorithm SHA256
```

---

**You’re set.** Back up periodically (e.g., monthly). When space gets weird after heavy writes, use §7 to rebuild from your tar and instantly reclaim the Windows disk space associated with VHDX growth.
