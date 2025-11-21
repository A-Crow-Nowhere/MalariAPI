# WSL Ubuntu Backup, Verify, and Restore Guide (PowerShell + WSL)

> Safely back up your entire WSL Ubuntu distro to a `.tar`, verify it, and (optionally) **rebuild/restore** from the tar to clean up virtual disk bloat — no Hyper‑V required. Copy‑paste friendly.

---

## 0) Notes & Terminology
- **Export (`wsl --export`)** → packs only your Linux files into a `.tar` (no free space).
- **Import (`wsl --import`)** → creates a fresh distro/VHDX from a `.tar` (right‑sized).
- **Recommended cleanup** after heavy writes: **export → unregister → import** (see §7).

Adjust `Ubuntu` to your exact name from `wsl -l -v`.

---

## 1) Find your exact distro name
```powershell
wsl -l -v
$distro = 'Ubuntu'   # replace with the NAME shown above
```

---

## 2) Choose a backup folder
```powershell
$dest = "$env:USERPROFILE\Downloads\wsl_backups"
New-Item -ItemType Directory -Path $dest -Force | Out-Null
```

---

## 3) Export (backup)
```powershell
wsl --shutdown
$stamp = Get-Date -Format yyyy-MM-dd_HH-mm
$tar   = Join-Path $dest "WSL-$($distro)-$stamp.tar"
wsl --export "$distro" "$tar"
Get-Item "$tar" | Format-List Name,Length,LastWriteTime,FullName
```

---

## 4) Quick validation
```powershell
tar -tf "$tar" > $null; $LASTEXITCODE   # 0 = OK
tar -xOf "$tar" etc/os-release | Select-Object -First 12
Get-FileHash "$tar" -Algorithm SHA256 | Tee-Object -FilePath "$tar.sha256.txt"
```

---

## 5) (Optional) Smaller backup
```powershell
# After export
Compress-Archive -Path "$tar" -DestinationPath "$tar.zip"
# OR stream-compress directly:
wsl --shutdown
$gz = Join-Path $dest "WSL-$($distro)-$stamp.tar.gz"
wsl --export "$distro" - | wsl -e gzip -c > "$gz"
wsl -e gzip -t "$gz"; $LASTEXITCODE     # 0 = OK
```

---

## 6) Strongest check: test import safely
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

## 7) **Recommended cleanup**: rebuild your distro (right‑size VHDX)
```powershell
# Assumes $distro, $dest set above
wsl --shutdown
$stamp = Get-Date -Format yyyy-MM-dd_HH-mm
$tar   = Join-Path $dest "WSL-$($distro)-$stamp.tar"
wsl --export "$distro" "$tar"

# Remove old instance (frees bloated VHDX)
wsl --unregister "$distro"

# Import from tar (creates fresh, compact VHDX)
$newHome = "$env:USERPROFILE\WSL\$distro"
New-Item -ItemType Directory -Path $newHome -Force | Out-Null
wsl --import "$distro" "$newHome" "$tar" --version 2

# Sanity check
wsl -d "$distro" -- uname -a
```

---

## 8) Set default user and HOME after import
**Inside WSL (Ubuntu):**
```bash
# Replace 'youruser' with your preferred username
USER=youruser
sudo adduser "$USER" || true
sudo usermod -aG sudo "$USER" || true
sudo mkdir -p /home/$USER
sudo usermod -d /home/$USER -m "$USER" 2>/dev/null || true

# Make WSL start as this user by default
printf "[user]\ndefault=%s\n" "$USER" | sudo tee /etc/wsl.conf
exit
```
**Back in PowerShell:**
```powershell
wsl --shutdown
wsl -d "$distro"   # opens as your user in /home/<user>
```

**Optional (Windows Terminal starting directory):** set to
```
\\wsl.localhost\<YourDistro>\home\<youruser>
```
or in `settings.json`:
```json
"startingDirectory": "\\\\wsl.localhost\\<YourDistro>\\home\\<youruser>"
```

---

## 9) **Shareable tar (beginner‑friendly)** — no accidental username carry‑over

### Make your tar generic (on your machine, before exporting)
**Inside WSL:** install a tiny first‑run helper and set default to root.
```bash
# /root/init-user.sh
sudo tee /root/init-user.sh >/dev/null <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
read -rp "New Linux username: " U
adduser "$U"
usermod -aG sudo "$U"
printf "[user]\ndefault=%s\n" "$U" > /etc/wsl.conf
echo "User $U created and set as default. Now run: wsl --shutdown, then reopen."
EOF
sudo chmod +x /root/init-user.sh

# Make root the default for the exported image
printf "[user]\ndefault=root\n" | sudo tee /etc/wsl.conf
```
**Back in PowerShell:** export your generic tar.
```powershell
wsl --shutdown
$stamp = Get-Date -Format yyyy-MM-dd_HH-mm
$tar   = Join-Path $dest "WSL-$($distro)-$stamp-generic.tar"
wsl --export "$distro" "$tar"
```

### Recipient’s simple restore (on their PC)
```powershell
# Import (creates a new distro from your tar)
wsl --import Ubuntu C:\WSL\Ubuntu C:\path\to\WSL-Ubuntu-generic.tar --version 2

# Guided user creation (prompts for username)
wsl -d Ubuntu -u root -- /root/init-user.sh

# Apply and reopen as the new user
wsl --shutdown
wsl -d Ubuntu
```
This flow avoids teaching user management; they just run one helper script.

---

## 10) (Optional) Zero‑fill free space inside Linux (why/when)
Helps compaction/imaging by turning free blocks into zeros. Prefer §7 (rebuild) for cleanup.
```bash
sudo dd if=/dev/zero of=/zero bs=1M status=progress || true
sync && sudo rm -f /zero
sudo fstrim -av || true
```

---

## Troubleshooting
- **`Wsl/ERROR_PATH_NOT_FOUND`** when exporting: ensure the destination folder exists.
- **Wrong distro name**: must match `wsl -l -v` exactly (quotes help).
- **What gets exported?** Files inside Linux only; Windows pagefile and WSL’s external `swap.vhdx` are not included.

---

## Minimal one‑liner (backup)
```powershell
$distro='Ubuntu'; $dest="$env:USERPROFILE\Downloads\wsl_backups"; New-Item -ItemType Directory -Path $dest -Force|Out-Null; wsl --shutdown; $stamp=Get-Date -Format yyyy-MM-dd_HH-mm; $tar=Join-Path $dest "WSL-$($distro)-$stamp.tar"; wsl --export "$distro" "$tar"; tar -tf "$tar" > $null; Get-FileHash "$tar" -Algorithm SHA256
```

---

**Done.** Back up monthly; use §7 (rebuild) any time disk usage gets weird after heavy writes.
