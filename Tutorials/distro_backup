# WSL Ubuntu Backup & Verification Guide (PowerShell)

> Export your entire WSL Ubuntu distro to a `.tar` without touching your current install. Works on Windows 10/11 with WSL2.

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

## 6) (Optional) Strongest check — test import safely

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

## 7) (Optional) Reclaim Windows disk space (compact VHDX)

> Not required for export. Shrinks the underlying `ext4.vhdx`. Do this occasionally.

**Inside WSL (zero free space to aid compaction):**
```bash
# Run inside the Linux shell:
sudo dd if=/dev/zero of=/zero bs=1M || true
sync
sudo rm -f /zero
```

**Back in PowerShell:**
```powershell
wsl --shutdown

# Newer WSL builds:
wsl --manage "$distro" --compact

# Or with Hyper-V tools (adjust the path):
# $vhd = "$env:LOCALAPPDATA\Packages\<YourUbuntuPackage>\LocalState\ext4.vhdx"
# Optimize-VHD -Path $vhd -Mode Full
```

---

## Troubleshooting

- **`Wsl/ERROR_PATH_NOT_FOUND`** → The destination directory doesn’t exist. Create it first:
  ```powershell
  New-Item -ItemType Directory -Path $dest -Force | Out-Null
  Test-Path (Split-Path -Parent $tar)   # should be True
  ```
- **Wrong distro name** → Use the exact NAME from `wsl -l -v` (quotes help):
  ```powershell
  wsl --export "Ubuntu-22.04" "C:\backups\WSL-Ubuntu-22.04.tar"
  ```
- **What gets exported?** Everything inside the Linux filesystem tree. Windows pagefile and WSL’s external `swap.vhdx` are **not** included. A swap file **inside** Linux (e.g., `/swapfile`) **is** included.
- **Prefer local disk for the first export**, then move the file to external/OneDrive after verifying.

---

## Handy one-liner (minimal)

```powershell
$distro='Ubuntu'; $dest="$env:USERPROFILE\Downloads\wsl_backups"; New-Item -ItemType Directory -Path $dest -Force|Out-Null; wsl --shutdown; $stamp=Get-Date -Format yyyy-MM-dd_HH-mm; $tar=Join-Path $dest "WSL-$($distro)-$stamp.tar"; wsl --export "$distro" "$tar"; tar -tf "$tar" > $null; Get-FileHash "$tar" -Algorithm SHA256
```

---

**Done.** Run the export periodically (e.g., monthly) and keep a few dated copies.
