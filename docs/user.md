# 🦬 MalariAPI User Guide (Simple Install)

If you just want to **run MAPI** — no git setup, no developer steps — this is all you need.

---

## 🚀 Quick install

```bash
mkdir -p ~/MalariAPI
cd ~/MalariAPI

curl -fsSL https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/tools/install_mapi.sh -o install_mapi.sh
chmod +x install_mapi.sh

# Run with your local conda path
MINICONDA_HOME="$HOME/tools/miniconda3" ./install_mapi.sh
```

This automatically:
- Downloads the latest **MalariAPI** release from GitHub
- Installs it under `~/MalariAPI`
- Makes all modules and scripts executable
- Adds MAPI commands to your PATH
- Builds any required Conda environments

---

## ✅ Verify your setup

```bash
source ~/.bashrc      # or open a new terminal
mapi help
mapi list modules
```

You should see available MAPI modules like:

```
summarize_bam
easteregg
```

and can run one, e.g.:

```bash
mapi run modules summarize_bam --in sample.bam --out summary.txt
```

---

## 📦 Updating MAPI

You can re-run the installer any time to refresh:

```bash
cd ~/MalariAPI
./tools/install_mapi.sh
```

This will:
- Pull the latest version from GitHub
- Rebuild environments if needed
- Keep your PATH settings intact

---

---

## 🔋 Tips

- Run `mapi list modules` to see available tools
- Run `mapi help <module>` for usage info
- You can freely delete `~/MalariAPI` and reinstall — everything is reproducible
- Keep your conda installation (usually `~/tools/miniconda3`) separate from MAPI

---

## 💬 Need help?

Open an **Issue** on GitHub or contact a maintainer.  
All discussion and feature requests are welcome!

---

(c) 2025 Noah Brown and collaborators.

