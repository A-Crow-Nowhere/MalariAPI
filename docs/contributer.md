# 🤝 Contributing to MalariAPI

Thank you for contributing to **MalariAPI (MAPI)**!  
This document explains how to set up your environment, create new modules or pipelines, and submit changes safely.

---

## 🩉 Repo structure overview

```
MalariAPI/
├── bin/              → MAPI launcher + helper wrappers
├── modules/          → Single-step MAPI tools (e.g., summarize_bam.sh)
├── pipeline/         → Multi-step pipelines
├── packages/         → External R/Python integrations
├── scripts/          → Standalone conversion scripts
├── templates/        → Boilerplate YAMLs and bash templates
├── tools/            → Installer, validator, environment YAMLs
│   ├── yaml/         → Environment specs (installer builds from these)
│   └── install_mapi.sh
└── envs/             → Actual Conda runtime environments (auto-created)
```

---

## ⚙️ Developer setup (first time)

1. **Fork** the repository on GitHub:
   - Go to [A-Crow-Nowhere/MalariAPI](https://github.com/A-Crow-Nowhere/MalariAPI)
   - Click **Fork** → create your own copy under your GitHub account.

2. **Set up SSH authentication** (if not already done):

   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519
   cat ~/.ssh/id_ed25519.pub
   ```
   Then:
   - Add your public key to GitHub → **Settings → SSH and GPG keys → New SSH key**
   - Test it:
     ```bash
     ssh -T git@github.com
     ```

3. **Clone your fork and install dependencies:**

   ```bash
   git clone git@github.com:<your-username>/MalariAPI.git
   cd MalariAPI
   MINICONDA_HOME="$HOME/tools/miniconda3" ./tools/install_mapi.sh
   ```

---

## 🧠 Development workflow

### 1. Create a new feature branch
```bash
git checkout -b feature/<short-description>
# examples:
#   feature/new-module
#   fix/yaml-paths
#   docs/typo
```

### 2. Edit or add modules/pipelines
Use the templates under `templates/` to ensure validator compatibility.

```bash
cp templates/module_template.sh modules/my_new_tool.sh
cp templates/module_template.yml modules/yaml/my_new_tool.yml
```

Follow the in-file comments for structure and metadata.

### 3. Test locally
```bash
mapi list modules
mapi run modules my_new_tool --help
```

If your module needs a new environment, place a YAML under `tools/yaml/` and rerun the installer.

### 4. Commit and push
```bash
git add -A
git commit -m "Add: new module my_new_tool"
git push -u origin feature/my_new_tool
```

### 5. Open a Pull Request (PR)
- Go to your fork on GitHub → **Compare & pull request**
- Target: `main` branch of `A-Crow-Nowhere/MalariAPI`
- Include:
  - Purpose / context
  - Changes made
  - Example usage

A maintainer will review and merge once checks pass.

---

## 🔐 Branch and merge policy

| Branch | Purpose | Who can push | Merge method |
|--------|----------|--------------|---------------|
| `main` | Stable, release-ready code | Maintainers only | Squash or rebase PRs |
| `feature/*` | Active development | Developers | via PR |
| `fix/*` | Bug fixes | Developers | via PR |
| `docs/*` | Documentation changes | Developers | via PR |

Protected branches (`main`) **cannot** be pushed to directly except by maintainers.

---

## 🧪 Testing & validation

All new modules/pipelines should:
- Pass syntax checks (`bash -n`, YAML linter)
- Run without external dependencies outside defined environments
- Validate cleanly using the MAPI validator:
  ```bash
  ./tools/validate.sh bin/modules/my_new_tool.sh
  ```

Optional CI workflows can run these checks automatically (future feature).

---

## 🪪 Code of conduct
Be kind, clear, and constructive.  
We are building a shared ecosystem for reproducible malaria genomics — treat every contribution as part of that collective mission.

---

## 🏷️ Attribution

Contributors will be credited in the release notes and documentation.  
If your module is published or referenced in a paper, please add appropriate citation info to its YAML metadata.

---
