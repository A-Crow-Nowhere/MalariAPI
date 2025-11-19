# Welcome to MalariAPI (MAPI)
_An all in one, beginner friendly, and easy to use system for compuational investigations of malaria._ 
Because the organisms that cause malaria (and releated apicomplexan parasites) are so complex, tools that are designed for use in broader contexts are not neccisarily optimzed for malaria. Here, every module and pipeline has been pre-optimized for specific parameters, allowing for a precise and starndardized workflow. MAPI includes steps for setting up a ready-to-use bioninformatics workspace; easy enough for coding beginners, but robust and flexible enough to be implemented by verteran bioinformaticians. 
## Quick Links:
   ### Set up base MAPI
   1. [Set up a functional environment](docs/setup_MalariAPI.md) 
   2. [Install base MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/PopulateMAPI.md)
   3. [Beginner user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/user.md)
   ### Set up advanced MAPI
   4. [Contributer user guide for MAPI](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/contributer.md)
   5. [Set up job subisison and sycing to a high performance computing (HPC) cluster](https://github.com/A-Crow-Nowhere/MalariAPI/blob/main/docs/HPC_Guide.md)   
   6. [Creating modules and pipelines](docs/wrapping_tools.md)
   ### Helpful scripts
   7. [A guide on how to backup, and cleanup local distros of Ubuntu](docs/distro_backup.md)

## Quick Install MalariAPI (MAPI)

Below, we outined a streamlined setup process for MalariAPI, covering both local installation and optional HPC configuration.

## 1. (Optional) Fork the Repository

If you plan to contribute:

1. Visit: https://github.com/A-Crow-Nowhere/MalariAPI
2. Select "Fork" to create your personal copy.

This lets you push to your own fork while still pulling updates from the main repository.
^ You can think of a 'fork' as your own personal copy of MAPI that you can do whatever you want with

## Installation Path Selection

This guide helps you choose the correct installation commands depending on whether you plan to develop MalariAPI or simply use it.

---

### User Installation (No Development)

Select one of the following based on how you obtained MalariAPI.

| Case | Description | Commands |
|------|-------------|----------|
| **Case 2** | Clone the main repo and just use it |```bash
git clone https://github.com/A-Crow-Nowhere/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI

./tools/install_mapi.sh
``` |
| **Case 4** | Forked MalariAPI, cloned your fork, but not developing |```bash
git clone git@github.com:<your-github-username>/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI

./tools/install_mapi.sh
``` |
| **Case 6** | Cloned someone else's MalariAPI repo and just using it |```bash
git clone git@github.com:<other-owner>/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI

./tools/install_mapi.sh
``` 

---

### Developer Installation (Planning to Modify or Contribute)

Select one of the following if you intend to write code, change modules, or submit PRs.

| Case | Description | Commands |
|------|-------------|----------|
| **Case 1** | Clone the main repo and develop |```bash
git clone https://github.com/A-Crow-Nowhere/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI

./tools/install_mapi.sh \
  --git-name "Your Name" \
  --git-email you@example.edu
``` |
| **Case 3** | Fork on GitHub, clone your fork, and develop |```bash
git clone git@github.com:<your-github-username>/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI

./tools/install_mapi.sh \
  --repo-owner "<your-github-username>" \
  --upstream-owner "A-Crow-Nowhere" \
  --git-name "Your Name" \
  --git-email you@example.edu
``` |
| **Case 5** | Clone someone else's repo and develop |```bash
git clone git@github.com:<other-owner>/MalariAPI.git ~/MalariAPI
cd ~/MalariAPI

./tools/install_mapi.sh \
  --repo-owner "<other-owner>" \
  --upstream-owner "A-Crow-Nowhere" \
  --git-name "Your Name" \
  --git-email you@example.edu
``` |

---

# Final Step for All Cases

```bash
source ~/.bashrc     # or: source ~/.zshrc
mapi --help
```

You are now ready to use or develop MalariAPI.


## 4. Install MalariAPI on Your HPC (optional)

MAPI supports the use of a High Preformance Cluster while operating out of a local environment. This reduces complexity for new users, while ensuring stable running environments. 

From your local clone:

```bash
cd ~/MalariAPI

./tools/hpc_install  <the_name_of_your_hpc (or any single word)>  <host_name> <partition> <allocation>

eg:

./tools/hpc_install  <rivanna> <userID@login.hpc.virginia.edu> <standard> <mymalarialab>

```
Test the connection:

```bash
mapi rivanna status
```

A valid response confirms successful HPC setup, should look like a single line:
```
'             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)'
```

## 5. Git Helper Commands in MAPI (optional, but recomended) 

Using the MAPI installer automatically wires you in to be able to use github. MAPI provides simple wrappers for common git operations. 

### Update your local branch from upstream (i.e. from main):

```bash
mapi git update
```

### Upload your current branch to your fork (i.e. your devloper environment):

```bash
mapi git upload
```

### Switch or create a branch (if you have access to multiple branches):

```bash
mapi git switch feature/new_module
```

## 6. Next Steps

List available modules and pipelines:

```bash
mapi 
```


Build any additional Conda environments as needed with MAPI’s environment tools.

## Installation Complete

You now have a fully functional local and HPC installation of MalariAPI.
Refer to the Contributing section for details on module and pipeline development.


