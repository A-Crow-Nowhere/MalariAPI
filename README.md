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

## Choose your installation path
Pick the box that matches what you did and what you want.

<p align="center">
  <img src="MalariAPI/docs/installOptions_croppted.pdf" width="800">
</p>

<p align="center">
  <a href="MalariAPI/docs/installOptions_croppted.pdf">Download PDF (copy-paste friendly)</a>
</p>

### After any of the above

```bash
source ~/.bashrc    # or: source ~/.zshrc
mapi --help
rc    # or: source ~/.zshrc
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


