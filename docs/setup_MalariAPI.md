# Configure MalariAPI Environment in mobaxterm

### 🔧 Update and Install Useful Tools

Inside WSL:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential git curl unzip htop openssh-server micro samtools bedtools
mkdir tools genomes bin envs scratch
mkdir /tools/miniconda3

#setup miniconda
cd tools/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O .
bash ../miniconda.sh -b -u -p .
rm ../miniconda.sh
source ~/tools/miniconda3/bin/activate
conda init --all

conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

cd ~/envs
curl -L -o https://github.com/A-Crow-Nowhere/MalariAPI/blob/f44a9ebb6c712a8482e1dbdfcebc19fc0c2db217/yaml/base.yml


#For a variety of dependency issues samtools and bedtools will be installed in the highest level
#executable bin. 
#I will write any text editing commands with 'micro' but you can use vim or nano if you prefer. 
```
Basic framework of MalariAPI is setup now, see how to install tools in the next totorial
