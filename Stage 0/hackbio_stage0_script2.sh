#part 2 installation script
#!/bin/bash

# Step 1: Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh

# Step 2: Install Miniconda silently into ~/miniconda
bash miniconda.sh -b -p $HOME/miniconda

# Step 3: Initialize conda in this script
source $HOME/miniconda/etc/profile.d/conda.sh

# Step 4: Check conda version
conda --version

# Step 5: Create a conda environment named 'funtools'
conda create -n funtools -y

# Step 6: Activate funtools environment
conda activate funtools

# Step 7: Install figlet from conda-forge
conda install -c conda-forge figlet -y

# Step 8: Run figlet with your name
figlet "Mariam"

# Step 9: Add bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Step 10: Install bioinformatics tools
conda install -c bioconda bwa blast samtools bedtools spades bcftools fastp multiqc -y

-------------------


#Link to LinkedIn Video:
https://www.linkedin.com/feed/update/urn:li:activity:7367566727838171136/
