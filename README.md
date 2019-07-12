# RAN-Seq-analysis workflow
RNA-Seq analysis workflow using Snakemake
## Snakemake environment setup
By setting up the environment, all the required tools will be installed with only one command. Besides of that, since the versions of all tools can be specified, the work is ensured to be reproducible.
### Install Miniconda
Download Miniconda installer from here: https://docs.conda.io/en/latest/miniconda.html Install it to your laptop or server.
### Install Snakemake
Open your favorate terminal, and type in the following command to install Snakemake via Conda
```
conda install -c bioconda -c conda-forge snakemake
```
### Download the repository
Download the repository from here: https://github.com/zhxiaokang/RNA-Seq-analysis/archive/master.zip and unzip it.
### Set up the environment for the analysis
Open your terminal in the directory of the repository RNA-Seq-analysis/ and run the following command to set up the environment:
```
conda env create --name RNA-Seq-analysis --file ./configs/environment.yaml
```
Then activate the environment:

For Linux and Mac:
```
source activate RNA-Seq-analysis
```
For Windows:
```
activate RNA-Seq-analysis
```
Now you're ready to go!
