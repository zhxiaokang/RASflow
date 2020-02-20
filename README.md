# RASflow: RNA-Seq Analysis Snakemake Workflow
RNA-Seq analysis workflow using Snakemake
## Workflow
<img src="https://github.com/zhxiaokang/RNA-Seq-analysis/blob/master/workflow/workflow_chart.jpg" width="450">

## Quick start
### Installation
Clone the repository:

`git clone https://github.com/zhxiaokang/RASflow.git`

Create the environment:

`conda env create -n rasflow -f env.yaml`

Activate the environment:

`conda activate rasflow`

### Set up configuration
Modify the metafile describing your data `configs/metadata.tsv`.

Customize the workflow based on your need in `configs/config_main.yaml`.

### Run RASflow
`python main.py`

## Tutorial
A more detailed tutorial of how to use this workflow can be found here: [Tutorial](https://github.com/zhxiaokang/RASflow/blob/master/Tutorial.pdf)

## Evaluation
RASflow has been evaluated on 4 datasets including two model organisms (human and mouse) and a non-model organism (Atlantic cod). To keep this repository as light as possible, the evaluation of RASflow on real datasets is deposited here: [RASflow_realData](https://git.app.uib.no/Xiaokang.Zhang/rasflow_realdata)

## References
RASflow is available as a preprint on bioRxiv: [Zhang X, Jonassen I. RASflow: An RNA-Seq Analysis Workflow with Snakemake. bioRxiv. 2019 Nov 839191.](https://www.biorxiv.org/content/10.1101/839191v1)

It has been submitted to BMC Bioinformatics and is currently under review.
