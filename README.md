[![Python package](https://img.shields.io/pypi/v/adaMVP.svg?color=brightgreen&label=python-package)](https://pypi.org/project/adaMVP)

# A3D3a's MVP (adaMVP)
This is the official codebase for **adaMVP: Probabilistic graph-based model uncovers druggable vulnerabilities in major solid cancers.**

## What is A3D3a's MVP?
A3D3a’s MVP (Adaptive AI-Augmented Drug Discovery and Development Molecular Vulnerability Picker) is a novel graph-based, cooperativity-led Markov chain model, developed and maintained by [Bissan Al-lazikani lab](https://faculty.mdanderson.org/profiles/bissan_al_lazikani.html) at the University of Texas MD Anderson Cancer Center. The algorithm exploits cooperativity of weak signals within a cancer molecular network to enhance the signal of true molecular vulnerabilities. 

## Workflow of A3D3a's MVP
![workflow](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/docs/workflow.png)

## Installation

adaMVP works with Python >= 3.8. Please make sure you have the correct version of Python installed pre-installation.

We highly recommend using an isolated python environment using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html).
1. Create python>=3.8 environment
   - Using conda: `conda create -n ada_mvp python=3.8`
   - Using virtualenv: `python -m venv ada_mvp`

2. Activate environment
   - Using conda: `conda activate ada_mvp`
   - Using virtualenv: `source ada_mvp/bin/activate`

After setting the environment, you could install adaMVP via pip:

```bash
pip install adaMVP
```

## Preparing INPUTS
#### Preparing your seed genes input file (mandatory)
This csv file must include two column names 'Gene' and 'Freq'. The 'Gene' column should have a list of genes with official gene symbols (HGNC symbols), and the 'Freq' column should be a list of numeric numbers between 0 and 1 representing the altered freq of a gene within the populations. An example file can be downloaded at [example input file](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/input/TCGA_BRCA_DNA_altered_freq.csv)

## Run examples
#### Find first neighbors and build MVP model for prioritizing molecular vulerabilities
```shell
from adaMVP import mvp_build_graph as mbg

mbg.find_fn_and_pgm(save_directory = output_path,
            to_remove = ['TTN','MUC16'],
            altered_freq_file = 'TCGA_BRCA_DNA_altered_freq.csv',
            fn_num = 550,
            thre = 0.05,
            Wm = 0.5,
            alpha = 0.1,
            n_perm = 1000)
```
### Parameters (mandatory)
- `altered_freq_file`: input file with altered freq for each gene in a csv file
  
### Optional Parameters
- `save_directory`: directory path for saving output files
- `to_remove`: Genes to be filtered from the seed genes, default=[]
- `threshold`: Threshold of FDR for the permutation test for finding first neighbors of the seed genes, default=0.05
- `fn_num`: Maximum number of first neighbors to be brought into the network, sorted by the FDR and number of neighbors in the seeds, default=550
- `n_perm`: Number of iterations for the permutation test, default=10000
- `Wm`: weight parameter on the self-loop of nodes of the Markov chain model, default=0.5
- `alpha`: cooperativity factor describing transition between nodes, default=0.1

### Output files
Two output files, 'first_neighbors.csv' and 'markov_output.csv' will be saved in the assigned directory `output_path`. To understand the results, please check the [results_documentation](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/docs/results_documentation.md)

### Tutorial 
The tutorial for running the adaMVP pipeline can be found at 

[build graph and community detection](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/tutorial/graph_modeling_and_community_detection.ipynb)





