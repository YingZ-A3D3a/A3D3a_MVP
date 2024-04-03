# A3D3a's MVP (adaMVP)
This is the official codebase for **adaMVP: Probabilistic graph-based model uncovers druggable vulnerabilities in major solid cancers.**

## What is A3D3a's MVP (adaMVP)?
A3D3a’s MVP (Adaptive AI-Augmented Drug Discovery and Development Molecular Vulnerability Picker) is a novel graph-based, cooperativity-led Markov chain model. The algorithm exploits cooperativity of weak signals within a cancer molecular network to enhance the signal of true molecular vulnerabilities. 

## Installation

adaMVP works with Python >= 3.7. Please make sure you have the correct version of Python installed pre-installation.

We highly recommend using an isolated python environment using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) or [virtualenv](https://docs.python.org/3/library/venv.html).

After setting the environment, you could install adaMVP via pip:

```bash
pip install adaMVP
```

## Preparing INPUTS
#### Preparing your seed genes input file (mandatory)
This csv file must include two columns names 'Gene' and 'Freq'. The 'Gene' column should have a list of genes with HGNC symbols, and the 'Freq' column should be a list of numeric numbers between 0 and 1 representing the altered freq of a gene within the populations. An example file can be downloaded at [example input file](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/input/TCGA_BRCA_DNA_altered_freq.csv)

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
- `save_directory`: directory path for saving output files
- `altered_freq_file`: input file with altered freq for each gene in a csv file
  
### Optional Parameters
- `to_remove`: Genes to be filtered from the seed genes, default=[]
- `threshold`: Threshold of FDR for the permutation test for finding first neighbors of the seed genes, default=0.05
- `fn_num`: Maximum number of first neighbors to be brought into the network, sorted by the FDR and number of neighbors in the seeds, default=550
- `n_perm`: Number of iterations for the permutation test, default=10000
- `Wm`: weight parameter on the self-loop of nodes of the Markov chain model, default=0.5
- `alpha`: cooperativity factor describing transition between nodes, default=0.1

### Output files
Two output files, 'first_neighbors.csv' and 'markov_output.csv' will be saved in the assigned `output_path`. To understand the results, please check the [results_documentation](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/docs/results_documentation.md)





