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
This csv file must include two columns names 'Gene' and 'Freq'. The 'Gene' column should have a list of genes with HGNC symbols, and the 'Freq' column should be a list of numeric numbers between 0 and 1 representing the altered freq of a gene within the populations. An example file can be found and downloaded at [example input file](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/input/TCGA_BRCA_DNA_altered_freq.csv)



