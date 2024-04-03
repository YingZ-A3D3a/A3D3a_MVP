[![Python package](https://img.shields.io/pypi/v/adaMVP.svg?color=brightgreen&label=python-package)](https://pypi.org/project/adaMVP)

Documentation
============================================
A3D3a’s MVP (Adaptive AI-Augmented Drug Discovery and Development Molecular Vulnerability Picker) is a novel graph-based, cooperativity-led Markov chain model. 

## Output files
After inputing a seed genes file and running the code for finding first neighbors and building graphical model by MVP, the model will output two files:

### 1. first_neighbors.csv
example like [first_neighbors.csv](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/docs/example_output/first_neighbors_Apr03-11-14.csv)

It contains the candidate first neighbor genes of the seed gene list, and their likelihood of being first neighbors is estimated by permutation test.
   
  It has columns:
- `candidate_gene`: candidate first neighbor gene of the seed gene list
- `neighbors`: number of first neighbors of the candidate gene in the background interactome
- `neighbors_in_seed`: number of first neighbors of the candidate gene in the seed gene list
- `neighbors_divide_by_interactome_size`: number of first neighbors of the candidate gene in the background interactome divided by background interactome size
- `neighbors_in_seed_divide_by_seed_size`: number of first neighbors of the candidate gene in the seed gene list divided by the size of seed gene list
- `ratios`: the ratio of `neighbors_in_seed_divide_by_seed_size` to `neighbors_divide_by_interactome_size`
- `p_value`: under permutation test, how often the observed ratio could be achieved by chance. The permutation test is performed by permutating the interactome to generate new seed lists of the same size as the original seed list
- `fdr_bh`: p-values are corrected by multiple hypothesis testing using the Benjamini and Hochberg method
### 2. markov_output.csv
example like [markov_output.csv](https://github.com/YingZ-A3D3a/A3D3a_MVP/blob/main/docs/example_output/markov_output_Wm_0.5_alpha_0.1_Apr03-11-14.csv)

It contains the disease network which is composed of the seed genes and the first neighbors. It provides the rank of the importance of the disease genes based on the MVP model prioritization.

  It has columns:
- `genes`: the genes in the disease network
- `final_score`: the final score by the MVP model for each gene
- `initial_score`: the initial score which is the altered rate/freq of each gene among the disease population
- `degree_in_disease`: number of neighbors of a gene in the disease network normalized by the size of the disease network
- `degree_in_background`: number of neighbors of a gene in the background interactome normalized by the size of the background interactome
- `degree_in_disease_normalized`: ratio of degree_in_disease to degree_in_background
- `publications`: number of publications in pubmed for each gene. (data `gene2pubmed.gz` downloaded from [ncbi](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) on August 2022)
- `final_rank`: the final rank by the MVP model for each gene
- `source`: whether the gene is 'seed' or 'first neighbor'
