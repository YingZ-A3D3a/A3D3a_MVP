### study bias
### input adjacency matrix of a graph
import os
import pandas as pd

def publication_count(gm):
    # dir2 = os.path.join('PGM_PY/other_files','Study Bias')
    # current_directory = os.getcwd()
    # parent_directory = os.path.abspath(os.path.join(current_directory, os.pardir))
    # dir2 = os.path.join(parent_directory,'adaMVP','other_files')
    # filename = os.path.join(dir2,'publications_all_genes.txt')
    filename = os.path.join(os.path.dirname(__file__), 'data', 'publications_all_genes.txt')
    counts_df = pd.read_csv(filename, sep = '\t')
    # remove duplicates
    counts_df = counts_df.drop_duplicates(subset = ['gene_symbol'])
    counts_df.index = counts_df.gene_symbol

    # print(counts_df.shape)
    genes = gm.index
    genes = [i for i in genes if i in counts_df.index]    # genes in both dataset
    # print(len(genes))
    return counts_df