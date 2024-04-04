####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes find first neighbor genes
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
### build networkx graph from adjacency matrix
import PGM_PY.adj_mat_interactome as ami
import PGM_PY.sig_first_neighbors as sfn

####################################################################################
def find_first_neighbors(save_dir,t_types,to_remove):
    
    gm = ami.adj_mat()
    # save directory
    for t_type in t_types:
        print(t_type)
        file_dir = os.path.join(save_dir,t_type,'elastic net model')
        fnm = os.path.join(file_dir,'wes+rna_classifier_022123.csv')
        df = pd.read_csv(fnm)
        df = df.loc[~df.Gene.isin(to_remove),]
        genelist = df.Gene.values.tolist()
        df_neighbor_p = sfn.find_fn(genelist,gm)

        fnm1 = os.path.join(file_dir,'first_neighbors_and_p-value_wes+rna_022123.csv')
        df_neighbor_p.to_csv(fnm1,index = False)
