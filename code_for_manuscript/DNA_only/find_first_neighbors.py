####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes find first neighbors of seed genes
####################################################################################

import os
import numpy as np
import pandas as pd

### build networkx graph from adjacency matrix
import PGM_PY.adj_mat_interactome as ami
import PGM_PY.sig_first_neighbors as sfn

####################################################################################

def find_fn(save_dir,
            to_remove,
            t_types):

    """
    Find first neighbors for each cancer type by permutation test
    """
    gm = ami.adj_mat()

    for t_type in t_types:
        save_dir1 = os.path.join(save_dir,t_type)
        fnm = os.path.join(save_dir1,'TCGA_'+t_type+'_wes_freq.csv')
        df = pd.read_csv(fnm)
        df = df.loc[~df.Gene.isin(to_remove),]
        genelist = df.Gene.values.tolist()
        df_neighbor_p = sfn.find_fn(genelist,gm)

        fnm1 = os.path.join(save_dir1,'first_neighbors_and_p-value_wes_022123.csv')
        df_neighbor_p.to_csv(fnm1,index = False)
