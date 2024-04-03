################ Import packages ################
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
from tqdm import tqdm

def find_fn(s_list_comb, gm, n_perm):  # input a gene list g_list and the adjacency matrix of the interactome gm
    fn = []
    # graph_seed_list = graph_seed.index.tolist()
    print(f"seed list size: {len(s_list_comb)}")
    for i in s_list_comb:
        if i in gm.index:
            tmp1 = gm.index[gm.loc[i,:]>0].tolist()
            fn.extend(tmp1)
    fn = list(set(fn))
    fn = [i for i in fn if i not in s_list_comb]
    # print(f"number of first neighbors of the seed graph: {len(fn)}")

    ### find the first neighbors of the seed_list that are not in the seed_list
    ### filter and keep the first neighbors that have at least 2 neighbors in the seed_list
    # or the only neighbor is in the seed_list

    fn_kp = []     # genes to keep after filtering
    for i in fn:
        tmp1 = gm.index[gm.loc[i,:]>0].tolist()
        overlap = list(set(tmp1).intersection(s_list_comb))
        if len(overlap)>=2 or (len(tmp1)==1 and len(overlap)==1):
            fn_kp.append(i)
    # print(f"number of first neighbors of the seed graph that passed the filter: {len(fn_kp)}")

    ### first neighbor list of candidate first neighbor genes
    fn_fn = []
    for i in fn_kp:
        tmp1 = gm.index[gm.loc[i,:]>0].tolist()
        fn_fn.append(tmp1)
        
    ### proportion of neighbors of first neighbors that are in the seed_list
    def compute_ratio_in_s_list(g_list_perm, n_all, n_list):
        ratio_in_s_list = []
        for i in range(len(fn_kp)):
            tmp1 = fn_fn[i]
            overlap = list(set(tmp1).intersection(g_list_perm))
            fn_over_all = len(tmp1)/n_all
            fn_over_s_list = len(overlap)/n_list
            temp = fn_over_s_list/fn_over_all
            ratio_in_s_list.append(temp)
        return ratio_in_s_list

    ### permutation test, method 1, enrichment score (ratio over ratio)
    n_all = len(gm.index)
    n_list = len(s_list_comb)

    ratio_in_s_list_original = compute_ratio_in_s_list(s_list_comb, n_all, n_list)
    x_ref = np.array(ratio_in_s_list_original)
    x_ref = x_ref.reshape((len(fn_kp),1))
    x_count = np.zeros((len(fn_kp),1))
    print('Start permutation test for finding first neighbors:')
    for i in tqdm(range(n_perm)):
        indx_p = np.random.permutation(len(gm.index))
        g_list_p = [gm.index[indx_p[i]] for i in range(len(s_list_comb))]

        tmp_list = compute_ratio_in_s_list(g_list_p, n_all, n_list)
        x1 = np.array(tmp_list)
        x1 = x1.reshape((len(fn_kp),1))
        x_count += x1>x_ref

    p_values = x_count/n_perm
    p_values_ = p_values.reshape(len(fn_kp),)
    adj_p = multipletests(p_values_,method = 'fdr_bh')

    ### create dataframe of the p-value
    number_of_neighbor_in_seed = []
    number_of_neighbor = []
    number_of_neighbor_in_seed_to_seed_size = []
    number_of_neighbor_to_interactome_size = []
    IoUs = []
    ratios = []
    for i in fn_kp:
        tmp1 = gm.index[gm.loc[i,:]>0].tolist()
        overlap = list(set(tmp1).intersection(s_list_comb))
        union = list(set(tmp1).union(s_list_comb))
        number_of_neighbor.append(len(tmp1))
        number_of_neighbor_in_seed.append(len(overlap))
        number_of_neighbor_to_interactome_size.append(len(tmp1)/n_all)
        number_of_neighbor_in_seed_to_seed_size.append(len(overlap)/n_list)
        ratios.append((len(overlap)/n_list)/(len(tmp1)/n_all))
        # IoUs.append(len(overlap)/len(union))
        
    df_neighbor_p = pd.DataFrame({'candidate_gene':fn_kp,
                                'neighbors':number_of_neighbor,
                                'neighbors_in_seed':number_of_neighbor_in_seed,
                                'neighbors_divide_by_interactome_size':number_of_neighbor_to_interactome_size,
                                'neighbors_in_seed_divide_by_seed_size':number_of_neighbor_in_seed_to_seed_size,
                                'ratios':ratios,
                                'p_value':p_values_,
                                'fdr_bh':adj_p[1]})
    df_neighbor_p.sort_values(by = ['fdr_bh','neighbors_in_seed_divide_by_seed_size'],ascending = [True, False], inplace = True)
    df_neighbor_p.reset_index(drop = True, inplace = True)
    # print(df_neighbor_p.head())
    return df_neighbor_p