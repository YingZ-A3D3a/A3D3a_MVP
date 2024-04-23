import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
### MVP functions
from adaMVP import adj_mat_interactome as ami
from adaMVP import sig_first_neighbors as sfn
from adaMVP import graphical_models as gms
from adaMVP import markov_model_dataframe as mmd # degree of genes in the graph

def find_fn_and_pgm(altered_freq_file,
            save_directory = '.',
            to_remove = [],
            fn_num = 550,thre = 0.05,Wm = 0.5,alpha = 0.1,n_perm = 10000):
    """
    Find first neighbors for each cancer type by permutation test
    """
    df = pd.read_csv(altered_freq_file)
    df = df.loc[~df.Gene.isin(to_remove),]
    genelist = df.Gene.values.tolist()
    gm = ami.adj_mat()
    fn = sfn.find_fn(genelist,gm,n_perm)
    fnm1 = os.path.join(save_directory,f"first_neighbors_{time.strftime('%b%d-%H-%M')}.csv")
    fn.to_csv(fnm1,index = False)

    g_score = df.copy()
    g_score = g_score.loc[~g_score.Gene.isin(to_remove),]
    g_seed = g_score.Gene.values.tolist()
    g_score.index = g_score['Gene']
    genex = g_seed[0]
    # first neighbors
    fn = fn.loc[fn.fdr_bh<thre,]
    fn = fn.sort_values(by = ['fdr_bh','neighbors_in_seed_divide_by_seed_size'],ascending = [True, False])
    fn0 = fn.iloc[:fn_num,]
    cb = set(g_seed).union(fn0.candidate_gene.values) # combine seed with first neighbors
    ### save the gene list, specify seed or not
    in_seed = []
    cb = list(cb)
    for j in cb:
        if j in g_seed:
            in_seed.append(1)
        else:
            in_seed.append(0)
    g_df = pd.DataFrame({'gene':cb,'in_seed':in_seed})
    ### load the genelist
    tmp = g_df.loc[g_df.in_seed==0,'gene'].values
    g_all = g_df['gene'].values
    g_filter = set(g_all).intersection(gm.index)

    # get the number of altered patients of the seed list (patient count that altered)
    score_ini = {}
    for j in g_filter:
        if j in g_score.Gene.values:
            score_ini[j] = g_score.loc[j,'Freq']
    s_list = g_filter

    ### run the pgm model
    print('---------------------------------------------')
    print(f'fn:{fn_num},Wm:{Wm},alpha:{alpha}')
    final_prob_markov0 = gms.run_pipeline_unequal(gm,genex,s_list,score_ini,alpha,Wm,modelx='Markov')
    final_prob_markov = mmd.info_markov(final_prob_markov0,s_list,gm)
    source = []
    for i in final_prob_markov.genes.values:
        if i not in g_score.index:
            source.append('first neighbor')
        else:
            source.append('seed')
    final_prob_markov['source'] = source        
    ## save final rank and probability
    filenm = os.path.join(save_directory,f"markov_output_Wm_{str(Wm)}_alpha_{str(alpha)}_{time.strftime('%b%d-%H-%M')}.csv")
    final_prob_markov.to_csv(filenm, index = False)
