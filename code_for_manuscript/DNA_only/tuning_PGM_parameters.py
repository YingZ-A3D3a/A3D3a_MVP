####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes tune the parameters of Markov chain models
####################################################################################

import os
import numpy as np
import pandas as pd
import PGM_PY.adj_mat_interactome as ami # adjacency matrix
import PGM_PY.graphical_models as gms # run graphical model
import PGM_PY.graph_builder as gb # largest connected graph
import PGM_PY.markov_model_dataframe as mmd # degree of genes in the graph

####################################################################################

def tuning_pgm(save_dir,
               t_types,to_remove,thre,fn_set = np.arange(50,650,100)):
    """
    Tuning the parameters by grid search
    'thre' is threshold for the FDR of the first neighbors
    """
    gm = ami.adj_mat()

    ###### parameters ######

    ### try different pgm models
    Wms = [0.3,0.4,0.5]
    alphas = [0,0.05,0.1]
    pgm_all_cancers(fn_set,t_types,thre,Wms,alphas,gm,save_dir,to_remove)

    ### try different pgm models
    Wms = [0.6,0.7]
    alphas = [0,0.05]
    pgm_all_cancers(fn_set,t_types,thre,Wms,alphas,gm,save_dir,to_remove)

    ### try different pgm models
    Wms = [0.8]
    alphas = [0,0.02]
    pgm_all_cancers(fn_set,t_types,thre,Wms,alphas,gm,save_dir,to_remove)

    ### try different pgm models
    Wms = [0.9]
    alphas = [0,0.01]
    pgm_all_cancers(fn_set,t_types,thre,Wms,alphas,gm,save_dir,to_remove)

### run pgm model for each cancer type
def pgm_all_cancers(fn_set,t_types,thre,Wms,alphas,gm,save_dir,to_remove):
    for t_type in t_types:
        # seed genes
        save_dir1 = os.path.join(save_dir,t_type)
        save_dir2 = os.path.join(save_dir1,'markov chain models')
        os.makedirs(save_dir2, exist_ok = True)

        seed_score = os.path.join(save_dir1,'TCGA_'+t_type+'_wes_freq.csv') # seed genes initial score
        g_score = pd.read_csv(seed_score)
        g_score = g_score.loc[~g_score.Gene.isin(to_remove),]
        g_seed = g_score.Gene.values.tolist()
        g_score.index = g_score['Gene']
        genex = g_seed[0]

        # first neighbors
        fn_dir = os.path.join(save_dir1,'first_neighbors_and_p-value_wes_022123.csv') # first neighbors
        fn = pd.read_csv(fn_dir)
        fn = fn.loc[fn.fdr_bh<thre,]
        fn = fn.sort_values(by = ['fdr_bh','neighbors_in_seed_divide_by_seed_size'],ascending = [True, False])

        for fn_num in fn_set: # try different number of first neighbors
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
            graph_connect = gb.seed_network(s_list,gm)

            ### run the pgm model
            for Wm in Wms:
                for alpha in alphas:
                    print('---------------------------------------------')
                    print(f'{t_type},fn:{fn_num},Wm:{Wm},alpha:{alpha}')
                    final_prob_markov0 = gms.run_pipeline_unequal(gm,genex,s_list,score_ini,alpha,Wm,modelx='Markov')
                    final_prob_markov = mmd.info_markov(final_prob_markov0,s_list,gm)
                    source = []

                    for i in final_prob_markov.genes.values:
                        if i not in g_score.index:
                            source.append('first neighbor')
                        else:
                            source.append('WES/SV')
                    final_prob_markov['source'] = source        

                    ## save final rank and probability
                    filenm = os.path.join(save_dir2,'markov_output_wes_fn_'+str(fn_num)+'_Wm_'+str(Wm)+'_alpha_'+str(alpha)+'_022123.csv')
                    final_prob_markov.to_csv(filenm, index = False)
