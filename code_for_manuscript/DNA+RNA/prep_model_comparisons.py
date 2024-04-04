####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes run Markov Chain, Pagerank, Personalized Pagerank models with different parameters
# prepare for model comparison
####################################################################################

import os
import numpy as np
import pandas as pd
import PGM_PY.graph_builder as gb # largest connected graph
import PGM_PY.graphical_models as gms
import PGM_PY.markov_model_dataframe as mmd
import PGM_PY.adj_mat_interactome as ami ### load adjacency matrix

####################################################################################
def prep_model_comparisons(save_dir,t_types,fn_min, fn_max, fn_step,
                           alpha_min, alpha_max, alpha_step):

    # ### load Markov model parameters with best performance
    # dir_mat = os.path.join(dir0,'All_Cancer_types','WES+RNA','results','output_data_matrix')
    # fnm = os.path.join(dir_mat,'rank_model_parameters_all_cancers.csv')
    # df_rnk = pd.read_csv(fnm, index_col = 0)

    ### for Markov model with different fn, save nodes and edges
    # save nodes and edges for the largest connected graph
    # run pagerank and personalized pagerank models
    
    
    fn_nums = np.arange(fn_min,fn_max+fn_step,fn_step)
    models = ['ppr','pr']
    alpha1s = np.round(np.arange(alpha_min,alpha_max,alpha_step),3)
    gm = ami.adj_mat()
    
    for t_type in t_types:
        print(t_type)
        print('--------------------------------------------')
        for fn_num in fn_nums:
            
            md_dir = os.path.join(save_dir, t_type, 'markov chain models')
            select_md = 'fn_'+str(fn_num)+'_Wm_0.3_alpha_0'
            
            s_md = 'markov_output_wes+rna_'+select_md+'_022123.csv'
            df = pd.read_csv(os.path.join(md_dir,s_md))

            g_list = df.genes.values.tolist()
            graph_connect = gb.seed_network(g_list, gm)

            # get the locations of nonzero elements of the adjacency matrix
            tmp = np.array(graph_connect)
            geneA = []
            geneB = []
            loca = np.nonzero(tmp>0)

            for i in range(len(loca[0])):
                geneA.append(graph_connect.index[loca[0][i]])
                geneB.append(graph_connect.index[loca[1][i]])
            edge_df = pd.DataFrame({'geneA':geneA, 'geneB':geneB})
            edge_df.head()

            ### save the node (in the largest connected graph) and the altered freq as a frame
            nodes = [i for i in graph_connect.index]
            ini_scores = []
            score_ini = {}

            for gene1 in nodes:
                if gene1 in df.genes.values:
                    ini_scores.append(df.loc[df.genes==gene1,'initial_score'].values[0])
                    score_ini[gene1] = df.loc[df.genes==gene1,'initial_score'].values[0]
                else:
                    ini_scores.append(0)

            node_df = pd.DataFrame({'node':nodes,'ini_score':ini_scores})
            node_df.sort_values(by = 'ini_score', ascending = False, inplace = True)

            om_dir = os.path.join(save_dir, t_type, 'other models detailed')
            os.makedirs(om_dir, exist_ok = True)

            edge_fnm = os.path.join(om_dir,'fn_'+str(fn_num)+'_largest_connected_graph_edges.csv')
            edge_df.to_csv(edge_fnm, index = False)

            node_fnm = os.path.join(om_dir,'fn_'+str(fn_num)+'_largest_connected_graph_nodes.csv')
            node_df.to_csv(node_fnm, index = False)

            ### pagerank and personalized pagerank models
            
            genex = df.genes.values[0]
            for model in models:
                for alpha1 in alpha1s:
                    final_ppr0 = gms.run_pipeline(gm,genex,g_list,score_ini,alpha0 = alpha1, modelx=model)
                    final_ppr = mmd.info_markov(final_ppr0,g_list,gm)
                    source = []

                    for item in final_ppr.genes.values:
                        source.append(df.loc[df.genes==item,'source'].values[0])
                    final_ppr['source'] = source        

                    # save final rank and probability
                    filenm = os.path.join(om_dir,'pgm_'+model+'_wes+rna_fn_'+str(fn_num)+'_alpha_'+str(alpha1)+'_022123.csv')
                    final_ppr.to_csv(filenm, index = False)

    
