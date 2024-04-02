import numpy as np
import pandas as pd
import networkx as nx
### function from adaMVP
from adaMVP import study_bias as study_bias
from adaMVP import graph_builder as gb

## add degree information to the output of markov model, return a data frame
def info_markov(final_prob_markov,s_list,gm):
    counts_df = study_bias.publication_count(gm)
    graph_connect = gb.seed_network(s_list,gm)
    g1 = gb.build_graph(graph_connect)
    temp = nx.degree_centrality(g1)
    temp = {k: v for k,v in sorted(temp.items(), key=lambda i:i[1], reverse = True)}
    degree_centrality_df = pd.DataFrame({'gene':[i for i in temp.keys()], 'degree_in_disease':[i for i in temp.values()]})

    ### normalized by background degree
    N = len(gm.index)
    degree_bg = []
    for j in degree_centrality_df['gene'].values:
        tmp1 = gm.index[gm.loc[j,:]>0].tolist()
        degree_bg.append(len(tmp1)/N)
    degree_centrality_df['degree_in_background'] = degree_bg
    degree_norm = np.round(np.array(degree_centrality_df['degree_in_disease'])/np.array(degree_bg),5)
    degree_centrality_df['degree_in_disease_normalized'] = degree_norm
    degree_centrality_df.head()

    deg_centrality = degree_centrality_df.copy()
    deg_centrality.index = deg_centrality['gene']

    temp = final_prob_markov['genes'].values.tolist()
    genes_l = []
    for j in temp:
        tmp1 = j.split('|')
        if len(tmp1)>1:
            genes_l.append(tmp1[0])
            continue
        genes_l.append(j)  
    deg_in_bg = deg_centrality.loc[final_prob_markov['genes'].values.tolist(),'degree_in_background'].values.tolist()
    deg_in_dis = deg_centrality.loc[final_prob_markov['genes'].values.tolist(),'degree_in_disease'].values.tolist()
    deg_in_dis_n = deg_centrality.loc[final_prob_markov['genes'].values.tolist(),'degree_in_disease_normalized'].values.tolist()
    publication = [counts_df.loc[i,'count'] if i in counts_df.index else np.nan for i in genes_l]
    final_prob_markov['degree_in_disease'] = deg_in_dis
    final_prob_markov['degree_in_background'] = deg_in_bg
    final_prob_markov['degree_in_disease_normalized'] = deg_in_dis_n
    final_prob_markov['publications'] = publication
    final_prob_markov['final_rank'] = final_prob_markov.index+1
    final_prob_markov.index = final_prob_markov.index+1
    return final_prob_markov