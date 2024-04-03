import numpy as np
import pandas as pd
import networkx as nx
### function from adaMVP
from adaMVP import graph_builder as gb

######################################################################################################
### build networkx graph from adjacency matrix
def build_graph(adjacency_matrix):
    g = nx.from_pandas_adjacency(adjacency_matrix)
    # print(nx.info(g))
    return g

######################################################################################################
### graphical models
## markov model, unequal weight assumption, output the stationary distribution of all genes in the graph

def markov_model_unequal_weight(p_p1,gm_seed,Wm,alpha):
    tm = gm_seed.copy()
    p_p = np.array(p_p1)
    p_p = Wm*(p_p-min(p_p))/(max(p_p)-min(p_p))

    ### set transition matrix
    for i in range(len(tm)):
        ind = np.where(tm.iloc[i,]==1)[0].tolist()
        self_loop = p_p[ind]
        N = len(ind)
        b = (1-p_p[i]-np.sum(self_loop)*alpha)/N
        weight = b+np.array(self_loop)*alpha
        
        tm.iloc[i,i] = p_p[i]
        for j in range(len(ind)):
            tm.iloc[i,ind[j]] = weight[j]

    if (tm<0).sum().sum()>0: # if there is any negative elements
        print('there is negative elements in the transition matrix!')

    ### markov chain model
    pi = [1/len(tm) for i in tm]  # initial state
    pi = np.array(pi)
    TM = tm.to_numpy()
    
    max_steps = 10000
    
    for i in range(max_steps+1):
        temp = np.matmul(pi,TM)
        if np.array_equal(np.round(temp,8),np.round(pi,8)):
#            print(f"chain state converges at step {i}")
            break
        else:
            pi = temp
            # if i%10==0:
            #     print(f"=========step {i}: {temp[:5]} ...")
    if i==max_steps:
        print('markov chain did not converge!')
        
    ### output probability for all genes
    genes = tm.index.tolist()
    final_prob = pi.tolist()
    output = pd.DataFrame(list(zip(genes, final_prob)), columns = ['genes','final_score'])
    output['initial_score']=p_p1
    output_df = output.sort_values(by=['final_score'], ascending=False, ignore_index = True)
    return output_df

## personalized pagerank model
def ppr_model(p_p1,gm_seed,alpha0,modelx):
    p_p = {}
    genes_c = gm_seed.index.tolist()
    for i in range(len(genes_c)):
        p_p[genes_c[i]] = p_p1[i]
    g_c = build_graph(gm_seed)
    if modelx == 'ppr':
        pg_rank = nx.pagerank(g_c, alpha = alpha0, personalization = p_p)
    elif modelx == 'pr':
        pg_rank = nx.pagerank(g_c, alpha = alpha0)
    a = sorted(pg_rank.items(), key=lambda x: x[1], reverse = True)    
    genes = []
    scores = []
    ini_score = []
    for i in a:
        genes.append(i[0])
        scores.append(i[1])
        ini_score.append(p_p[i[0]])
    df_pg_rank = pd.DataFrame({'genes':genes,'predicted_score':scores,'initial_score':ini_score})
    return df_pg_rank

### Markov model pipelines
# set the weight of nodes in the seed list proportional to the mutation rate

def pipeline_uequal_markov(s_list, gm, genex, score_ini, Wm, alpha):
# graph_seed: largest connected seed network
# gm: adjacency mat of full interactome, genex: gene of interest,
# score_ini:altered rate of all genes in the seed, Wm: max weight on self loop

    graph_newnet = gb.seed_network(s_list,gm)  # adjacency matrix of disease network
    network_size = graph_newnet.shape[0]
    p_p0 = load_mutation_rate(graph_newnet, score_ini)
    final_prob = markov_model_unequal_weight(p_p0,graph_newnet,Wm,alpha)
    if genex in final_prob['genes'].values:
        n_neigh = get_genex_neighbors(graph_newnet,genex)
        rank = final_prob.index[final_prob['genes']==genex].values[0]+1
    else:
        n_neigh = np.nan
        rank = np.nan
    return [network_size,n_neigh,rank,final_prob]

def pipeline_ppr(s_list,gm,genex,score_ini,alpha0,modelx):
# graph_seed: largest connected seed network
# gm: adjacency mat of full interactome, genex: gene of interest,
# score_ini:altered rate of all genes in the seed, Wm: max weight on self loop

    graph_newnet = gb.seed_network(s_list,gm)
    network_size = graph_newnet.shape[0]
    n_neigh = get_genex_neighbors(graph_newnet,genex)
    p_p0 = load_mutation_rate(graph_newnet, score_ini)
    final_prob = ppr_model(p_p0,graph_newnet,alpha0,modelx)
    rank = final_prob.index[final_prob['genes']==genex].values[0]+1
    return [network_size,n_neigh,rank,final_prob]

### run pipeline of ppr
def run_pipeline(gm,genex,s_list,score_ini,alpha0,modelx):  # modelx can be pr or ppr
    outs = pipeline_ppr(s_list,gm,genex,score_ini,alpha0,modelx)
    print(f"model: {modelx}, network size: {outs[0]}, completed!")
    return outs[3]


### run pipeline of unequal weight markov
def run_pipeline_unequal(gm,genex,s_list,score_ini,alpha,Wm=0.6,modelx='Markov'):
    outs = pipeline_uequal_markov(s_list,gm,genex,score_ini,Wm,alpha)
    print(f"model:{modelx}, network size: {outs[0]}, completed!")
    return outs[3]

### some functions: number of first neighbors, load mutation rate (patient count that altered)
# get the number of first neighbors of a selected gene
def get_genex_neighbors(gm_seed,genex):  # input adjacency matrix
    n_neigh = gm_seed.loc[genex,].sum().astype('int')
    return n_neigh

# get the mutation rate of the seed list (patient count that altered)
def load_mutation_rate(gm_seed, score_ini):
    # load patient mutation percentage table
    p_p = [] # percent of patients with mutated genes
    genes = gm_seed.index.tolist()
    for gene_A in genes:
        if gene_A in score_ini:
            p_p.append(score_ini[gene_A])
        else:
            p_p.append(0)
    return p_p
