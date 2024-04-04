import numpy as np
import pandas as pd
import networkx as nx

######################################################################################################
### build networkx graph from adjacency matrix
def build_graph(adjacency_matrix):
    g = nx.from_pandas_adjacency(adjacency_matrix)
    # print(nx.info(g))
    return g

######################################################################################################
### get the largest connected component of a graph
# input a list of disease genes and the adjacency matrix of the background interactome
# return the adjacency matrix of the largest connected graph of the disease genes
def seed_network(s_list,adj_tome):
    all_nodes = adj_tome.index
    s_list_cp = [i for i in s_list if i in all_nodes]
    gm_seed = adj_tome.loc[s_list_cp,s_list_cp]
    G = build_graph(gm_seed)
    largest_cc = max(nx.connected_components(G), key=len) ## get the largest connected component
    largest_cc = sorted(list(largest_cc))
    gm_seed = adj_tome.loc[largest_cc,largest_cc]
# remove self loop
    for i in gm_seed.index:
        gm_seed.loc[i,i] = 0
    return gm_seed