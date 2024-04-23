### prepare the unique interactome (one edge between any two nodes) for graphical modeling
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def adj_mat():
    ### load interactome
    interactome = os.path.join(os.path.dirname(__file__), 'data', 'unique_interactome_CanSAR_KEGG_Hi_union_021524.txt')
    interacts = pd.read_csv(interactome, sep = '\t')

    # get the edges
    edges = []
    for i in range(len(interacts))a:
        edges.append(tuple(interacts.iloc[i,:].values.tolist()))
    # print(f"number of edges: {len(edges)}")

    ### convert adjacency list to matrix
    # get the nodes

    geneA = list(set(interacts['geneA'].values))
    geneB = list(set(interacts['geneB'].values))
    geneA.extend(geneB)
    nodes = sorted([str(i) for i in list(set(geneA))])
    # print(f"number of nodes: {len(nodes)}")
    print('Interactome loaded!')

    ### convert lists of edges to adjacency matrix

    g = nx.Graph()
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    gm = nx.to_pandas_adjacency(g)
    
    return gm
