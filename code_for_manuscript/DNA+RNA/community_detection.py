####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes perform community detection
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import seaborn as sns
import networkx as nx
import PGM_PY.graph_builder as gb
import PGM_PY.adj_mat_interactome as ami
import PGM_PY.graphical_models as gms # run graphical model
import PGM_PY.markov_model_dataframe as mmd
import math
import networkx.algorithms.community as nx_comm
import random
### enrichment analysis
import gseapy
import time

####################################################################################
def community_detection_pipeline(dir0,save_dir,t_types,resolution = 0.7, seed0 = 12345,fdr_thre = 0.05,plot_topN = 60):

    ### load Markov model parameters with best performance
    fig_dir = os.path.join(save_dir,'figures')
    mat_dir = os.path.join(save_dir,'output_data_matrix')
    
    ### additional data directory (for druggability annotation)
    add_data_dir = os.path.join(dir0,'druggability annotation')
    save_dir1 = os.path.join(fig_dir,'network_graph')
    os.makedirs(save_dir1, exist_ok = True)
    
    ### load druggability and depmap information
    fnm = os.path.join(add_data_dir,'Depmap_essentiality_and_GDSC_targets_for_additional_information_of_druggability_top20_ranked_detailed_new.csv')
    info_df = pd.read_csv(fnm)
    
    # load additional druggability annotation for BRCA subtypes
    fnm = os.path.join(add_data_dir,'additional genes BRCA subtypes with druggabilities.csv')
    info_df1 = pd.read_csv(fnm)

    # concatenate druggability annotation together
    info_df = pd.concat([info_df,info_df1],axis = 0)
    print(info_df.shape)
    info_df.head()
    
    gm = ami.adj_mat() # adjacency matrix of the interactome

    ### load complete depmap essentials for every cancer type
    fnm = os.path.join(mat_dir,'depmap essential genes for every tumor type.csv')
    dep_df = pd.read_csv(fnm)
    
    ### load GDSC targets of top ranked genes for every cancer type
    fnm = os.path.join(mat_dir,'compare_targets_in_GDSC_top_ranked_MVP_and_top_100_altered.csv')
    gdsc_df = pd.read_csv(fnm)
    gdsc_df.head()
    
    ### load GDSC targets of all genes in the network for every cancer type
    fnm = os.path.join(mat_dir,'compare_targets_in_GDSC_all_network_genes.csv')
    gdsc_df1 = pd.read_csv(fnm)
    gdsc_df1.head()

    ### load best Markov model, plot the top ranked genes
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    select_md = df_rnk.index[0]
    
    for t_type in t_types:
        print(t_type)
        print('---------------------------------------------')
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes+rna_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        t_subtype = t_type
        plot_the_top_ranked_disease_network(df,gm,info_df,dep_df,gdsc_df,t_type,t_subtype,fig_dir,top_k=plot_topN)

    comm_hallmark_all_t_types = community_detection(mat_dir,gdsc_df1,dep_df,fig_dir,save_dir,t_types,gm,resolution,seed0) # all cancers
    
def nudge(pos, x_shift, y_shift):
    return {n:(x + x_shift, y + y_shift) for n,(x,y) in pos.items()}

### plot graph with druggability info

"""
To run this function, it requires a info_df file that contains the druggability class of top ranked genes
a dep_df file that contains the depmap essential genes of each cancer type
and a gdsc_df file that contains the GDSC targets of top ranked genes
"""
def plot_the_top_ranked_disease_network(df,gm,info_df,dep_df,gdsc_df,t_type,t_subtype,fig_dir,top_k=60,
                                        topN=20,
                                        edge_width = 1,line_width = 0.08,
                                        gray = '#C4C3D0',orange = '#FFD700',
                                        lightgreen = '#98FB98', green = '#03AC13',
                                        lightblue = '#87CEEB',blue = '#0030BF',
                                        yellow = '#CC7722',
                                        purple = '#9867C5',
                                        fig_size = [4, 3.5]):
    
    # plot top_k genes with topN genes with druggability info
    
    final_prob_markov = df.copy()
    top_g = final_prob_markov.genes[:topN].values.tolist()  # topN genes

    # top_k = final_prob_markov.shape[0]
    genes_top = final_prob_markov.genes[:top_k].values.tolist()
    g_adj = gb.seed_network(genes_top,gm)
    g = gb.build_graph(g_adj)
    print(f'{len(g.nodes)} nodes are connected among top {top_k} nodes')
    
    nodes_1 = [] # seed
    nodes_2 = [] # seed, RNA
    nodes_3 = [] # first neighbor

    node_colorcode = []

    for k in g.nodes:
        if final_prob_markov.loc[final_prob_markov.genes==k,'source'].values[0]=='WES/SV':
            nodes_1.append(k)
        elif final_prob_markov.loc[final_prob_markov.genes==k,'source'].values[0]=='RNA':
            nodes_2.append(k)
        elif final_prob_markov.loc[final_prob_markov.genes==k,'source'].values[0]=='first neighbor':   # other first neighbors
            nodes_3.append(k)
            
    if len(nodes_1)+len(nodes_2)+len(nodes_3)<len(g.nodes):
        logging.error('Please check the number of node types.')
        return

    plt.rcParams['figure.figsize'] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    ax = plt.axes()
    # pos = nx.spring_layout(g, seed = 500, k = 0.95)
    pos = nx.kamada_kawai_layout(g)
    # size
    slope = 0.8
    intercept = 1

    max_wt = max(final_prob_markov.iloc[:top_k,7].values) # final rank
    wt1 = [slope*(max_wt-final_prob_markov.loc[final_prob_markov.genes==k,'final_rank'].values[0])+intercept for k in nodes_1]
    wt2 = [slope*(max_wt-final_prob_markov.loc[final_prob_markov.genes==k,'final_rank'].values[0])+intercept for k in nodes_2]
    wt3 = [slope*(max_wt-final_prob_markov.loc[final_prob_markov.genes==k,'final_rank'].values[0])+intercept for k in nodes_3]
    
    # node face color
    cr1 = []
    face_c = {'A':lightblue,'~A':lightblue,'B':lightgreen,'~B':lightgreen,'C':orange,'~C':orange,'D':purple,'~D':purple}
    for gene in nodes_1:
        if gene in final_prob_markov.iloc[:topN,0].values: # if k is among the top N genes
            class0 = info_df.loc[info_df.gene_name==gene, 'Class'].values[0]
            if class0 in face_c:
                cr1.append(face_c[class0])
            else:
                cr1.append('#CCCCCC')
        else:
            cr1.append('w')
    
    cr2 = []
    for gene in nodes_2:
        if gene in final_prob_markov.iloc[:topN,0].values: # if k is among the top N genes
            class0 = info_df.loc[info_df.gene_name==gene, 'Class'].values[0]
            if class0 in face_c:
                cr2.append(face_c[class0])
            else:
                cr2.append('#CCCCCC')
        else:
            cr2.append('w')
            
    cr3 = [] # for nodes_3
    for gene in nodes_3:
        if gene in final_prob_markov.iloc[:topN,0].values: # if k is among the top N genes
            class0 = info_df.loc[info_df.gene_name==gene, 'Class'].values[0]
            if class0 in face_c:
                cr3.append(face_c[class0])
            else:
                cr3.append('#CCCCCC')
        else:
            cr3.append('w')

    # depmap essential, including common essentials
    ess = dep_df.loc[dep_df['tumor type']==t_type,'depmap essentials'].values.tolist()
    
    # GDSC targets
    if t_type in set(gdsc_df['Tumor type']):
        tmp = gdsc_df.loc[(gdsc_df['Tumor type']==t_type)&(gdsc_df['Method']=='MVP'),'GDSC targets among top 100 with low z-scored LN_IC50'].values[0]
        gdsc_t = tmp.split(',')
    else:
        gdsc_t = []
    
    # node edge color
    edgecolors1 = []
    for k in nodes_1:
        if (k in ess) and (k in gdsc_t):
            edgecolors1.append(green)
        else:   
            if k in ess:
                edgecolors1.append(yellow)
            elif k in gdsc_t:
                edgecolors1.append(blue)
            else:
                edgecolors1.append('black')

    edgecolors2 = []
    for k in nodes_2:
        if (k in ess) and (k in gdsc_t):
            edgecolors2.append(green)
        else:   
            if k in ess:
                edgecolors2.append(yellow)
            elif k in gdsc_t:
                edgecolors2.append(blue)
            else:
                edgecolors2.append('black')
                
    edgecolors3 = []
    for k in nodes_3:
        if (k in ess) and (k in gdsc_t):
            edgecolors3.append(green)
        else:   
            if k in ess:
                edgecolors3.append(yellow)
            elif k in gdsc_t:
                edgecolors3.append(blue)
            else:
                edgecolors3.append('black')
        nx.draw_networkx_nodes(g,pos,nodelist=nodes_1,edgecolors=edgecolors1,linewidths=edge_width,node_color = cr1,node_shape='o',node_size=wt1) # WES          
    nx.draw_networkx_nodes(g,pos,nodelist=nodes_2,edgecolors=edgecolors2,linewidths=edge_width,node_color = cr2,node_shape='d',node_size=wt2) # RNA
    nx.draw_networkx_nodes(g,pos,nodelist=nodes_3,edgecolors=edgecolors3,linewidths=edge_width,node_color = cr3,node_shape='s',node_size=wt3) # first neighbor     

    pos_nodes = nudge(pos, 0, 0.07)  
    nx.draw_networkx_edges(g,pos,width = line_width,edge_color='#BBBBBB',alpha = 0.6)
    labels = nx.draw_networkx_labels(g,pos=pos_nodes,font_size = 5,font_family = 'sans-serif')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
   
    plt.title(t_subtype, fontsize = 5)
    plt.tight_layout()
    
    fnm = os.path.join(fig_dir,f'{t_subtype} top {top_k} disease network with {topN} genes with druggability.svg')
    plt.savefig(fnm,dpi = 300, bbox_inches='tight')
    fnm = os.path.join(fig_dir,f'{t_subtype} top {top_k} disease network with {topN} genes with druggability.png')
    plt.savefig(fnm,dpi = 300, bbox_inches='tight')
    plt.show()

def community_detection(mat_dir,gdsc_df1,dep_df,fig_dir,save_dir,t_types,gm,resolution,seed0):
    ### community detection
    # load best Markov model, set the initial score of each of topN genes to zero and re-run the Markov model
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    select_md = df_rnk.index[0]
    temp = select_md.split('_')
    Wm = float(temp[3])
    alpha = float(temp[5])

    comm_hallmark_all_t_types = pd.DataFrame()
    comm_all = []
    t_types_all = []
    g_comm_all = []
    for t_type in t_types:
        print(t_type)
        print('---------------------------------------------')
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes+rna_'+select_md+'_022123.csv'
        final_prob_markov = pd.read_csv(os.path.join(md_dir,s_md))
        df0 = final_prob_markov.copy()
        df0.index = df0.genes
        genes_top = final_prob_markov.genes
        Ncg = len(genes_top)
        g_adj = gb.seed_network(genes_top,gm)
        g = gb.build_graph(g_adj)

        comms = nx_comm.louvain_communities(g, resolution = resolution, seed=seed0)
        # sorting on basis of size of list
        comms.sort(key = lambda x:len(x), reverse = True)

        enri_mat = pd.DataFrame() 
        for i in range(len(comms)):
            print(f"community {i+1}, size={len(comms[i])}")
            comm = comms[i]
            t_types_all.extend([t_type for j in comm])
            comm_all.extend([i+1 for j in comm])
            g_comm_all.extend(list(comm))
            
            print(list(comm)[:5])
            Hallmark_depmap = gseapy.enrichr(gene_list=list(comm),
                gene_sets=['MSigDB_Hallmark_2020'],
                organism='Human',
                cutoff=0.5)
            tmp = Hallmark_depmap.res2d
            time.sleep(0.5)
            tmp['community'] = [i+1 for k in range(tmp.shape[0])]
            enri_mat = pd.concat([enri_mat,tmp],axis = 0)

            ### plot top ranked genes in each community
            genex = list(comm)[0]
            ### run graphical model on a community
            score_ini = {}
            for item in comm:
                score_ini[item] = final_prob_markov.loc[final_prob_markov.genes==item,'initial_score'].values[0]
            final_prob_markov0 = gms.run_pipeline_unequal(gm,genex,comm,score_ini,alpha,Wm,modelx='Markov')
            
            s_list = comm
            final_prob_markov_sub = mmd.info_markov(final_prob_markov0,s_list,gm)
            
            final_prob_markov_sub['source'] = df0.loc[final_prob_markov_sub.genes,'source'].values.tolist()
            comm_id = i+1
            t_subtype = t_type
            plot_the_top_ranked_genes_in_community(final_prob_markov_sub,gm,gdsc_df1,dep_df,t_type,t_subtype,comm_id,fig_dir)
            
        enri_mat['tumor_type'] = [t_type for k in range(enri_mat.shape[0])]
        comm_hallmark_all_t_types = pd.concat([comm_hallmark_all_t_types,enri_mat],axis = 0)
    ### save community matrix
    comm_all_df = pd.DataFrame({'tumor_type':t_types_all,'community':comm_all,'gene':g_comm_all})
    fnm = os.path.join(mat_dir,'community_genes_all_cancers.csv')
    comm_all_df.to_csv(fnm,index = False)
    ### save hallmark pathways for every community for all cancer types
    fnm = os.path.join(mat_dir,'hallmark_pathways_community_all_cancers.csv')
    comm_hallmark_all_t_types.to_csv(fnm,index = False)
    ### make clustermap
    for cancer1 in t_types:
        make_clustermap(comm_hallmark_all_t_types,fig_dir,cancer1)
    return comm_hallmark_all_t_types

### make clustermap for comparing the enriched pathways of a selected cancer type
def make_clustermap(comm_hallmark_all_t_types,fig_dir,cancer1,fdr_thre = 0.05):
    tmp = comm_hallmark_all_t_types.loc[(comm_hallmark_all_t_types['Adjusted P-value']<fdr_thre),]
    tmp = tmp.loc[tmp.tumor_type.isin([cancer1]),]
    hallmarks = tmp.copy()
    hallmarks['tumor_community'] = hallmarks['tumor_type']+': comm_'+hallmarks['community'].astype(int).astype(str)
    hallmarks['-log10(adjusted p-value)'] = -np.log10(hallmarks['Adjusted P-value'].values.astype('float'))
    h_df = hallmarks.pivot(index = 'tumor_community',columns = 'Term',values = '-log10(adjusted p-value)')
    h_df = h_df.fillna(0)

    plt.rcParams["figure.figsize"] = (1.5,6)
    plt.rcParams['svg.fonttype'] = 'none'
    g = sns.clustermap(h_df.T,cmap = 'GnBu',xticklabels=True, yticklabels=True, cbar_pos=(1, .8, .02, .06),
                    cbar_kws={'shrink': 0.1},linewidth = 0.2,linecolor = 'gray',vmin=0, vmax=5,
                    figsize=(1.5,(0.6*(len(h_df.index)))))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 5, rotation = 90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 5)
    g.ax_heatmap.set(ylabel='hallmark pathways', xlabel='tumor community')
    title = '-log10 (adjusted p-value)'
    plt.title(title, fontsize = 5)
    fnm = os.path.join(fig_dir,f'hallmark pathways enriched each community clustermap_{cancer1}.svg')
    # plt.tight_layout()
    plt.savefig(fnm,dpi = 300, bbox_inches='tight')
    plt.show()

def plot_the_top_ranked_genes_in_community(df,gm,gdsc_df1,dep_df,t_type,t_subtype,comm_id,fig_dir,top_k=60,
                                           edge_width = 0.8,line_width = 0.08,
                                           green = '#03AC13',
                                           blue = '#0030BF',
                                           yellow = '#CC7722',
                                           node_size = 10,fig_size = [5, 4.6]):
    tmp = t_type.split('_')
    t_type0 = tmp[0]
    final_prob_markov = df.copy()

    genes_top = final_prob_markov.genes[:top_k].values.tolist()
    size_comm = final_prob_markov.shape[0]
    g_adj = gb.seed_network(genes_top,gm)
    g = gb.build_graph(g_adj)
    
    nodes_1 = []
    nodes_2 = []
    nodes_3 = []

    node_colorcode = []

    for k in g.nodes:
        if final_prob_markov.loc[final_prob_markov.genes==k,'source'].values[0]=='WES/SV':   # WES
            nodes_1.append(k)
        elif final_prob_markov.loc[final_prob_markov.genes==k,'source'].values[0]=='RNA':   # first neighbors
            nodes_2.append(k)
        elif final_prob_markov.loc[final_prob_markov.genes==k,'source'].values[0]=='first neighbor':   # first neighbors
            nodes_3.append(k)

    plt.rcParams['figure.figsize'] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    ax = plt.axes()
    # pos = nx.spring_layout(g, seed = 500, k = 0.95)
    pos = nx.kamada_kawai_layout(g)
    # size
    slope = 0
    intercept = 5

    # node facecolor
    max_cr = top_k
    cr1 = [(final_prob_markov.loc[final_prob_markov.genes==k,'final_rank'].values[0]) for k in nodes_1]
    cr2 = [(final_prob_markov.loc[final_prob_markov.genes==k,'final_rank'].values[0]) for k in nodes_2]
    cr3 = [(final_prob_markov.loc[final_prob_markov.genes==k,'final_rank'].values[0]) for k in nodes_3]

    vmin, vmax = 1, max_cr
    cmap = 'Wistia_r'

    # depmap essential
    ess = dep_df.loc[dep_df['tumor type']==t_type0,'depmap essentials'].values.tolist()
    
    # GDSC targets
    if t_type in set(gdsc_df1['Tumor type']):
        tmp = gdsc_df1.loc[gdsc_df1['Tumor type']==t_type0,'GDSC targets with low z-scored LN_IC50'].values[0]
        gdsc_t = tmp.split(',')
    else:
        gdsc_t = []
        
    # node edge color
    edgecolors1 = []
    for k in nodes_1:
        if (k in ess) and (k in gdsc_t):
            edgecolors1.append(green)
        else:   
            if k in ess:
                edgecolors1.append(yellow)
            elif k in gdsc_t:
                edgecolors1.append(blue)
            else:
                edgecolors1.append('black')

    edgecolors2 = []
    for k in nodes_2:
        if (k in ess) and (k in gdsc_t):
            edgecolors2.append(green)
        else:   
            if k in ess:
                edgecolors2.append(yellow)
            elif k in gdsc_t:
                edgecolors2.append(blue)
            else:
                edgecolors2.append('black')
                
    edgecolors3 = []
    for k in nodes_3:
        if (k in ess) and (k in gdsc_t):
            edgecolors3.append(green)
        else:   
            if k in ess:
                edgecolors3.append(yellow)
            elif k in gdsc_t:
                edgecolors3.append(blue)
            else:
                edgecolors3.append('black')

    nx.draw_networkx_nodes(g,pos,nodelist=nodes_1,edgecolors=edgecolors1,linewidths=edge_width,node_color = cr1,node_shape='o',node_size=node_size,cmap = cmap,vmin=vmin,vmax=vmax) # WES/SV
    nx.draw_networkx_nodes(g,pos,nodelist=nodes_2,edgecolors=edgecolors2,linewidths=edge_width,node_color = cr2,node_shape='d',node_size=node_size,cmap = cmap,vmin=vmin,vmax=vmax) # RNA
    nx.draw_networkx_nodes(g,pos,nodelist=nodes_3,edgecolors=edgecolors3,linewidths=edge_width,node_color = cr3,node_shape='s',node_size=node_size,cmap = cmap,vmin=vmin,vmax=vmax) # first neighbor
    
    pos_nodes = nudge(pos, 0, 0.05)
    nx.draw_networkx_edges(g,pos,width=.1,edge_color='#BBBBBB',alpha = 0.6)
    labels = nx.draw_networkx_labels(g,pos_nodes,font_color = 'k', font_weight = 'normal', font_size = 5,font_family = 'Arial')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks = np.arange(10,vmax,10),orientation = 'vertical',aspect = 5,fraction=0.04, pad=0.03, label = 'rank of genes')
    tick_font_size = 5
    cbar.ax.tick_params(labelsize=tick_font_size)
    cbar.set_label(label='rank of genes',size=5)
    plt.title(t_subtype+' Community '+str(comm_id),fontsize = 5)
    # plt.tight_layout()
    fnm = os.path.join(fig_dir,f'{t_subtype} top {top_k} in community {comm_id} disease network with genes.svg')
    plt.savefig(fnm,dpi = 300, bbox_inches='tight')
    fnm = os.path.join(fig_dir,f'{t_subtype} top {top_k} in community {comm_id} disease network with genes.png')
    plt.savefig(fnm,dpi = 300, bbox_inches='tight')
    plt.show()
