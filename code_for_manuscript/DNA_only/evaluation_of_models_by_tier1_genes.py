####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes optimize the MVP (Markov chain) models by tier 1 genes
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import ranksums

####################################################################################

def eval_mc_tier1(dir0,save_dir,t_types,topN):
    """
    Main script:
    * plot tier 1 genes
    * compute IOU of tier 1 genes and top 20 ranked genes for each cancer type
    * plot rank of models by average IOU
    * plot top ranked genes of optimal model for every cancer type
    * compute the overlap of top altered genes and tier 1 genes and compared that of top predicted and tier 1 genes
    """

    ### cgc tier 1 genes, filtered by hotspot mutation for each cancer type
    cgc_dir = os.path.join(dir0, 'databases/COSMIC','cancer_specific_tier1_hotspot_annotated')
    ### output figure directory
    fig_dir = os.path.join(save_dir,'figures') # figure directory
    os.makedirs(fig_dir, exist_ok=True)  # succeeds even if directory exists.
    ### output matrix directory
    mat_dir = os.path.join(save_dir,'output_data_matrix')
    os.makedirs(mat_dir, exist_ok=True)  # succeeds even if directory exists.

    ### plot tier 1 genes
    plot_tier1_genes_count(cgc_dir,t_types,fig_dir)
    ### compute IOU
    iou_df = compute_IOU(save_dir,cgc_dir,t_types,mat_dir,fig_dir,topN)
    ### rank of models
    rank_models(iou_df,mat_dir,fig_dir)

    ### plot top ranked genes of optimal model for every cancer type
    select_md = top_ranked_genes_all_cancer(t_types,save_dir,fig_dir,mat_dir,topN)
    ### compute the intersection and union of top altered and top predicted genes
    compare_top_altered_with_top_predicted(t_types,select_md,dir0,save_dir,fig_dir,mat_dir,topN)
    ### compute the overlap of top altered genes and tier 1 genes and compared that of top predicted and tier 1 genes
    intersect_top_altered_and_tier1(t_types,select_md,dir0,save_dir,cgc_dir,fig_dir,topN)

    
def plot_tier1_genes_count(cgc_dir,t_types,fig_dir):
    """ plot tier 1 gene count for each cancer type """

    num_t1 = []
    for t_type in t_types:
        cgc_fnm = os.path.join(cgc_dir,t_type+'_tier1_hotspot_filtered.csv')
        cgc_df = pd.read_csv(cgc_fnm)
        cgc_df = cgc_df.rename(columns = {'Unnamed: 0':'Gene'})
        gt = cgc_df.Gene.values.tolist() # 'ground truth'
        num_t1.append(len(gt))
    plt.rcParams["figure.figsize"] = (8,6)
    sns.barplot(y = t_types, x = num_t1, orient = 'h')
    plt.xlabel('number of tier 1 genes')
    fnm = os.path.join(fig_dir,'number of tier1 genes all cancers.png')
    plt.savefig(fnm,dpi = 256)

def compute_IOU(save_dir,cgc_dir,t_types,mat_dir,fig_dir,topN,fig_size = (5,5)):
    """ compute IOU for all models """

    ### model directory
    md_dir = os.path.join(save_dir, t_types[0], 'markov chain models')
    files = os.listdir(md_dir)
    N_m = len(files) # number of models
    m_df = pd.DataFrame(np.zeros([N_m,len(t_types)]))
    iou_df = pd.DataFrame(np.zeros([N_m,len(t_types)]))
    m_nm = []
    for item in files:
        tmp = item.split('_')
        m_nm.append('fn_'+tmp[4]+'_Wm_'+tmp[6]+'_alpha_'+tmp[8])  # model parameters
    m_df.index = m_nm
    m_df.columns = t_types
    iou_df.index = m_nm
    iou_df.columns = t_types
    for t_type in t_types:
        cgc_fnm = os.path.join(cgc_dir,t_type+'_tier1_hotspot_filtered.csv')
        cgc_df = pd.read_csv(cgc_fnm)
        cgc_df = cgc_df.rename(columns = {'Unnamed: 0':'Gene'})
        cgc_df['rank'] = cgc_df.index+1
        gt = cgc_df.Gene.values.tolist() # 'ground truth'           
        for item in files:
            md_dir = os.path.join(save_dir, t_type, 'markov chain models')
            files = os.listdir(md_dir)
            tmp = item.split('_')
            md_pr = 'fn_'+tmp[4]+'_Wm_'+tmp[6]+'_alpha_'+tmp[8]
            df = pd.read_csv(os.path.join(md_dir,item))            
            # compute the overlap between top ranked genes and tier 1 genes
            top_g = df.iloc[:topN,0].values.tolist()
            overlap = len(set(gt).intersection(top_g))
            iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
            m_df.loc[md_pr,t_type] = overlap
            iou_df.loc[md_pr,t_type] = iou
    ### save iou dataframe
    fnm = os.path.join(mat_dir,'iou_predicted_tier1_model_parameters_all_cancers.csv')
    iou_df.to_csv(fnm)
    ### make heatmap
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    g = sns.clustermap(iou_df,dendrogram_ratio=0.1, cbar_pos=(0.02, 0.9, .02, .05),figsize=fig_size)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 7, rotation = 90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6, rotation = 0)
    # g.ax_heatmap.set_xlabel('Cancer type', fontsize=10)
    g.ax_heatmap.set_ylabel('Model parameters', fontsize = 10)
    g.ax_cbar.tick_params(labelsize = 7)
    title = f'IoU'
    plt.title(title, fontsize = 10)
    fnm = os.path.join(fig_dir,'iou_all_cancers_tier1_hotspot_annotated.png')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(fig_dir,'iou_all_cancers_tier1_hotspot_annotated.svg')
    plt.savefig(fnm,dpi = 300)
    return iou_df

def rank_models(iou_df,mat_dir,fig_dir,top_m = 15, fig_size = (5,3)):
    """ 
    Rank of models for every cancer type
    Plot top 15 models or all models
    """
    rank_df = iou_df.rank(ascending = False)
    mean_rnk = rank_df.mean(axis = 1)
    std_rnk = rank_df.std(axis = 1)
    median_rnk = rank_df.median(axis = 1)
    min_rnk = rank_df.min(axis = 1)
    max_rnk = rank_df.max(axis = 1)
    df_rnk = pd.concat([mean_rnk,std_rnk],axis = 1)
    df_rnk.columns = ['mean_rank','std_rank']
    df_rnk['median_rank'] = median_rnk
    df_rnk['min_rank'] = min_rnk
    df_rnk['max_rank'] = max_rnk
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    ### save rank file
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk.to_csv(fnm)
    ### plot IoU rank of top N models
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=5, markersize=3)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=2)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 10)
    plt.title(f'Top {top_m} Models', fontsize = 14)
    plt.ylabel('Rank of Model', fontsize = 14)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of top {top_m} models all cancers.png')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(fig_dir,f'rank of top {top_m} models all cancers.svg')
    plt.savefig(fnm,dpi = 300)
    
    ### plot IoU rank of all models
    top_m = df_rnk.shape[0]
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)

    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=2, markersize=2)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=0.7)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 6)
    plt.title(f'All Models', fontsize = 12)
    plt.ylabel('Rank of Model', fontsize = 10)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of all models all cancers.png')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(fig_dir,f'rank of all models all cancers.svg')
    plt.savefig(fnm,dpi = 300)

def top_ranked_genes_all_cancer(t_types,save_dir,fig_dir,mat_dir,topN,fig_size = (15,10)):
    """ make heatmap for top 20 ranked genes for every cancer type """

    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    select_md = df_rnk.index[0]
    all_top_g = []
    top_g_ctype = {}
    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        top_g = df.iloc[:topN,0].values.tolist()
        top_g_ctype[t_type] = top_g
        all_top_g.extend(top_g)
    all_top_g = list(set(all_top_g))   
    df_top = pd.DataFrame(np.zeros([len(all_top_g),len(t_types)]))
    df_top.index = all_top_g
    df_top.columns = t_types
    for t_type in t_types:
        tmp = top_g_ctype[t_type]
        df_top.loc[tmp,t_type] = 1
    ### make heatmap
    plt.rcParams["figure.figsize"] = fig_size
    # Define colors
    colors = ((0.9, 0.9, 0.9), "#00A087")
    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
    g = sns.clustermap(df_top, method = 'ward',figsize=(15,(.2*(len(df_top.index)))),
                    cmap = cmap, cbar_pos=(0.08, 1, .01, .01),cbar_kws={"shrink": 0.01},
                    yticklabels=True, xticklabels=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 14)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'top20 genes of optimal model for all cancers.png')
    plt.savefig(fnm,dpi = 256, bbox_inches='tight')
    return select_md

def set_size(w,h,ax=None,axis_linewidth = 0.75):
    """ 
    Set format of a plot
    w, h: width, height in inches 
    """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(axis_linewidth)
    ax.spines['left'].set_linewidth(axis_linewidth)
    return ax

def intersection_barplot(fig_dir,t_types,ov_list,un_list,height0,topN,ftype = '',fig_size = [2.5,2],
                         axis_linewidth = 0.75,width = 0.5,alpha0 = 0.6,linewidth0 = 0.75, ylim0 = 0,
                         colors = ['#91D1C2FF', '#B09C85FF']):
    """
    Make intersection barplots of two sets with the same size 
    ov_list: overlap list
    un_list: union list
    height0: the size of each set, two sets have the same size
    """
    # plot data in grouped manner of bar type 
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax = set_size(fig_size[0]-0.5,fig_size[1]-0.5,ax = ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(axis_linewidth)
    ax.spines['left'].set_linewidth(axis_linewidth)

    N = len(ov_list)
    x = np.arange(0,N)

    height = [height0 for i in x]
    bottom1 = []
    overlap_percent = []
    
    for i in range(len(ov_list)):
        a = un_list[i]
        b = ov_list[i]
        bottom1.append((a-b)/2)
        overlap_percent.append(100*b/topN)
    print(f'mean overlap percent (%): {np.mean(overlap_percent)}')
    print(f'std overlap percent (%): {np.std(overlap_percent)}')
    # bottom1 = [-i/2 for i in ov_list]
    bottom2 = [0 for i in un_list]
    g1 = plt.bar(x, height, width, bottom1, color = colors[0], alpha = alpha0, edgecolor = 'k', linewidth = linewidth0)
    g2 = plt.bar(x, height, width, bottom2, color = colors[1], alpha = alpha0, edgecolor = 'k', linewidth = linewidth0)

    plt.legend([g1,g2],['MVP','Top altered'],fontsize = 6)

    plt.xticks(x,t_types,fontsize = 6, rotation = 90)
    plt.yticks(fontsize = 6)
    if ylim0!=0:
        plt.ylim([0,ylim0])
    # plt.title(f'Top {topN} genes',fontsize = 7)
    plt.ylabel('Count',fontsize = 7)

    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'Intersection_and_union_model_and_top_{topN}_altered_{ftype}_barplot.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(fig_dir,f'Intersection_and_union_model_and_top_{topN}_altered_{ftype}_barplot.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
    
def compare_top_altered_with_top_predicted(t_types,select_md,dir0,save_dir,fig_dir,mat_dir,topN):
    """ compare top 20 altered genes vs. top 20 predicted genes """
    ### make intersection and union plots
    ov_list = [] # overlap
    un_list = [] # union
    colors = ['#91D1C2FF', '#B09C85FF']
    pcg_dir = os.path.join(dir0, 'databases/TCGA_GTEx','TcgaTargetGTEX_gene_biotype_by_biomaRt.csv')
    fx = pd.read_csv(pcg_dir)
    pcg = fx.loc[fx.gene_biotype=='protein_coding','hgnc_symbol'].values.tolist()
    
    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        list_m = df.iloc[:topN,0].values.tolist() # top-ranked genes by model

        fq_dir = os.path.join(save_dir, t_type, 'TCGA_'+t_type+'_wes_freq_long.csv')
        df_q = pd.read_csv(fq_dir)
        df_q = df_q.loc[df_q.Gene.isin(pcg),]
        df_q = df_q.sort_values(by = ['Freq'],ascending = False)
        df_q = df_q.reset_index(drop = True)
        list_q = df_q.iloc[:topN,0].values.tolist() # top-ranked genes by model

        ov = len(set(list_m).intersection(list_q)) # overlap
        un = len(set(list_m).union(list_q)) # union
        ov_list.append(ov)
        un_list.append(un)
    intersection_barplot(fig_dir,t_types,ov_list,un_list,height0 = topN, topN = topN, ftype = 'MVP and top altered',ylim0 = 80)

def intersect_top_altered_and_tier1(t_types,select_md,dir0,save_dir,cgc_dir,fig_dir,topN):
    """ Compare intersection of top ranked genes and Tier 1 genes vs. top altered genes and Tier 1 genes"""
    # load protein coding genes
    pcg_dir = os.path.join(dir0, 'databases/TCGA_GTEx','TcgaTargetGTEX_gene_biotype_by_biomaRt.csv')
    fx = pd.read_csv(pcg_dir)
    pcg = fx.loc[fx.gene_biotype=='protein_coding','hgnc_symbol'].values.tolist()
    
    iou_df = pd.DataFrame(np.zeros([len(t_types),2]))
    iou_df.index = t_types
    iou_df.columns = ['MVP','Top altered']
    
    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        top_g = df.iloc[:topN,0].values.tolist()
        fq_dir = os.path.join(save_dir, t_type, 'TCGA_'+t_type+'_wes_freq_long.csv')
        df_q = pd.read_csv(fq_dir)
        df_q = df_q.loc[df_q.Gene.isin(pcg),]
        df_q = df_q.sort_values(by = 'Freq',ascending = False)
        df_q = df_q.reset_index(drop = True)
        top_f = df_q.iloc[:topN,0].values.tolist()
        
        cgc_fnm = os.path.join(cgc_dir,t_type+'_tier1_hotspot_filtered.csv')
        cgc_df = pd.read_csv(cgc_fnm)
        cgc_df = cgc_df.rename(columns = {'Unnamed: 0':'Gene'})
        cgc_df['rank'] = cgc_df.index+1
        gt = cgc_df.Gene.values.tolist() # 'ground truth'
        
        iou_g = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
        iou_f = len(set(gt).intersection(top_f))/len(set(gt).union(top_f))
        iou_df.loc[t_type,'MVP'] = iou_g
        iou_df.loc[t_type,'Top altered'] = iou_f
    ### compare the intersection of top predicted genes with tier 1 genes and that of top altered with tier 1 genes
    plt.rcParams["figure.figsize"] = (1.5,2)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    colors = ['#91D1C2FF', '#B09C85FF']
    my_pal = {"MVP": colors[0], "Top altered":colors[1]}
    sns.violinplot(data = iou_df, palette=my_pal,linewidth=1)
    
    plt.xticks(rotation = 90,fontsize = 7,fontname = "Helvetica")
    plt.yticks(fontsize = 7,fontname = "Helvetica")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,'compare_MVP_and_top_altered_violinplot.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(fig_dir,'compare_MVP_and_top_altered_violinplot.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
    
    ### statistical test
    a = iou_df.iloc[:,0]
    b = iou_df.iloc[:,1]
    statistic = ranksums(a,b,alternative = 'greater')[0]
    p_value = ranksums(a,b,alternative = 'greater')[1]
    print(f'mean of IoU of tier 1 genes with top ranked genes by MVP: {np.mean(a)}')
    print(f'mean of IoU of tier 1 genes with top ranked genes by freq: {np.mean(b)}')
    print(f'ratio of the mean: {np.mean(a)/np.mean(b)}')
    print(f'statistic={statistic},p-value={p_value}')
