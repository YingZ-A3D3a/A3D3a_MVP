####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes evaluate the MC models by tier1 genes
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums
from matplotlib.colors import LinearSegmentedColormap

####################################################################################
def eval_mc_tier1(dir0,save_dir,t_types,topN,top_m = 15,fig_size = (5,5)):

    ### cgc tier 1 genes, filtered by hotspot mutation for each cancer type
    cgc_dir = os.path.join(dir0, 'databases/COSMIC','cancer_specific_tier1_hotspot_annotated')
    fig_dir = os.path.join(save_dir,'figures')
    os.makedirs(fig_dir,exist_ok = True)

    ### output matrix directory
    mat_dir = os.path.join(save_dir,'output_data_matrix')
    os.makedirs(mat_dir,exist_ok = True)

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

    ### iou
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

    ### rank of models for every cancer type
    rank_df = iou_df.rank(ascending = False)
    rank_df.head()
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
    df_rnk.head()

    ### save rank file
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk.to_csv(fnm)
    # top 15 models by IoU
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = (10,8)
    plt.figure()
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=6)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=2)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 10)
    plt.title(f'Top {top_m} Models', fontsize = 14)
    plt.ylabel('Rank of Model', fontsize = 14)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of top {top_m} models all cancers.png')
    plt.savefig(fnm,dpi = 256)

    # all models ranked by IoU
    top_m = df_rnk.shape[0]
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = (10,8)
    plt.figure()
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=6)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=2)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 6)
    plt.title(f'All Models', fontsize = 14)
    plt.ylabel('Rank of Model', fontsize = 14)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of all models all cancers.png')
    plt.savefig(fnm,dpi = 256)

    ### top N ranked genes for all cancer types, heatmap
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    select_md = df_rnk.index[0]
    all_top_g = []
    top_g_ctype = {}
    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes+rna_'+select_md+'_022123.csv'
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
    plt.rcParams["figure.figsize"] = (15,10)
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
    plt.show()

    ### top N altered genes vs. top N predicted genes
    all_top_g = [] # all top predicted genes
    top_g_ctype = {} # top predicted genes
    top_f_ctype = {} # top altered genes
    overlap_ctype = {} # overlap between top predicted genes and top altered genes
    union_ctype = {} # union between top predicted genes and top altered genes
    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes+rna_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))

        top_g = df.iloc[:topN,0].values.tolist()
        top_g_ctype[t_type] = top_g
        all_top_g.extend(top_g)
        
        fq_dir = os.path.join(save_dir, t_type, 'elastic net model', 'wes+rna_classifier_022123.csv')
        df_q = pd.read_csv(fq_dir)
        top_f = df_q.iloc[:topN,0].values.tolist()
        overlap_ctype[t_type] = list(set(top_g).intersection(top_f))
        union_ctype[t_type] = list(set(top_g).union(top_f))

    ### number of overlapped genes
    num_u = []
    num_o = []
    for t_type in t_types:
        num_u.append(len(union_ctype[t_type]))
        num_o.append(len(overlap_ctype[t_type]))
    plt.rcParams["figure.figsize"] = (8,6)
    x = np.arange(0,2*len(t_types),2)
    width = 1
    # plot data in grouped manner of bar type
    plt.rcParams["figure.figsize"] = (8,6)
    plt.figure()
    plt.bar(x, num_u, width, color = '#4DBBD5')
    plt.bar(x, num_o, width, color = '#E64B35')
    plt.ylim([0,2*topN+5])
    plt.xticks(x,t_types,rotation = 90)
    plt.legend([f"union of top {topN} predicted and top {topN} altered", 
                "intersection of top {topN} predicted and top {topN} altered"])
    # sns.barplot(y = t_types, x = num_o, orient = 'h')
    plt.xlabel('intersection between top 20 predicted and top 20 altered genes')
    print(f'mean: {np.mean(num_o)}, std: {np.std(num_o)}')
    fnm = os.path.join(fig_dir,'intersection and union between top 20 predicted and top 20 altered genes.png')
    plt.savefig(fnm,dpi = 256)
    plt.show()

    ### intersection of top altered genes and tier 1 genes
    iou_df = pd.DataFrame(np.zeros([len(t_types),2]))
    iou_df.index = t_types
    iou_df.columns = ['MC','top_altered']

    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes+rna_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        top_g = df.iloc[:topN,0].values.tolist()
        fq_dir = os.path.join(save_dir, t_type, 'elastic net model', 'wes+rna_classifier_022123.csv')
        df_q = pd.read_csv(fq_dir)
        top_f = df_q.iloc[:topN,0].values.tolist()
        cgc_fnm = os.path.join(cgc_dir,t_type+'_tier1_hotspot_filtered.csv')
        cgc_df = pd.read_csv(cgc_fnm)
        cgc_df = cgc_df.rename(columns = {'Unnamed: 0':'Gene'})
        cgc_df['rank'] = cgc_df.index+1
        gt = cgc_df.Gene.values.tolist() # 'ground truth'
        iou_g = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
        iou_f = len(set(gt).intersection(top_f))/len(set(gt).union(top_f))
        iou_df.loc[t_type,'MC'] = iou_g
        iou_df.loc[t_type,'top_altered'] = iou_f

    plt.rcParams["figure.figsize"] = (4,5)
    iou_df.boxplot(grid = False, rot = 90)
    plt.ylabel('IoU of tier 1 genes and top 20 genes')
    fnm = os.path.join(fig_dir,'compare_MC_and_top_altered.png')
    plt.xticks(rotation = 90, fontsize = 8)
    plt.yticks(fontsize = 8)
    plt.tight_layout()
    plt.savefig(fnm,dpi = 256)
    plt.show()

    ### statistical test
    a = iou_df.iloc[:,0]
    b = iou_df.iloc[:,1]
    statistic = ranksums(a,b,alternative = 'greater')[0]
    p_value = ranksums(a,b,alternative = 'greater')[1]
    print(f'statistic={statistic},p-value={p_value}')
