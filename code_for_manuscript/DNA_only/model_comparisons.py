####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes compare different models using tier 1 genes
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums

####################################################################################
def model_comparisons_pipeline(dir0,save_dir,t_types,topNs):

    ### cgc tier 1 genes, filtered by hotspot mutation for each cancer type
    cgc_dir = os.path.join(dir0, 'databases/COSMIC','cancer_specific_tier1_hotspot_annotated')
    fig_dir = os.path.join(save_dir,'figures')
    mat_dir = os.path.join(save_dir,'output_data_matrix')

    ### evaluate pr,ppr models
    model0s = ['pr','ppr']
    for model0 in model0s:
        eval_pr_ppr(model0,t_types,save_dir,cgc_dir,fig_dir,mat_dir)
    ### evaluate diffusion models
    eval_diffusion(t_types,save_dir,cgc_dir,fig_dir,mat_dir)

    ### load model parameters with best performance
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk_MC = pd.read_csv(fnm, index_col = 0)
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_pr.csv')
    df_rnk_pr = pd.read_csv(fnm, index_col = 0)
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_ppr.csv')
    df_rnk_ppr = pd.read_csv(fnm, index_col = 0)
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_raw_diffusion.csv')
    df_rnk_raw = pd.read_csv(fnm, index_col = 0)
    for topN in topNs:
        model_compare(t_types,df_rnk_MC,df_rnk_pr,df_rnk_ppr,df_rnk_raw,cgc_dir,save_dir,fig_dir,topN = topN)


def eval_pr_ppr(model0,t_types,save_dir,cgc_dir,fig_dir,mat_dir,top_m = 15,topN = 20,fig_size = (5,3),
               fn_nums = np.arange(50,650,100),alphas = np.round(np.arange(0.05,1,0.05),3)):
    """
    Plot the top_m ranked models by average IoU
    IoU was computed between the topN genes and Tier 1 genes
    """
    ### model directory, pr/ppr
    md_dir = os.path.join(save_dir, t_types[0], 'other models detailed')
    files = os.listdir(md_dir)
    N_m = len(fn_nums)*len(alphas)
    m_df = pd.DataFrame(np.zeros([N_m,len(t_types)]))
    iou_df = pd.DataFrame(np.zeros([N_m,len(t_types)]))
    m_nm = []
    for fn_num in fn_nums:
        for alpha in alphas:
            m_nm.append('pgm_'+model0+'_wes_fn_'+str(fn_num)+'_alpha_'+str(alpha))  # model parameters
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
        for item in m_nm:
            md_dir = os.path.join(save_dir, t_type, 'other models detailed')
            file = item+'_022123.csv'
            df = pd.read_csv(os.path.join(md_dir,file))
            
            # compute the overlap between top ranked genes and tier 1 genes
            top_g = df.iloc[:topN,0].values.tolist()
            overlap = len(set(gt).intersection(top_g))
            iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
            m_df.loc[item,t_type] = overlap
            iou_df.loc[item,t_type] = iou
    ### save iou dataframe
    fnm = os.path.join(mat_dir,'iou_predicted_tier1_model_parameters_all_cancers_by_'+model0+'.csv')
    iou_df.to_csv(fnm)
    ### iou
    g = sns.clustermap(iou_df,dendrogram_ratio=0.1, cbar_pos=(0.01, 0.9, .02, .05),figsize=(5,5))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 7, rotation = 90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6, rotation = 0)
    g.ax_heatmap.set_ylabel('Model parameters', fontsize = 10)
    g.ax_cbar.tick_params(labelsize = 7)
    title = f'IoU'
    # plt.tight_layout()
    plt.title(title, fontsize = 10)
    fnm = os.path.join(fig_dir,'iou_all_cancers_tier1_hotspot_annotated_by_'+model0+'.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
    
    ### rank of models for every cancer type
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
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_'+model0+'.csv')
    df_rnk.to_csv(fnm)
    # plot rank of top_m ranked model parameters
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = fig_size
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=5, markersize=3)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=2)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 6)
    plt.title(f'Top {top_m} Models', fontsize = 12)
    plt.ylabel('Rank of Model', fontsize = 10)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of top {top_m} models all cancers_by_'+model0+'.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
    
    # plot rank of all model parameters
    top_m = df_rnk.shape[0]
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = fig_size
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=2, markersize=2)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=0.7)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 3.5)
    plt.yticks(fontsize = 6)
    plt.title(f'All Models', fontsize = 12)
    plt.ylabel('Rank of Model', fontsize = 10)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of all models all cancers_by_'+model0+'.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()

def eval_diffusion(t_types,save_dir,cgc_dir,fig_dir,mat_dir,top_m = 15,topN = 20,fig_size = (5,3),
                  fn_nums = np.arange(50,650,100),sigma_2s = np.arange(5,35,5)):
    ### model directory, raw diffusion
    md_dir = os.path.join(save_dir, t_types[0], 'other models detailed')
    files = os.listdir(md_dir)
    N_m = len(fn_nums)*len(sigma_2s)
    m_df = pd.DataFrame(np.zeros([N_m,len(t_types)]))
    iou_df = pd.DataFrame(np.zeros([N_m,len(t_types)]))
    m_nm = []
    for fn_num in fn_nums:
        for sigma_2 in sigma_2s:
            m_nm.append('diffusion_model_wes_fn_'+str(fn_num)+'_raw_sigma2='+str(sigma_2))  # model parameters
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
        for item in m_nm:
            md_dir = os.path.join(save_dir, t_type, 'other models detailed')
            file = item+'.csv'
            df = pd.read_csv(os.path.join(md_dir,file))
            
            # compute the overlap between top ranked genes and tier 1 genes
            top_g = df.iloc[:topN,1].values.tolist()
            overlap = len(set(gt).intersection(top_g))
            iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
            m_df.loc[item,t_type] = overlap
            iou_df.loc[item,t_type] = iou
    ### save iou dataframe
    fnm = os.path.join(mat_dir,'iou_predicted_tier1_model_parameters_all_cancers_by_raw_diffusion.csv')
    iou_df.to_csv(fnm)
    ### iou
    g = sns.clustermap(iou_df,dendrogram_ratio=0.1, cbar_pos=(0.01, 0.9, .02, .05),figsize=(6,5))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 7, rotation = 90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6, rotation = 0)
    g.ax_heatmap.set_ylabel('Model parameters', fontsize = 10)
    g.ax_cbar.tick_params(labelsize = 7)
    title = f'IoU'
    # plt.tight_layout()
    plt.title(title, fontsize = 10)
    fnm = os.path.join(fig_dir,'iou_all_cancers_tier1_hotspot_annotated_by_raw_diffusion.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
    
    ### rank of models for every cancer type
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
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_raw_diffusion.csv')
    df_rnk.to_csv(fnm)
    # plot rank of top_m ranked model parameters
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = fig_size
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=5, markersize=3)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=2)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 6)
    plt.title(f'Top {top_m} Models', fontsize = 12)
    plt.ylabel('Rank of Model', fontsize = 10)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of top {top_m} models all cancers_by_raw_diffusion.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
    
    # plot rank of all model parameters
    top_m = df_rnk.shape[0]
    df_rnk_top = df_rnk.iloc[:top_m,]
    # create stacked errorbars:
    means = np.array(df_rnk_top['mean_rank'].values)
    std = np.array(df_rnk_top['std_rank'].values)
    mins = np.array(df_rnk_top['min_rank'].values)
    maxs = np.array(df_rnk_top['max_rank'].values)
    plt.rcParams["figure.figsize"] = fig_size
    plt.errorbar(np.arange(top_m), means, std, fmt='ok', lw=2, markersize=2)
    plt.errorbar(np.arange(top_m), means, [means - mins, maxs - means],
                fmt='.k', ecolor='orange', lw=0.7)
    plt.xlim(-1, top_m)
    plt.xticks(np.arange(top_m),df_rnk_top.index,rotation = 90, fontsize = 3.5)
    plt.yticks(fontsize = 6)
    plt.title(f'All Models', fontsize = 12)
    plt.ylabel('Rank of Model', fontsize = 10)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'rank of all models all cancers_by_raw_diffusion.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()

def model_compare(t_types,df_rnk_MC,df_rnk_pr,df_rnk_ppr,df_rnk_raw,cgc_dir,save_dir,fig_dir,topN = 20,fig_size = (2,1.5)):
### compare Markov model, pr, ppr, and diffusion models
    df_rnk_MC = df_rnk_MC.sort_values(by = 'mean_rank')
    select_md_MC = df_rnk_MC.index[0]
    df_rnk_pr = df_rnk_pr.sort_values(by = 'mean_rank')
    select_md_pr = df_rnk_pr.index[0]
    df_rnk_ppr = df_rnk_ppr.sort_values(by = 'mean_rank')
    select_md_ppr = df_rnk_ppr.index[0]
    df_rnk_raw = df_rnk_raw.sort_values(by = 'mean_rank')
    select_md_raw = df_rnk_raw.index[0]
    iou_df = pd.DataFrame(np.zeros([len(t_types),4]))
    iou_df.index = t_types
    iou_df.columns = ['MVP','PR','PPR','raw']
    for t_type in t_types:
        cgc_fnm = os.path.join(cgc_dir,t_type+'_tier1_hotspot_filtered.csv')
        cgc_df = pd.read_csv(cgc_fnm)
        cgc_df = cgc_df.rename(columns = {'Unnamed: 0':'Gene'})
        cgc_df['rank'] = cgc_df.index+1
        gt = cgc_df.Gene.values.tolist() # 'ground truth'        
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md_MC+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))    
        # compute the overlap between top ranked genes and tier 1 genes
        top_g = df.iloc[:topN,0].values.tolist()
        overlap = len(set(gt).intersection(top_g))
        iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
        iou_df.loc[t_type,'MVP'] = iou        
        ### other models directory ###
        om_dir = os.path.join(save_dir, t_type, 'other models detailed')        
        ### pr model ###
        s_md = select_md_pr+'_022123.csv'
        df = pd.read_csv(os.path.join(om_dir,s_md))    
        # compute the overlap between top ranked genes and tier 1 genes
        top_g = df.iloc[:topN,0].values.tolist()
        overlap = len(set(gt).intersection(top_g))
        iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
        iou_df.loc[t_type,'PR'] = iou          
        ### ppr model ###
        s_md = select_md_ppr+'_022123.csv'
        df = pd.read_csv(os.path.join(om_dir,s_md))    
        # compute the overlap between top ranked genes and tier 1 genes
        top_g = df.iloc[:topN,0].values.tolist()
        overlap = len(set(gt).intersection(top_g))
        iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
        iou_df.loc[t_type,'PPR'] = iou
        ### diffusion models ###
        s_md = select_md_raw+'.csv'
        df = pd.read_csv(os.path.join(om_dir,s_md))    
        # compute the overlap between top ranked genes and tier 1 genes
        top_g = df.iloc[:topN,1].values.tolist()
        overlap = len(set(gt).intersection(top_g))
        iou = len(set(gt).intersection(top_g))/len(set(gt).union(top_g))
        iou_df.loc[t_type,'raw'] = iou
    # violin plot
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    sns.violinplot(data = iou_df,linewidth=1)
    fnm = os.path.join(fig_dir,f'top {topN} compare_different_graphical_models_violinplot_optimal.svg')
    plt.xticks(rotation = 0, fontsize = 7,fontname = "Helvetica")
    plt.yticks(fontsize = 7,fontname = "Helvetica")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.tight_layout()
    plt.savefig(fnm,dpi = 300)
    # statistical test
    a = iou_df.iloc[:,0]
    statistic = []
    p_values = []
    for i in range(1,iou_df.shape[1]):
        b = iou_df.iloc[:,i]
        statistic.append(ranksums(a,b)[0])
        p_values.append(ranksums(a,b)[1])
        print(f'p-value between MVP and {iou_df.columns[i]} for {topN} is {ranksums(a,b)[1]}')
