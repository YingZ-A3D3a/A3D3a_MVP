####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes compare study bias between models
####################################################################################

import os
import numpy as np
from matplotlib import cm
from numpy import linspace
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import seaborn as sns
from scipy.stats import ranksums

####################################################################################

def study_bias(save_dir,t_types,topNs):
    
    fig_dir = os.path.join(save_dir,'figures')
    mat_dir = os.path.join(save_dir,'output_data_matrix')

    ### load model parameters with best performance
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk_MC = pd.read_csv(fnm, index_col = 0)
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_pr.csv')
    df_rnk_pr = pd.read_csv(fnm, index_col = 0)
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_ppr.csv')
    df_rnk_ppr = pd.read_csv(fnm, index_col = 0)
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers_by_raw_diffusion.csv')
    df_rnk_raw = pd.read_csv(fnm, index_col = 0)
    plt.rcParams.update({'figure.max_open_warning': 0})
    for topN in topNs:
        study_bias_compare(t_types,df_rnk_MC,df_rnk_pr,df_rnk_ppr,df_rnk_raw,save_dir,fig_dir,topN)

def study_bias_compare(t_types,df_rnk_MC,df_rnk_pr,df_rnk_ppr,df_rnk_raw,save_dir,fig_dir,topN):
    ### best Markov model, study bias check, spearman correlation
    df_rnk_MC = df_rnk_MC.sort_values(by = 'mean_rank')
    select_md_MC = df_rnk_MC.index[0]
    df_rnk_pr = df_rnk_pr.sort_values(by = 'mean_rank')
    select_md_pr = df_rnk_pr.index[0]
    df_rnk_ppr = df_rnk_ppr.sort_values(by = 'mean_rank')
    select_md_ppr = df_rnk_ppr.index[0]
    df_rnk_raw = df_rnk_raw.sort_values(by = 'mean_rank')
    select_md_raw = df_rnk_raw.index[0]
    corr_df = pd.DataFrame(np.zeros([len(t_types),4]))
    corr_df.index = t_types
    corr_df.columns = ['MVP','PR','PPR','raw']
    # x: log2(publications), y: final score, z: final rank
    x_mc,y_mc,r_mc = {},{},{}
    x_pr,y_pr,r_pr = {},{},{}
    x_ppr,y_ppr,r_ppr = {},{},{}
    x_raw,y_raw,r_raw = {},{},{}
    for t_type in t_types:
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md_MC+'_022123.csv'
        pgm_output = pd.read_csv(os.path.join(md_dir,s_md))
        indx = pgm_output[pgm_output['publications'].isna()].index
        pgm_output.loc[indx,'publications'] = np.median(pgm_output.loc[:,'publications'].dropna()) # filled the gene with NA in publications with the median publication
        data = pgm_output.copy()
        data = data.iloc[:topN,]
        x = np.log2(data.publications.values.tolist())
        y = data.final_score.values.tolist()
        z = data.final_rank.values.tolist()
        x_mc[t_type]=x
        y_mc[t_type]=y
        r_mc[t_type]=z
        corr_df.loc[t_type,'MVP'] = spearmanr(x,z)[0]        
        ### other models directory ###
        om_dir = os.path.join(save_dir, t_type, 'other models detailed')        
        ### pr model ###    
        s_md = select_md_pr+'_022123.csv'
        pgm_output = pd.read_csv(os.path.join(om_dir,s_md))
        indx = pgm_output[pgm_output['publications'].isna()].index
        pgm_output.loc[indx,'publications'] = np.median(pgm_output.loc[:,'publications'].dropna()) # filled the gene with NA in publications with the median publication
        data = pgm_output.copy()
        data = data.iloc[:topN,]
        x = np.log2(data.publications.values.tolist())
        y = data.predicted_score.values.tolist()
        z = data.final_rank.values.tolist()
        x_pr[t_type]=x
        y_pr[t_type]=y
        r_pr[t_type]=z
        corr_df.loc[t_type,'PR'] = spearmanr(x,z)[0]        
        ### ppr model ###    
        s_md = select_md_ppr+'_022123.csv'
        pgm_output = pd.read_csv(os.path.join(om_dir,s_md))
        indx = pgm_output[pgm_output['publications'].isna()].index
        pgm_output.loc[indx,'publications'] = np.median(pgm_output.loc[:,'publications'].dropna()) # filled the gene with NA in publications with the median publication
        data = pgm_output.copy()
        data = data.iloc[:topN,]
        x = np.log2(data.publications.values.tolist())
        y = data.predicted_score.values.tolist()
        z = data.final_rank.values.tolist()
        x_ppr[t_type]=x
        y_ppr[t_type]=y
        r_ppr[t_type]=z
        corr_df.loc[t_type,'PPR'] = spearmanr(x,z)[0]        
        ### diffusion models ###
        s_md = select_md_raw+'.csv'
        df = pd.read_csv(os.path.join(om_dir,s_md))    
        data = df.copy()
        data = data.rename(columns={"rank": "final_rank"})
        data = data.iloc[:topN,]
        tmp = data.node_id.values.tolist()
        pub = []
        for gene in tmp:
            tmp1 = pgm_output.loc[pgm_output.genes==gene,'publications'].values[0]
            pub.append(tmp1)
        x = np.log2(pub)
        y = data.node_score.values.tolist()
        z = data.final_rank.values.tolist()
        x_raw[t_type]=x
        y_raw[t_type]=y
        r_raw[t_type]=z
        corr_df.loc[t_type,'raw'] = spearmanr(x,z)[0]
    make_scatter_plot(t_types,topN,x_mc,r_mc,x_raw,r_raw,x_pr,r_pr,x_ppr,r_ppr,fig_dir)
    make_violin_plot(corr_df,topN,fig_dir)

def make_scatter_plot(t_types,topN,x_mc,r_mc,x_raw,r_raw,x_pr,r_pr,x_ppr,r_ppr,fig_dir):
    start = 0.0
    stop = 1.0
    number_of_lines= len(t_types)
    cm_subsection = linspace(start, stop, number_of_lines)
    colors = [ cm.tab20(x) for x in cm_subsection ]
    ### make a plot for each cancer
    for i in range(len(t_types)):
        plt.rcParams["figure.figsize"] = (3,3)
        plt.rcParams['svg.fonttype'] = 'none'
        fig, axs = plt.subplots(2,2,sharex=True, sharey=True)
        t_type = t_types[i]        
        x0 = x_mc[t_type]
        y0 = r_mc[t_type]
        sp_mc = spearmanr(x0,y0)[0]
        x0 = x_raw[t_type]
        y0 = r_raw[t_type]
        sp_raw = spearmanr(x0,y0)[0]
        x0 = x_pr[t_type]
        y0 = r_pr[t_type]
        sp_pr = spearmanr(x0,y0)[0]
        x0 = x_ppr[t_type]
        y0 = r_ppr[t_type]
        sp_ppr = spearmanr(x0,y0)[0]        
        axs[0, 0].scatter(x_mc[t_type],r_mc[t_type],color = colors[10],label = t_type, s=2)
        axs[0, 1].scatter(x_raw[t_type],r_raw[t_type],color = colors[10],label = t_type, s=2)
        axs[1, 0].scatter(x_pr[t_type],r_pr[t_type],color = colors[10],label = t_type, s=2)
        axs[1, 1].scatter(x_ppr[t_type],r_ppr[t_type],color = colors[10],label = t_type, s=2)
    # add legend to the plot
        axs[0, 0].set_title('MVP'+'_corr_'+str(np.round(sp_mc,3)),fontsize = 5,fontname = "Helvetica")
        axs[0, 0].set_xlim(xmin = 2, xmax = 15)
        axs[0, 0].set_xticks([5,10,15])
        axs[0, 0].set_yticks([0,50,100])
        axs[0, 0].spines['top'].set_visible(False)
        axs[0, 0].spines['right'].set_visible(False)
        axs[0, 0].spines['bottom'].set_linewidth(0.75)
        axs[0, 0].spines['left'].set_linewidth(0.75)
        axs[0, 1].set_title('raw'+'_corr_'+str(np.round(sp_raw,3)),fontsize = 5,fontname = "Helvetica")
        axs[0, 1].set_xlim(xmin = 2, xmax = 15)
        axs[0, 1].set_xticks([5,10,15])
        axs[0, 1].set_yticks([0,50,100])
        axs[0, 1].spines['top'].set_visible(False)
        axs[0, 1].spines['right'].set_visible(False)
        axs[0, 1].spines['bottom'].set_linewidth(0.75)
        axs[0, 1].spines['left'].set_linewidth(0.75)
        axs[1, 0].set_title('PR'+'_corr_'+str(np.round(sp_pr,3)),fontsize = 5,fontname = "Helvetica")
        axs[1, 0].set_xlim(xmin = 2, xmax = 15)
        axs[1, 0].set_xticks([5,10,15])
        axs[1, 0].set_yticks([0,50,100])
        axs[1, 0].spines['top'].set_visible(False)
        axs[1, 0].spines['right'].set_visible(False)
        axs[1, 0].spines['bottom'].set_linewidth(0.75)
        axs[1, 0].spines['left'].set_linewidth(0.75)
        axs[1, 1].set_title('PPR'+'_corr_'+str(np.round(sp_ppr,3)),fontsize = 5,fontname = "Helvetica")
        axs[1, 1].set_xlim(xmin = 2, xmax = 15)
        axs[1, 1].set_xticks([5,10,15])
        axs[1, 1].set_yticks([0,50,100])
        axs[1, 1].spines['top'].set_visible(False)
        axs[1, 1].spines['right'].set_visible(False)
        axs[1, 1].spines['bottom'].set_linewidth(0.75)
        axs[1, 1].spines['left'].set_linewidth(0.75)       
        # loop over all axes
        for a in axs.flatten():
            a.tick_params(axis='both', which='major', labelsize=7)
            a.tick_params(axis='both', which='minor', labelsize=7)
        fnm = os.path.join(fig_dir,f'top {topN} {t_type} study_bias_final_rank_scatter_plot_model_comparisons.svg')
        plt.tight_layout()
        fig.text(0.5, 0.02, 'log2(number of publication)', fontsize = 7, ha='center')
        fig.text(0.02, 0.5, 'Rank of gene', fontsize = 7, va='center', rotation='vertical')
        plt.savefig(fnm,dpi = 300)

def make_violin_plot(corr_df,topN,fig_dir):
    plt.rcParams["figure.figsize"] = (2,1.5)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    sns.violinplot(data = corr_df,linewidth=1)
    fnm = os.path.join(fig_dir,f'top {topN} study_bias_compare_different_graphical_models_violinplot_optimal.svg')
    plt.xticks(rotation = 0,fontsize = 7,fontname = "Helvetica")
    plt.yticks(fontsize = 7,fontname = "Helvetica")
    plt.ylabel('Spearman correlation', fontsize = 7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.tight_layout()
    plt.savefig(fnm,dpi = 300)
    ### statistical test
    a = corr_df.iloc[:,0]
    statistic = []
    p_values = []
    for i in range(1,corr_df.shape[1]):
        b = corr_df.iloc[:,i]
        statistic.append(ranksums(a,b)[0])
        p_values.append(ranksums(a,b)[1])
