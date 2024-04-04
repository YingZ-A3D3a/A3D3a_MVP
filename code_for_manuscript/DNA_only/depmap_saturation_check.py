####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes perform detailed depmap validation, including Tier1 genes, altered genes in each modality (mutation, cnv, sv, RNA), enrichment analysis, etc.
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr,spearmanr
from scipy.stats import ranksums
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection
from collections import Counter
import logging

####################################################################################
### model name
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉") # function to convert to subscript
example_string = "A3D3a"
Model_nm = example_string.translate(SUB)+"'s MVP"

def depmap_validation_detailed(dir0,save_dir,t_types,
        mt_types,max_n = 600,topN = 100,types = [Model_nm,'Top altered'],
         colors = ['#91D1C2FF', '#B09C85FF']):

    fig_dir = os.path.join(save_dir,'figures')
    save_dir1 = os.path.join(fig_dir,'depmap_validation_detailed')
    os.makedirs(save_dir1,exist_ok = True)
    
    mat_dir = os.path.join(save_dir,'output_data_matrix')
    ### load Markov model parameters with best performance
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    
    ### depmap
    dir_d = os.path.join(dir0,'databases/Depmap')
    ef = os.path.join(dir_d, 'CRISPRGeneEffect_1222.csv')
    ce = os.path.join(dir_d, 'AchillesCommonEssentialControls_1222.csv')
    model = os.path.join(dir_d, 'Model_1222.csv')
    ### genes with Chronos score
    g_ef = pd.read_csv(ef,index_col = 0)
    ind = g_ef.columns
    ind1 = [i.split(' ')[0] for i in ind]
    g_ef.columns = ind1
    ### common essential genes
    g_ce = pd.read_csv(ce)
    ce_g = []
    for item in g_ce.Gene.values:
        tmp = item.split(' ')
        ce_g.append(tmp[0])
    print(f"number of common essential genes: {len(ce_g)}")
    overlap_ce = list(set(g_ef.columns).intersection(ce_g))
    print(f"number of common essential genes that intersect with the gene effect matrix: {len(overlap_ce)}")
    ### load cell line info
    cline = pd.read_csv(model, index_col = 0)
    print(cline.shape)

    # load gene effect data
    ge_list_all = {} # gene list for each cell line for every type of cancer
    ge_list_all_rm_ce = {} # gene list for each cell line for every type of cancer, remove common essentials
    num_cell_lines = {}
    for i in range(len(mt_types)):
        ### select a cancer type of interest
        coe = mt_types[i]
        t_type = t_types[i]
        if type(coe)==list:
            mds = cline.loc[cline.OncotreeCode.isin(coe),].index.tolist()
            mds = list(set(mds).intersection(g_ef.index))
            print(f'{"+".join(coe)} ({t_type}): {len(mds)}') 
        else:
            mds = cline.loc[cline.OncotreeCode.values==coe,].index.tolist()
            mds = list(set(mds).intersection(g_ef.index))
            print(f'{coe} ({t_type}): {len(mds)}')
        num_cell_lines[t_type] = len(mds)
        ge_list = []
        ge_list_rm_ce = [] # remove common essentials        
        for i in range(len(mds)):
            if mds[i] in g_ef.index:
                glist1 = g_ef.loc[mds[i],]
                glist1 = glist1.sort_values(ascending = True)
                g_list1 = glist1.index[(glist1<-1)]
                ge_list.append(g_list1)
                tmp = [i for i in g_list1 if i not in overlap_ce]
                ge_list_rm_ce.append(tmp)                
        ge_list_all[t_type] = ge_list
        ge_list_all_rm_ce[t_type] = ge_list_rm_ce

        #################### number of essential genes ####################
    ### essential genes and pooled essential genes
    ng_rm_ce = {} # number of essential genes for each cell line for every type of cancer, remove common essentials
    for t_type in t_types:
        ng_rm_ce[t_type] = [len(i) for i in ge_list_all_rm_ce[t_type]]

    ng = {} # number of essential genes for each cell line for every type of cancer
    for t_type in t_types:
        ng[t_type] = [len(i) for i in ge_list_all[t_type]]

    ### distribution of essential genes, common essential genes removed
    ### pooled across all cell lines of a cancer type
    ### distribution of essential genes, common essential genes removed
    ng_rm_ce_pool = [] # number of essential genes pooled all cell lines for every type of cancer, remove common essentials
    g_rm_ce_pool = {} # essential genes pooled all cell lines for every type of cancer, remove common essentials
    for t_type in t_types:
        tmp = []
        for i in ge_list_all_rm_ce[t_type]:
            tmp.extend(i)
        ng_rm_ce_pool.append(len(set(tmp)))
        g_rm_ce_pool[t_type] = list(set(tmp))


    ng_pool = [] # number of essential genes pooled all cell lines for every type of cancer
    g_pool = {} # essential genes pooled all cell lines for every type of cancer

    for t_type in t_types:
        tmp = []
        for i in ge_list_all[t_type]:
            tmp.extend(i)
        ng_pool.append(len(set(tmp)))
        g_pool[t_type] = list(set(tmp))
        
    ### protein coding genes
    dir2 = 'databases/TCGA_GTEx'
    fnm = os.path.join(dir0,dir2,'TcgaTargetGTEX_gene_biotype_by_biomaRt.csv')
    fx = pd.read_csv(fnm)
    pcg = fx.loc[fx.gene_biotype=='protein_coding','hgnc_symbol'].values.tolist()
    print(f'number of protein coding genes {len(pcg)}')
    
    ### best Markov model, check depmap genes in the disease network
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    select_md = df_rnk.index[0]
    in_depmap_list_all = []
    in_depmap_list_all_rm_ce = []

    for t_type in t_types:
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        df = df.iloc[:max_n,]

        fq_dir = os.path.join(save_dir, t_type, 'TCGA_'+t_type+'_wes_freq_long.csv')
        df_q = pd.read_csv(fq_dir)
        df_q = df_q.loc[df_q.Gene.isin(pcg),]
        df_q = df_q.sort_values(by = ['Freq'],ascending = False)
        df_q = df_q.reset_index(drop = True)
        df_q = df_q.iloc[:max_n,]
        freq_q = df_q.Freq.values.tolist()

        ### compare genes in depmap among top predicted vs. top altered
        in_depmap = []
        in_depmap_q = []
        for i in range(df.shape[0]):
            tmp = df.iloc[:(i+1),0].values.tolist()
            in_depmap.append(len(set(tmp).intersection(g_pool[t_type])))

        for i in range(df_q.shape[0]):
            tmp = df_q.iloc[:(i+1),0].values.tolist()
            in_depmap_q.append(len(set(tmp).intersection(g_pool[t_type])))

        in_depmap_list = [in_depmap,in_depmap_q]
        freq_list = freq_q
        make_plots(in_depmap_list,t_type,freq_list,save_dir1 = save_dir1, types = types,topN = topN,colors = colors)
        in_depmap_list_all.append(in_depmap_list)

        ### compare genes in depmap among top predicted vs. top altered, remove common essentials
        in_depmap_rm_ce = []
        in_depmap_q_rm_ce = []
        for i in range(df.shape[0]):
            tmp = df.iloc[:(i+1),0].values.tolist()
            in_depmap_rm_ce.append(len(set(tmp).intersection(g_rm_ce_pool[t_type])))

        for i in range(df_q.shape[0]):
            tmp = df_q.iloc[:(i+1),0].values.tolist()
            in_depmap_q_rm_ce.append(len(set(tmp).intersection(g_rm_ce_pool[t_type])))

        in_depmap_list_rm_ce = [in_depmap_rm_ce,in_depmap_q_rm_ce]
        make_plots(in_depmap_list_rm_ce,t_type,freq_list,save_dir1 = save_dir1, types = types,topN = topN,colors = colors,ftype = 'remove_common_essential')
        in_depmap_list_all_rm_ce.append(in_depmap_list_rm_ce)
        
    ### make boxplot, including common essentials or excluding common essentials
    ################# compare AUC #################
    area_list = []
    for item in range(len(in_depmap_list_all)): # each item correspond to a cancer type
        in_depmap_list = in_depmap_list_all[item] # in_depmap_list is a list with at least two elements
        a = area_between_curves(in_depmap_list,topN = topN)
        area_list.append(a)
    area_mat = np.array(area_list)
    df_all = pd.DataFrame(area_mat)
    df_all.columns = types
    make_boxplot(df_all,colors = colors,types = types,topN = topN,save_dir1 = save_dir1,ftype = 'AUC_including_common_essentials')

    area_list = []
    for item in range(len(in_depmap_list_all_rm_ce)):
        in_depmap_list = in_depmap_list_all_rm_ce[item]
        a = area_between_curves(in_depmap_list,topN = topN)
        area_list.append(a)
    area_mat = np.array(area_list)
    df_all = pd.DataFrame(area_mat)
    df_all.columns = types
    make_boxplot(df_all,colors = colors,types = types,topN = topN,save_dir1 = save_dir1,ftype = 'AUC_excluding_common_essentials')

    ################# compare gene count #################
    c_list = []
    for item in range(len(in_depmap_list_all)): # each item correspond to a cancer type
        in_depmap_list = in_depmap_list_all[item] # in_depmap_list is a list with at least two elements
        tmp = []
        for list1 in in_depmap_list:
            a = list1[topN-1]
            tmp.append(a)
        c_list.append(tmp)
    c_mat = np.array(c_list)
    df_all = pd.DataFrame(c_mat)
    df_all.columns = types
    make_boxplot(df_all,colors = colors,types = types,topN = topN,save_dir1 = save_dir1, ylabel = 'Count of genetic dependencies',ftype = 'including_common_essentials')

    c_list = []
    for item in range(len(in_depmap_list_all_rm_ce)):
        in_depmap_list = in_depmap_list_all_rm_ce[item]
        tmp = []
        for list1 in in_depmap_list:
            a = list1[topN-1]
            tmp.append(a)
        c_list.append(tmp)
    c_mat = np.array(c_list)
    df_all = pd.DataFrame(c_mat)
    df_all.columns = types
    make_boxplot(df_all,colors = colors,types = types,topN = topN,save_dir1 = save_dir1, ylabel = 'Count of genetic dependencies',ftype = 'excluding_common_essentials')

### make errorbar plot, including common essentials or excluding common essentials
    make_plot_errorbar(in_depmap_list_all, topN = max_n, types = types, colors = colors, save_dir1 = save_dir1, ftype = 'including_common_essentials')

    make_plot_errorbar(in_depmap_list_all_rm_ce, topN = max_n, types = types, colors = colors, save_dir1 = save_dir1, ftype = 'excluding_common_essentials')
    
#################################### hypergeometric test ####################################
    g_entire = set(pcg).union(g_ef.columns)
    M = len(g_entire)
    print(f'Entire gene set size: {M}')
    
    ### Including_common_essentials
    pvals = [] 
    """
    pvals is a list of multiple lists of p-values
    the length of pvals should be the same as the length of colors
    """
    for j in range(len(in_depmap_list_all[0])): 
        pval_y = []
        for i in range(len(t_types)):
            t_type = t_types[i]
            n = len(g_pool[t_type])
            N = topN
            tmp = in_depmap_list_all[i][j]
            y = tmp[topN-1]
            pvaly = hypergeom.sf(y-1, M, n, N)
            pval_y.append(pvaly)
        adj_p = list(fdrcorrection(pval_y)[1])
        pvals.append(adj_p)   # adjusted p-values

    pvals_mat = -np.log10(np.transpose(np.array(pvals)))
    df_all = pd.DataFrame(pvals_mat)
    df_all.columns = types
    make_boxplot(df_all,colors,types,topN,ftype = 'including_common_essentials_p-value',save_dir1 = save_dir1,ylabel = '-log10 (p-value)',log_scale = False)

    # ### barplot
    make_barplot(pvals,colors,save_dir1 = save_dir1,t_types = t_types,types = types)

    ### Excluding_common_essentials
    pvals = [] 
    """
    pvals is a list of multiple lists of p-values
    the length of pvals should be the same as the length of colors
    """
    for j in range(len(in_depmap_list_all_rm_ce[0])): 
        pval_y = []
        for i in range(len(t_types)):
            t_type = t_types[i]
            n = len(g_rm_ce_pool[t_type])
            N = topN
            tmp = in_depmap_list_all_rm_ce[i][j]
            y = tmp[topN-1]
            pvaly = hypergeom.sf(y-1, M, n, N)
            pval_y.append(pvaly)
        adj_p = list(fdrcorrection(pval_y)[1])
        pvals.append(adj_p)   # adjusted p-values

    pvals_mat = -np.log10(np.transpose(np.array(pvals)))
    df_all = pd.DataFrame(pvals_mat)
    df_all.columns = types
    make_boxplot(df_all,colors,types,topN,ftype = 'excluding_common_essentials_p-value',save_dir1 = save_dir1, ylabel = '-log10 (p-value)',log_scale = False)

    # ### barplot
    make_barplot(pvals,colors,save_dir1 = save_dir1,t_types = t_types,types = types,ftype = 'excluding_common_essentials')


################################# Defined functions below #################################
    
def set_size(w,h,ax=None,axis_linewidth = 0.75):
    """ w, h: width, height in inches """
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

def make_plots(in_depmap_list,t_type,freq_list,types,topN,colors,save_dir1,ftype='',
           vertical_line = [], plot_freq = True, axis_linewidth = 0.5,linewidth = 0.75,
           ylabel_left = 'Essential gene count',
           ylabel_right = 'Altered rate',
           xlabel = 'Rank of gene',
           dataset = 'depmap'): 
    """
    Length of 'types' should be equal to the length of 'in_depmap_list', and 'colors'
    't_type' is the tumor type
    'freq_list' is a single list, which correspond to the altered freq of the last item in the 'in_depmap_list'
    'vertical_line' is a list of single value
    """
    if len(types)!=len(colors):
        logging.error('Please input N colors if the number of types is N.')
        return
    if len(in_depmap_list)!=len(types):
        logging.error('Length of types should be the same as the length of in_depmap_list.')
        return
    if len(vertical_line)>1:
        logging.error('Vertical_line is a list of single value.')
        return

    plt.rcParams["figure.figsize"] = (2.5,2)
    plt.rcParams['svg.fonttype'] = 'none'
    fig,ax = plt.subplots()
    fig.tight_layout()
    ax = set_size(2,1.5,ax = ax,axis_linewidth = axis_linewidth)
    # ax.spines['right'].set_visible(True)
    for i in range(len(types)):
        if (i==len(types)-1):
            if len(vertical_line)==0:
                ax.plot(range(1,len(in_depmap_list[i])+1),in_depmap_list[i],color = colors[i],linewidth=linewidth,label = types[i])
            else:
                ax.plot(range(1,vertical_line[0]+1),in_depmap_list[i][:vertical_line[0]],color = colors[i],linewidth=linewidth,label = types[i])
        else:
            ax.plot(range(1,len(in_depmap_list[i])+1),in_depmap_list[i],color = colors[i],linewidth=linewidth,label = types[i])
    plt.legend(fontsize = 5)
    ax.set_ylabel(ylabel = ylabel_left,fontsize = 7)
    ax.set_xlabel(xlabel,fontsize = 7)
    ax2=ax.twinx()
    ax2 = set_size(2,1.5,ax = ax2)
    if plot_freq:
        ax2.spines['right'].set_visible(True)
        ax2.set_ylabel(ylabel = ylabel_right,fontsize = 7)
        ax2.tick_params(axis="y", labelsize=6)
    else:
        ax2.spines['right'].set_visible(False)
        ax2.tick_params(right=False, labelright=False)
    if len(vertical_line)!=0:
        j = vertical_line[0]
        ax2.vlines(x=j, ymin = 0, ymax = max(freq_list), color = 'gray', linestyle='dashed',linewidth=linewidth)
        if plot_freq:
            ax2.plot(range(1,j+1),freq_list[:j],color = colors[i], linestyle='dashed',linewidth=linewidth)
    else:
        if plot_freq:
            ax2.plot(range(1,len(in_depmap_list[i])+1),freq_list,color = colors[i], linestyle='dashed',linewidth=linewidth)
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    if len(ftype)==0:
        plt.title(t_type,fontsize = 6)
    else:
        plt.title(t_type+'_'+ftype,fontsize = 6)
        # plt.tight_layout()
    labels = '_'.join(types)
    fnm = os.path.join(save_dir1,f'{dataset}_saturation_{labels}_{t_type}_'+ftype+'.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(save_dir1,f'{dataset}_saturation_{labels}_{t_type}_'+ftype+'.png')
    plt.savefig(fnm,dpi = 300)
    plt.close('all')
    # plt.show()

    fig,ax = plt.subplots()
    fig.tight_layout()
    ax = set_size(2,1.5,ax = ax)
    # ax.spines['right'].set_visible(True)
    min_n = np.min([topN,min([len(k) for k in in_depmap_list])])

    for i in range(len(types)):
        if (i==len(types)-1):
            if len(vertical_line)==0 or vertical_line[0]>min_n:
                ax.plot(range(1,min_n+1),in_depmap_list[i][:min_n],color = colors[i],linewidth=linewidth,label = types[i])
            else:
                ax.plot(range(1,vertical_line[0]+1),in_depmap_list[i][:vertical_line[0]],color = colors[i],linewidth=linewidth,label = types[i])
        else:
            ax.plot(range(1,min_n+1),in_depmap_list[i][:min_n],color = colors[i],linewidth=linewidth,label = types[i])
    ax.set_ylabel(ylabel = ylabel_left,fontsize = 7)
    ax.set_xlabel(xlabel,fontsize = 7)
    plt.legend(fontsize = 5)
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax2=ax.twinx()
    ax2 = set_size(2,1.5,ax = ax2)
    if plot_freq:
        ax2.spines['right'].set_visible(True)
        ax2.set_ylabel(ylabel = ylabel_right,fontsize = 7)
        ax2.tick_params(axis="y", labelsize=6)
    else:
        ax2.spines['right'].set_visible(False)
        ax2.tick_params(right=False, labelright=False)
    if len(vertical_line)!=0:
        j = vertical_line[0]
        if j<=topN:
            ax2.vlines(x=j, ymin = 0, ymax = max(freq_list), color = 'gray', linestyle='dashed',linewidth=linewidth)
            if plot_freq:
                ax2.plot(range(1,j+1),freq_list[:j],color = colors[i], linestyle='dashed',linewidth=linewidth)
        else:
            if plot_freq:
                ax2.plot(range(1,min_n+1),freq_list[:min_n],color = colors[i], linestyle='dashed',linewidth=linewidth)
    else:
        if plot_freq:
            ax2.plot(range(1,min_n+1),freq_list[:min_n],color = colors[i], linestyle='dashed',linewidth=linewidth)
    if len(ftype)==0:
        plt.title(t_type,fontsize = 6)
    else:
        plt.title(t_type+'_'+ftype,fontsize = 6)
    # plt.tight_layout()
    fnm = os.path.join(save_dir1,f'{dataset}_saturation_{labels}_of_top_{topN}_{t_type}_'+ftype+'.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(save_dir1,f'{dataset}_saturation_{labels}_of_top_{topN}_{t_type}_'+ftype+'.png')
    plt.savefig(fnm,dpi = 300)
    plt.close('all')
    # plt.show()

def area_between_curves(in_depmap_list,topN = 100):
    """
    Compute area under curves
    'in_depmap_list' maybe a list of multiple lists or a single list.
    """
    if len(in_depmap_list)==0:
        print('no list of in depmap counts!')
        return
    if any(isinstance(el, list) for el in in_depmap_list): # if in_depmap_list is a list of lists
        min_n = np.min([topN,min([len(i) for i in in_depmap_list])])
        a = []
        for i in range(len(in_depmap_list)):
            l = in_depmap_list[i][:min_n]
            a.append(0.5*l[-1]+np.sum(l[0:(-1)]))
        return a
    else:
        min_n = np.min([topN,len(in_depmap),len(in_depmap_list)])
        a = []
        l = in_depmap_list[:min_n]
        a.append(0.5*l[-1]+np.sum(l[0:(-1)]))
        return a

def make_boxplot(df_all,colors,types,topN,ftype,save_dir1,ylabel = 'Area under the curve',log_scale = False):
    """The number of columns of df_all should be the same as the length of types"""
    if len(df_all.columns)!=len(types):
        logging.error("The number of columns of 'df_all' should be the same as the length of 'types'!")
        return
    if len(df_all.columns)!=len(colors):
        logging.error("The number of columns of 'df_all' should be the same as the length of 'colors'!")        
        return

    fig, ax = plt.subplots()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["figure.figsize"] = (1.2,1.4)
    ax = set_size(1,1.2,ax = ax)

    N = len(df_all.columns)
    my_pal = {}
    for i in range(N):
        my_pal[types[i]]=colors[i]
    sns.boxplot(data = df_all,palette=my_pal,linewidth = 0.5, flierprops={"marker": "x", "markersize":0})

    for i in range(N):
        x = df_all.iloc[:,i]
        plt.scatter([i+np.random.normal(scale = 0.01) for j in x],x,s = .2,c = 'k')

    for i in range(N-1):
        for j in range(i+1,N):
            x = df_all.iloc[:,i]
            y = df_all.iloc[:,j]
            z = ranksums(x,y)    
            tmp1 = 1.1*max(x)
            tmp2 = 1.1*max(y)
            tmp = max(tmp1,tmp2)
            plt.vlines(x=i,ymin = tmp,ymax=tmp)
            plt.vlines(x=j,ymin = tmp,ymax=tmp)
            plt.hlines(y=2*(i+j)+tmp, xmin = i, xmax = j,color = 'k',linewidth = 0.5)

            xlocation = (i+j)/2-0.2
            ylocation = 2*(i+j)+tmp
            if (z[1]<0.05)&(z[1]>=0.01):
                plt.text(x=xlocation,y=ylocation,s='*',fontsize = 'xx-small')
            elif (z[1]<0.01)&(z[1]>=0.001):
                plt.text(x=xlocation,y=ylocation,s='**',fontsize = 'xx-small')
            elif z[1]<0.001:
                plt.text(x=xlocation,y=ylocation,s='***',fontsize = 'xx-small')
            else:
                plt.text(x=xlocation,y=ylocation,s='n.s.',fontsize = 'xx-small')

    plt.title(f'Top {topN} genes',fontsize = 6)
    plt.ylabel(ylabel,fontsize = 6)
    plt.xticks(rotation = 90, fontsize = 5)
    plt.yticks(fontsize = 5)
    plt.tight_layout()
    if log_scale == True:
        ax.set_yscale('log')
    fnm = os.path.join(save_dir1,f'compare_{types[0]}_and_top_{topN}_altered_{types[1]}_{ftype}.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(save_dir1,f'compare_{types[0]}_and_top_{topN}_altered_{types[1]}_{ftype}.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()

### make barplot for adjusted p-values
def make_barplot(pvals,colors,save_dir1,t_types,types,topN = 100,width = 0.4,fig_size = (2.5,1.7),axis_linewidth = 0.75,linewidth = 0.75,ftype = 'including_common_essential'):
    """
    pvals is a list of multiple lists of p-values
    the length of pvals should be the same as the length of colors
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

    N = len(pvals)
    x = np.arange(0,N*len(pvals[0]),N)
    for i in range(N):
        tmp = pvals[i]
        plt.bar(x+(width+0.2)*i, -np.log10(tmp), width, color = colors[i])

    plt.hlines(y = -np.log10(0.05), xmin = 0, xmax = N*len(pvals[0]), color = 'gray', linestyle='dashed',
               linewidth = linewidth)
    plt.xticks(x+0.25,t_types,fontsize = 5, rotation = 90)
    plt.yticks(fontsize = 5)
    # ax.set_yscale('log')
    plt.title(f'Top {topN} genes',fontsize = 6)
    plt.ylabel('-log10 (adjusted p-value)',fontsize = 6)
    plt.xticks(rotation = 90, fontsize = 5)
    plt.yticks(fontsize = 5)
    plt.tight_layout()
    labels = '_'.join(types)
    fnm = os.path.join(save_dir1,f'compare_model_and_top_{topN}_altered_{labels}_hypergeometric_test_{ftype}_adj_pvalue_barplot.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(save_dir1,f'compare_model_and_top_{topN}_altered_{labels}_hypergeometric_test_{ftype}_adj_pvalue_barplot.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()

def make_plot_errorbar(in_depmap_list_all, topN, types, colors, save_dir1, ftype='',
               axis_linewidth = 0.5,linewidth = 0.75,
               ylabel_left = 'Count of genetic dependencies',
               xlabel = 'Rank of gene'):
    """
    'in_depmap_list_all' is a list of N lists, N is the number of tumor types. 
    Each of the N list has m elements, m is the same as the length of 'colors' and 'types'
    """
    plt.rcParams["figure.figsize"] = (2,1.7)
    plt.rcParams['svg.fonttype'] = 'none'
    fig,ax = plt.subplots()
    fig.tight_layout()
    ax = set_size(1.7,1.4,ax = ax,axis_linewidth = axis_linewidth)

    for i in range(len(colors)):
        in_depmap_wes = []
        for item in range(len(in_depmap_list_all)):
            in_depmap_list = in_depmap_list_all[item]
            tmp = in_depmap_list[i]
            in_depmap_wes.append(tmp)
        df_s = np.array(in_depmap_wes)
        df_mean = np.mean(df_s,axis = 0)[:topN]
        df_std = np.std(df_s,axis = 0)[:topN]
        plt.plot(range(1,len(df_mean)+1),df_mean,color = colors[i])
        plt.errorbar(range(1,len(df_mean)+1), df_mean, yerr = df_std,elinewidth = 0.1,color = colors[i],alpha = 0.5,label = types[i])
            # plt.tight_layout()
    plt.legend(loc = 'upper left',fontsize = 6)
    # ax.legend(g,labels = types,fontsize = 5)
    ax.set_ylabel(ylabel = ylabel_left,fontsize = 7)
    ax.set_xlabel(xlabel,fontsize = 7)
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=6)

    tmp = '_'.join(types)
    fnm = os.path.join(save_dir1,f'plot_error_bar_of_top_{topN}_{tmp}_{ftype}.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(save_dir1,f'plot_error_bar_of_top_{topN}_{tmp}_{ftype}.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
