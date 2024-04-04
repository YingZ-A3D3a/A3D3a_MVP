####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes set initial score of driver/onco genes to zero
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from matplotlib import rc
import textalloc as ta
### build networkx graph from adjacency matrix
import PGM_PY.adj_mat_interactome as ami
import PGM_PY.graphical_models as gms # run graphical model
import PGM_PY.graph_builder as gb # largest connected graph
import PGM_PY.markov_model_dataframe as mmd # degree of genes in the graph

####################################################################################
def set_initial_to_zero(dir0,save_dir,topN = 60,
                           t_names = ['Bladder Cancer','Breast Cancer',
          'Colorectal Cancer','Head and Neck Squamous Cell Carcinoma','Renal Cell Carcinoma',
          'Glioma','Non-Small Cell Lung Cancer','Ovarian Cancer',
          'Prostate Cancer','Thyroid Cancer',['Endometrial Cancer','Uterine Serous Carcinoma','Uterine Papillary Serous Carcinoma'],
          'Melanoma'],ct_types = ['BLCA','BRCA',
            'COADREAD','HNSC','KIRC',
            'LGG','LUAD','OV',
            'PRAD','THCA','UCEC',
            'SKCM']):  # 't_names': tumor types in oncoKB, 'ct_types':corresponding TCGA types

    ### load Markov model parameters with best performance
    fig_dir = os.path.join(save_dir,'figures')
    mat_dir = os.path.join(save_dir,'output_data_matrix')

    gm = ami.adj_mat()
    ### load oncoKB results
    dir2 = os.path.join(dir0,'databases/oncoKB')
    fnm = os.path.join(dir2,'oncokb_biomarker_drug_associations.tsv')
    df = pd.read_csv(fnm,sep = '\t')

    tmp = df['Cancer Types']
    t_type = 'Melanoma'
    idx = []
    for i in range(len(tmp)):
        item = tmp[i]
        tmp1 = item.split('/')
        tmp2 = item.split(', ')
        if (t_type in tmp1) or (t_type in tmp2):
            idx.append(i)
    
    ### save the genes level 1-3 from oncokb
    ### get the oncogene from oncokb
    onco = df['Cancer Types']
    df_all = pd.DataFrame()

    for i in range(len(ct_types)):
        t_type = ct_types[i]
        t_nm = t_names[i]
        
        idx = []
        for j in range(len(onco)):
            item = onco[j]
            tmp1 = item.split('/')
            tmp2 = item.split(', ')
            if type(t_nm) is list:
                for t_nm_k in t_nm:
                    if (t_nm_k in tmp1) or (t_nm_k in tmp2):
                        idx.append(j)
            else:
                if (t_nm in tmp1) or (t_nm in tmp2):
                    idx.append(j)
        df1 = df.iloc[idx,]
        df1 = df1.loc[df1.Level.isin(['1','2','3']),]
        df_all = pd.concat([df_all,df1],axis = 0)
    df_all.reset_index(drop = True, inplace = True)
    fnm = os.path.join(dir2, 'oncogenes_level_1_to_3.csv')
    df_all.to_csv(fnm, index = False)
    
    ### add CanSAR curated list    
    target_dir = os.path.join(dir0,'databases/Cansar disease interventions target and type lists Bissan')
    tumor_types = ['BLCA','BRCA','COADREAD','HNSC',
                'KIRC','LGG','LUAD','OV','PRAD','SKCM','THCA','UCEC']
    types_in_cansar = ['Bladder Cancer','Breast Cancer','Bowel Cancer','Head and Neck Cancer',
                    'Kidney Cancer','Glioma Tumour','NSCLC Cancer','Ovarian Cancer','Prostate Cancer','Melanoma Cancer','Thyroid Cancer','Endometrial Cancer']
    
    g_list = {}
    for i in range(len(types_in_cansar)):
        tmp = types_in_cansar[i]
        fnm = os.path.join(target_dir,'Disease-'+tmp+'-Interventions.csv')
        tmp0 = pd.read_csv(fnm)
        tmp0 = tmp0.loc[tmp0.type.isin(['cytotoxic','targetted']),]
        gene_list = tmp0['targets'].values.tolist()
        g_list[tumor_types[i]] = gene_list

    ### load best Markov model, set the initial score of each of topN genes to zero and re-run the Markov model
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    select_md = df_rnk.index[0]
    tmp = select_md.split('_')
    Wm = float(tmp[3])
    alpha = float(tmp[5])

    ### get the oncogene from oncokb
    onco = df['Cancer Types']        
    for i in range(len(ct_types)):
        t_type = ct_types[i]
        t_nm = t_names[i]
        print(t_type)
        print('---------------------------------------------')
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df_mc = pd.read_csv(os.path.join(md_dir,s_md))
        g_all = df_mc['genes'].values
        g_filter = set(g_all).intersection(gm.index)
        s_list = g_filter
        score_ini = {}
        for item in s_list:
            score_ini[item] = df_mc.loc[df_mc.genes==item,'initial_score'].values[0]
            
        idx = []
        for j in range(len(onco)):
            item = onco[j]
            tmp1 = item.split('/')
            tmp2 = item.split(', ')
            if type(t_nm) is list:
                for t_nm_k in t_nm:
                    if (t_nm_k in tmp1) or (t_nm_k in tmp2):
                        idx.append(j)
            else:
                if (t_nm in tmp1) or (t_nm in tmp2):
                    idx.append(j)
        df1 = df.iloc[idx,]
        top_genes = list(set(df1.loc[df1.Level.isin(['1','2','3']),'Gene'].values.tolist()))
        top_genes = [i for i in top_genes if i!='Other Biomarkers']
        if t_type in g_list:
            top_genes_all = top_genes.copy()
            top_genes_all.extend(g_list[t_type]) # add manually curated list
            top_genes = set(top_genes_all)
        top_genes = [i for i in top_genes if i in df_mc.genes.values[:topN]]
        
        ### rank matrix stores 
        # (1) the initial score
        # (2) MC model pre-rank
        # (3) MC model post-rank (after set the initial score of a gene to zero)
        # (4) Disease network size
        # (5) degree centrality of a gene in the disease network
        
        rank_mat = pd.DataFrame(np.zeros([len(top_genes),5]))
        rank_mat.index = top_genes
        rank_mat.columns = ['Initial score','MVP pre-rank','MVP post-rank','Network size','Normalized degree']
        rank_mat.loc[:,'Network size'] = [df_mc.shape[0] for j in range(len(top_genes))]
            
        for i in range(len(top_genes)):
            geneoe = top_genes[i]
            # set initial score to zero
            score_ini1 = score_ini.copy()
            score_ini1[geneoe] = 0
            final_prob_markov0 = gms.run_pipeline_unequal(gm,geneoe,s_list,score_ini1,alpha,Wm,modelx='Markov')
            final_prob_markov0 = mmd.info_markov(final_prob_markov0,s_list,gm)
            rank_mat.loc[geneoe,'Initial score'] = df_mc.loc[df_mc.genes==geneoe,'initial_score'].values[0]
            rank_mat.loc[geneoe,'MVP pre-rank'] = df_mc.loc[df_mc.genes==geneoe,'final_rank'].values[0]
            rank_mat.loc[geneoe,'MVP post-rank'] = final_prob_markov0.loc[final_prob_markov0.genes==geneoe,'final_rank'].values[0]
            rank_mat.loc[geneoe,'Normalized degree'] = df_mc.loc[df_mc.genes==geneoe,'degree_in_disease'].values[0]

        fnm = os.path.join(mat_dir,f'{t_type}_set_initial_score_of_each_of_oncogenes_among_rank_of_top_{topN}_to_zero_MC_post_rank.csv')
        rank_mat.to_csv(fnm)
        
    ### make a single csv file
    df_comb = pd.DataFrame()
    for i in range(len(ct_types)):
        t_type = ct_types[i]
        t_nm = t_names[i]
        fnm = os.path.join(mat_dir,f'{t_type}_set_initial_score_of_each_of_oncogenes_among_rank_of_top_{topN}_to_zero_MC_post_rank.csv')
        tmp = pd.read_csv(fnm)
        tmp = tmp.rename(columns = {'Unnamed: 0':'Oncogene'})
        tmp['Tumor type'] = [t_type for i in range(tmp.shape[0])]
        cols = tmp.columns.tolist()
        cols = [cols[-1]]+cols[:-1]
        tmp = tmp[cols]
        df_comb = pd.concat([df_comb,tmp],axis = 0)
    df_comb.reset_index(inplace = True, drop = True)
    df_comb['Initial score'] = np.round(df_comb['Initial score'],2)
    df_comb['Normalized degree'] = np.round(df_comb['Normalized degree'],2)
    fnm = os.path.join(mat_dir,f'all_cancers_set_initial_score_of_each_of_oncogenes_among_rank_of_top_{topN}_to_zero_MC_post_rank.csv')
    df_comb.to_csv(fnm)

    ### visualize the rank change
    deltas = []
    gs_all = []
    delta_type = {}
    gs_type = {}
    z_type = {}
    post_rank = {}

    for t_type in ct_types:
        fnm = os.path.join(mat_dir,f'{t_type}_set_initial_score_of_each_of_oncogenes_among_rank_of_top_{topN}_to_zero_MC_post_rank.csv')
        df1 = pd.read_csv(fnm, index_col = 0)
        x = df1['MVP pre-rank'].values.tolist()
        y = df1['MVP post-rank'].values.tolist()
        delta = (df1['MVP post-rank']-df1['MVP pre-rank']).values.tolist()
        d = df1['Normalized degree'].values.tolist()
        z = df1['Initial score'].values.tolist()
        gs = [i for i in df1.index.tolist()]
        z_type[t_type] = z
        delta_type[t_type] = delta
        post_rank[t_type] = y
        gs_type[t_type] = gs
        deltas.extend(delta)
        gs_all.extend(gs)

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rcParams["figure.figsize"] = (4.5,4)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)

    for i in range(len(ct_types)):
        t_type = ct_types[i]
        x = [j+1 for j in delta_type[t_type]]
        y = [2*i+np.random.normal(scale = 0.3) for k in delta_type[t_type]]
        texts = gs_type[t_type]
        sizes = []
        for j in post_rank[t_type]:
            if j<=60:
                sizes.append(15)
            else:
                sizes.append(6)
        plt.scatter(x,y,c=z_type[t_type], edgecolors='black',linewidths=.2, alpha = 1, s = sizes, vmin = 0, vmax = 0.8, cmap='gist_earth_r')
        # Add labels to markers
        labels = [ax.text(np.round(x[j],2), np.round(y[j],2), label, fontsize = 5) for j, label in enumerate(texts)]
        # labels = [ax.text(x[i], y[i], label, ha='center', va='bottom', fontsize = 5) for i, label in enumerate(texts)]
        
        # Adjust the labels to avoid collisions
        adjust_text(labels,arrowprops=dict(arrowstyle='-', color='red', alpha=0), 
                    force_text = (1, 1))
    plt.yticks([2*i for i in range(len(ct_types))],ct_types,fontsize = 6)
    plt.xticks([i for i in np.arange(0,600,50)],fontsize = 6)

    plt.xlabel('Rank change after setting initial score to 0',fontsize = 6)
    cbar = plt.colorbar(label = 'Initial score')
    cbar.ax.tick_params(labelsize=5) # Set the font size of the color bar
    cbar.set_label(label = 'Initial score', size=5) # Set the font size of the color bar

    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'set initial to zero for oncogenes among {topN} rank change.svg')
    plt.savefig(fnm,dpi = 300)
    fnm = os.path.join(fig_dir,f'set initial to zero for oncogenes among {topN} rank change.png')
    plt.savefig(fnm,dpi = 300)
    plt.show()
