####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes perform GDSC validation
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle

####################################################################################
### load the results from CanSAR
def GDSC_validation(dir0,save_dir,t_types,topN = 100):
          
    fig_dir = os.path.join(save_dir,'figures')
    mat_dir = os.path.join(save_dir,'output_data_matrix')
    dir3 = os.path.join(dir0,'databases/GDSC/select_tumor_types')
    df0 = load_top_ranked_with_essential_info(mat_dir,topN)
    ### for MVP model topN genes
    output = {}
    gene_s = []
    t_type_s = []
    score_s = []
    cellline_s = []
    drug_s = []
    gene_nm_dict = {
        'PIK3CA':'PI3Kalpha',
        'PIK3CB':'PI3Kbeta',
        'PIK3CD':'PI3Kdelta',
        'PIK3CG':['PI3Kgamma','PK3CG'],
        'MAP2K1':'MEK1',
        'MAP2K2':'MEK2',
        'MAPK1':'ERK2',
        'MAPK3':'ERK1',
        'CREBBP':'CBP',
        'MTOR':['mTOR','MTOR'],
        'MAPK11':'p38beta',
        'MAPK14':'p38alpha',
        'HSP90AA1':'HSP90',
        'HSP90AA2':'HSP90',
        'HSP90AB1':'HSP90',
        'PTK2':'FAK',
        'ERK5':'MAPK7',
        'MAP2K5':'MEK5',
        'PDPK1':'PDK1 (PDPK1)',
        'IDH1':'IDH1 (R132H)',
        'KRAS':'KRAS (G12C)'
    }
    
    for t_type in t_types:
        fnm = t_type+'_GDSC1.csv'
        if fnm in os.listdir(dir3):
            output[t_type] = []
            dfa = pd.read_csv(os.path.join(dir3,t_type+'_GDSC1.csv'))
            dfb = pd.read_csv(os.path.join(dir3,t_type+'_GDSC2.csv'))
            dfab = pd.concat([dfa,dfb],axis = 0)           
            genes = df0.loc[df0.cancer_type==t_type,'top_genes'].values.tolist()            
            for i in range(len(dfab.PUTATIVE_TARGET.values)):
                if type(dfab.PUTATIVE_TARGET.values[i])==str:
                    tmp = dfab.PUTATIVE_TARGET.values[i].split(', ')
                    if len(tmp)<=4:
                        for g in genes:
                            g1 = g
                            if g in gene_nm_dict:
                                g1 = gene_nm_dict[g]
                            if type(g1)==list:
                                for element in g1:
                                    if element in tmp:
                                        gene_s.append(g)
                                        score_s.append(dfab.iloc[i,11])
                                        t_type_s.append(t_type)
                                        cellline_s.append(dfab.iloc[i,1])
                                        drug_s.append(dfab.iloc[i,3])
                                        output[t_type].append(dfab.iloc[i,1]+':'
                                                              +dfab.iloc[i,3]+':'
                                                             +dfab.iloc[i,4]+':'
                                                             +str(dfab.iloc[i,11]))
                                        continue
                            else:
                                if g1 in tmp:
                                    gene_s.append(g)
                                    score_s.append(dfab.iloc[i,11])
                                    t_type_s.append(t_type)
                                    cellline_s.append(dfab.iloc[i,1])
                                    drug_s.append(dfab.iloc[i,3])
                                    output[t_type].append(dfab.iloc[i,1]+':'
                                                        +dfab.iloc[i,3]+':'
                                                        +dfab.iloc[i,4]+':'
                                                        +str(dfab.iloc[i,11]))
    target_score = pd.DataFrame({'gene':gene_s,'z_score_IC50':score_s,
                                'tumor_type':t_type_s,'cell_line':cellline_s,
                                'drug':drug_s})
    
    ### save the target score for the top ranked genes
    fnm = os.path.join(mat_dir,f'z_score_IC50_top_{topN}_ranked_genes_by_MVP.csv')
    target_score.to_csv(fnm,index = False)
    
    ### protein coding genes
    dir2 = 'databases/TCGA_GTEx'
    fnm = os.path.join(dir0,dir2,'TcgaTargetGTEX_gene_biotype_by_biomaRt.csv')
    fx = pd.read_csv(fnm)
    pcg = fx.loc[fx.gene_biotype=='protein_coding','hgnc_symbol'].values.tolist()
    print(f'number of protein coding genes {len(pcg)}')
    ### for top altered N genes
    output1 = {}
    gene_s = []
    t_type_s = []
    score_s = []
    cellline_s = []
    drug_s = []
    for t_type in t_types:
        fnm = t_type+'_GDSC1.csv'    
        if fnm in os.listdir(dir3):
            fq_dir = os.path.join(save_dir, t_type, 'TCGA_'+t_type+'_wes_freq_long.csv')
            df_q = pd.read_csv(fq_dir)
            df_q = df_q.loc[df_q.Gene.isin(pcg),]
            df_q = df_q.reset_index(drop = True)
            genes = df_q.iloc[:topN,0].values.tolist()
            output1[t_type] = []
            dfa = pd.read_csv(os.path.join(dir3,t_type+'_GDSC1.csv'))
            dfb = pd.read_csv(os.path.join(dir3,t_type+'_GDSC2.csv'))
            dfab = pd.concat([dfa,dfb],axis = 0)                            
            for i in range(len(dfab.PUTATIVE_TARGET.values)):
                if type(dfab.PUTATIVE_TARGET.values[i])==str:
                    tmp = dfab.PUTATIVE_TARGET.values[i].split(', ')
                    if len(tmp)<=4:
                        for g in genes:
                            g1 = g
                            if g in gene_nm_dict:
                                g1 = gene_nm_dict[g]
                            if type(g1)==list:
                                for element in g1:
                                    if element in tmp:
                                        gene_s.append(g)
                                        score_s.append(dfab.iloc[i,11])
                                        t_type_s.append(t_type)
                                        cellline_s.append(dfab.iloc[i,1])
                                        drug_s.append(dfab.iloc[i,3])
                                        output1[t_type].append(dfab.iloc[i,1]+':'
                                                              +dfab.iloc[i,3]+':'
                                                             +dfab.iloc[i,4]+':'
                                                             +str(dfab.iloc[i,11]))
                                        continue
                            else:
                                if g1 in tmp:
                                    gene_s.append(g)
                                    score_s.append(dfab.iloc[i,11])
                                    t_type_s.append(t_type)
                                    cellline_s.append(dfab.iloc[i,1])
                                    drug_s.append(dfab.iloc[i,3])
                                    output1[t_type].append(dfab.iloc[i,1]+':'
                                                          +dfab.iloc[i,3]+':'
                                                         +dfab.iloc[i,4]+':'
                                                         +str(dfab.iloc[i,11]))
    target_score1 = pd.DataFrame({'gene':gene_s,'z_score_IC50':score_s,
                                'tumor_type':t_type_s,'cell_line':cellline_s,
                                'drug':drug_s})
    
    ### save the target score for the top altered genes
    fnm = os.path.join(mat_dir,f'z_score_IC50_top_{topN}_altered_genes.csv')
    target_score1.to_csv(fnm,index = False)
    
    ### plot heatmap of IC50 of gdsc targets
    heatmap_top_ranked_by_MVP(t_types,dir3,fig_dir,output,topN)
    heatmap_top_altered(t_types,dir3,fig_dir,output1,topN)
    compare_GDSC_top_ranked_vs_top_altered(dir3,output,output1,t_types,topN,fig_dir)
    compare_IC50(mat_dir,topN,fig_dir)
    
    ### GDSC targets among all genes in the disease network
    fnm = os.path.join(mat_dir,'rank_model_parameters_all_cancers.csv')
    df_rnk = pd.read_csv(fnm, index_col = 0)
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    select_md = df_rnk.index[0]
    
    output = {}
    gene_s = []
    t_type_s = []
    score_s = []
    cellline_s = []
    drug_s = []

    for t_type in t_types:
        fnm = t_type+'_GDSC1.csv'
        if fnm in os.listdir(dir3):
            output[t_type] = []
            dfa = pd.read_csv(os.path.join(dir3,t_type+'_GDSC1.csv'))
            dfb = pd.read_csv(os.path.join(dir3,t_type+'_GDSC2.csv'))
            dfab = pd.concat([dfa,dfb],axis = 0)
            print(t_type)
            print('---------------------------------------------')
            ### MVP model ###
            md_dir = os.path.join(save_dir, t_type, 'markov chain models')
            s_md = 'markov_output_wes_'+select_md+'_022123.csv'
            df = pd.read_csv(os.path.join(md_dir,s_md))
            genes = df.genes.values.tolist()

            for i in range(len(dfab.PUTATIVE_TARGET.values)):
                if type(dfab.PUTATIVE_TARGET.values[i])==str:
                    tmp = dfab.PUTATIVE_TARGET.values[i].split(', ')
                    if len(tmp)<=4:
                        for g in genes:
                            g1 = g
                            if g in gene_nm_dict:
                                g1 = gene_nm_dict[g]
                            if type(g1)==list:
                                for element in g1:
                                    if element in tmp:
                                        gene_s.append(g)
                                        score_s.append(dfab.iloc[i,11])
                                        t_type_s.append(t_type)
                                        cellline_s.append(dfab.iloc[i,1])
                                        drug_s.append(dfab.iloc[i,3])
                                        output[t_type].append(dfab.iloc[i,1]+':'
                                                              +dfab.iloc[i,3]+':'
                                                             +dfab.iloc[i,4]+':'
                                                             +str(dfab.iloc[i,11]))
                                        continue
                            else:
                                if g1 in tmp:
                                    gene_s.append(g)
                                    score_s.append(dfab.iloc[i,11])
                                    t_type_s.append(t_type)
                                    cellline_s.append(dfab.iloc[i,1])
                                    drug_s.append(dfab.iloc[i,3])
                                    output[t_type].append(dfab.iloc[i,1]+':'
                                                          +dfab.iloc[i,3]+':'
                                                         +dfab.iloc[i,4]+':'
                                                         +str(dfab.iloc[i,11]))
    target_score = pd.DataFrame({'gene':gene_s,'z_score_LN_IC50':score_s,
                                  'tumor_type':t_type_s,'cell_line':cellline_s,
                                 'drug':drug_s})

    ### save the target score for the top ranked genes
    fnm = os.path.join(mat_dir,f'z_score_LN_IC50_all_network_genes.csv')
    target_score.to_csv(fnm,index = False)
    
    ### targets_in_GDSC_all_network_genes
    tumor_types = []
    num_gene_in_gdsc = []
    num_gene_low_z_score = []
    genes_with_low_z_score = []

    z_df = target_score.copy()
    t_types1 = sorted(set(z_df['tumor_type']))
    for t_type in t_types1:
        # top ranked by MVP
        tumor_types.append(t_type)
        p1 = set(z_df.loc[z_df.tumor_type==t_type,'gene'].values.tolist())  
        num_gene_in_gdsc.append(len(p1))
        sub = z_df.loc[(z_df.tumor_type==t_type)&(z_df.z_score_LN_IC50<-1.5),]
        p2 = set(sub.loc[sub.tumor_type==t_type,'gene'].values.tolist())   
        genes_with_low_z_score.append(','.join(sorted(p2)))
        num_gene_low_z_score.append(len(p2))     

    summary_table = pd.DataFrame({'Tumor type':tumor_types,
                                  f'GDSC target count':num_gene_in_gdsc,
                                 f'GDSC target count with low z-scored LN_IC50':num_gene_low_z_score,
                                 f'GDSC targets with low z-scored LN_IC50':genes_with_low_z_score})
    ### save dataframe
    fnm = os.path.join(mat_dir,f'compare_targets_in_GDSC_all_network_genes.csv')
    summary_table.to_csv(fnm, index = False)
    

def heatmap_top_ranked_by_MVP(t_types,dir3,fig_dir,output,topN=20,fig_size=(16,16)):
    ### topN ranked by MVP
    for t_type in t_types:
        fnm = t_type+'_GDSC1.csv'
        if fnm in os.listdir(dir3):
            tmp = output[t_type]
            cell_lines = []
            drug_and_target = []
            z_IC50 = []
            for i in range(len(tmp)):
                tmp1 = tmp[i].split(':')
                cell_lines.append(tmp1[0])
                drug_and_target.append(tmp1[1]+':'+tmp1[2])
                z_IC50.append(float(tmp1[3]))

            zIC50_df = pd.DataFrame({'cell line of '+t_type:cell_lines, 
                                    'drug and targets':drug_and_target,
                                    'z-score of IC50':z_IC50})
            zIC50_df = zIC50_df.groupby(['cell line of '+t_type,'drug and targets']).mean().reset_index()
            dfzx = zIC50_df.pivot(index = 'cell line of '+t_type,columns = 'drug and targets',values = 'z-score of IC50')            
            if dfzx.shape[0]>0:
                plt.rcParams["figure.figsize"] = fig_size
                dfzx = dfzx.fillna(0)
                g = sns.clustermap(dfzx.T,center = 0,cmap = 'PiYG', cbar_pos=(0.08, 1, .01, .01),
                            cbar_kws={'shrink': 0.2},linewidth = 0.5,linecolor = 'gray',vmin=-3, vmax=3)
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 8, rotation = 90)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 9)
                g.ax_heatmap.set_xlabel('cell line of '+t_type, fontsize=18)
                g.ax_heatmap.set_ylabel(f'drug and targets (top {topN} ranked genes by MVP)', fontsize = 18)                
                title = f'sensitivity in z-score'
                plt.title(title, fontsize = 14)
                fnm = os.path.join(fig_dir,t_type+'_WES_SV_as_seed_sensitivity in z-score_top'+str(topN)+'_heatmap_z-score_IC50.png')
                plt.tight_layout()
                plt.savefig(fnm,dpi = 300, bbox_inches='tight')
                plt.show()
                # plt.close('all')

def heatmap_top_altered(t_types,dir3,fig_dir,output1,topN=20):
    ### topN altered
    for t_type in t_types:
        fnm = t_type+'_GDSC1.csv'
        if fnm in os.listdir(dir3):
            tmp = output1[t_type]
            cell_lines = []
            drug_and_target = []
            z_IC50 = []
            for i in range(len(tmp)):
                tmp1 = tmp[i].split(':')
                cell_lines.append(tmp1[0])
                drug_and_target.append(tmp1[1]+':'+tmp1[2])
                z_IC50.append(float(tmp1[3]))
            zIC50_df = pd.DataFrame({'cell line of '+t_type:cell_lines, 
                                    'drug and targets':drug_and_target,
                                    'z-score of IC50':z_IC50})
            zIC50_df = zIC50_df.groupby(['cell line of '+t_type,'drug and targets']).mean().reset_index()
            dfzx = zIC50_df.pivot(index = 'cell line of '+t_type,columns = 'drug and targets',values = 'z-score of IC50')            
            if dfzx.shape[1]>0:
                if dfzx.shape[1]>1:
                    plt.rcParams["figure.figsize"] = (30,15)
                    dfzx = dfzx.fillna(0)
                    g = sns.clustermap(dfzx.T,center = 0,cmap = 'PiYG', cbar_pos=(0.08, 1, .01, .01),
                                cbar_kws={'shrink': 0.2},linewidth = 0.5,linecolor = 'gray',vmin=-3, vmax=3)
                    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 8, rotation = 90)
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 9)
                    g.ax_heatmap.set_xlabel('cell line of '+t_type, fontsize=18)
                    g.ax_heatmap.set_ylabel(f'drug and targets (top {topN} altered genes)', fontsize = 18)
                else:
                    plt.rcParams["figure.figsize"] = (10,8)
                    dfzx = dfzx.fillna(0)
                    g = sns.heatmap(dfzx.T,center = 0,cmap = 'PiYG', 
                                cbar_kws={'shrink': 0.2},linewidth = 0.5,linecolor = 'gray',vmin=-3, vmax=3)
                    g.set_xticklabels(g.get_xticklabels(), fontsize = 10, rotation = 90)
                    g.set_yticklabels(g.get_yticklabels(), fontsize = 10)
                    g.set_xlabel('cell line of '+t_type, fontsize=18)
                    g.set_ylabel('drug and targets (top 20 altered genes)', fontsize = 18)               
                title = f'sensitivity in z-score'
                plt.title(title, fontsize = 14)
                fnm = os.path.join(fig_dir,t_type+'_WES_SV_as_seed_sensitivity in z-score_top'+str(topN)+'_altered_heatmap_z-score_IC50.png')
                plt.tight_layout()
                plt.savefig(fnm,dpi = 256, bbox_inches='tight')
                plt.show()
                # plt.close('all')

def load_top_ranked_with_essential_info(mat_dir,topN = 20):
    file_dir = os.path.join(mat_dir,f'top_{topN}_genes_MVP_with_depmap_essential_info_including common essentials.csv')
    df = pd.read_csv(file_dir)
    return df
    
def make_violinplot(x,y,t_type,topN,fig_dir):
    data = {f'top {topN} ranked by MVP':x,f'top {topN} altered':y}
    data_df = pd.DataFrame([data[f'top {topN} ranked by MVP'], data[f'top {topN} altered']]).transpose()

    # Label the columns of the DataFrame
    data_df = data_df.set_axis([f'top {topN} ranked by MVP',f'top {topN} altered'], axis=1)

    # Violin plot
    plt.rcParams["figure.figsize"] = (8,6)
    plt.title(t_type, fontsize = 24)
    
    colors = ['#00A087FF', '#B09C85FF']
    my_pal = {f'top {topN} ranked by MVP': colors[0], f'top {topN} altered':colors[1]}
    sns.violinplot(data=data_df,
                   inner="quartile",linewidth = 2,color = 'w',scale = 'count')
    sns.swarmplot(data=data_df, size = 2,palette = my_pal)
    fnm = os.path.join(fig_dir,t_type+'_WES and SV as seed_GDSC_violinplot_comparison.png')
    plt.ylabel('z-score of LN_IC50',fontsize = 22)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 18)
    plt.savefig(fnm,dpi = 256)
    plt.close('all')
    
def compare_GDSC_top_ranked_vs_top_altered(dir3,output,output1,t_types,topN,fig_dir):
    ### compare GDSC results for topN ranked by MVP vs. topN altered
    for t_type in t_types:
        fnm = t_type+'_GDSC1.csv'
        if fnm in os.listdir(dir3):
        # topN ranked by MVP
            tmp = output[t_type]
            cell_lines = []
            drug_and_target = []
            z_IC50 = []

            for i in range(len(tmp)):
                tmp1 = tmp[i].split(':')
                cell_lines.append(tmp1[0])
                drug_and_target.append(tmp1[1]+':'+tmp1[2])
                z_IC50.append(float(tmp1[3]))

            zIC50_df = pd.DataFrame({'cell line of '+t_type:cell_lines, 
                                    'drug and targets':drug_and_target,
                                    'z-score of IC50':z_IC50})
            zIC50_df = zIC50_df.groupby(['cell line of '+t_type,'drug and targets']).mean().reset_index()
            dfzx = zIC50_df.pivot(index = 'cell line of '+t_type,columns = 'drug and targets',values = 'z-score of IC50')

         # topN altered  
            tmp = output1[t_type]
            cell_lines = []
            drug_and_target = []
            z_IC50 = []

            for i in range(len(tmp)):
                tmp1 = tmp[i].split(':')
                cell_lines.append(tmp1[0])
                drug_and_target.append(tmp1[1]+':'+tmp1[2])
                z_IC50.append(float(tmp1[3]))

            zIC50_df = pd.DataFrame({'cell line of '+t_type:cell_lines, 
                                    'drug and targets':drug_and_target,
                                    'z-score of IC50':z_IC50})
            zIC50_df = zIC50_df.groupby(['cell line of '+t_type,'drug and targets']).mean().reset_index()
            dfzx1 = zIC50_df.pivot(index = 'cell line of '+t_type,columns = 'drug and targets',values = 'z-score of IC50')

            x = dfzx.to_numpy().flatten()
            y = dfzx1.to_numpy().flatten()
            make_violinplot(x,y,t_type,topN,fig_dir)

def compare_IC50(mat_dir,topN,fig_dir):
    ### compare the z-score of IC50 top ranked and top altered genes
    fnm = os.path.join(mat_dir,f'z_score_IC50_top_{topN}_ranked_genes_by_MVP.csv')
    fnm1 = os.path.join(mat_dir,f'z_score_IC50_top_{topN}_altered_genes.csv')
    z_df = pd.read_csv(fnm)
    z_df1 = pd.read_csv(fnm1)
    t_types1 = sorted(set(z_df['tumor_type']))
    tumor_types = []
    num_gene_in_gdsc = []
    types = []
    num_gene_low_z_score = []
    genes_with_low_z_score = []
    percent_gene_low_score = []

    for t_type in t_types1:
        # top ranked by MVP 
        tumor_types.append(t_type)
        p1 = set(z_df.loc[z_df.tumor_type==t_type,'gene'].values.tolist())  
        num_gene_in_gdsc.append(len(p1))
        sub = z_df.loc[(z_df.tumor_type==t_type)&(z_df.z_score_IC50<-1.5),]
        types.append(f'MVP') 
        p2 = set(sub.loc[sub.tumor_type==t_type,'gene'].values.tolist())   
        genes_with_low_z_score.append(','.join(sorted(p2)))
        num_gene_low_z_score.append(len(p2))
        percent_gene_low_score.append(100*len(p2)/topN)
        
        # top altered
        if t_type not in z_df1.tumor_type.values:        
            num_gene_in_gdsc.append(0)
            num_gene_low_z_score.append(0)
            percent_gene_low_score.append(0)
            genes_with_low_z_score.append('')
        else:
            p1 = set(z_df1.loc[z_df1.tumor_type==t_type,'gene'].values.tolist())
            num_gene_in_gdsc.append(len(p1))
            sub1 = z_df1.loc[(z_df1.tumor_type==t_type)&(z_df1.z_score_IC50<-1.5),]
            p2 = set(sub1.loc[sub1.tumor_type==t_type,'gene'].values.tolist())
            num_gene_low_z_score.append(len(p2))
            genes_with_low_z_score.append(','.join(sorted(p2)))
            percent_gene_low_score.append(100*len(p2)/topN)
            
        tumor_types.append(t_type)
        types.append(f'Top altered')
    summary_table = pd.DataFrame({'Tumor type':tumor_types,'Method':types,
                                f'GDSC target count among top {topN} ranked':num_gene_in_gdsc,
                                f'GDSC target count with low z-scored LN_IC50 among top {topN} ranked':num_gene_low_z_score,
                                f'Percent of top {topN} genes with low z-scored LN_IC50 (%)':percent_gene_low_score,
                                f'GDSC targets among top {topN} with low z-scored LN_IC50':genes_with_low_z_score})
    plt.rcParams["figure.figsize"] = (9,6)
    clrs = ['#00A087FF','#B09C85FF']
    fig, ax = plt.subplots()
    sns.barplot(data=summary_table, x="Tumor type", 
                y=f'GDSC target count with low z-scored LN_IC50 among top {topN} ranked', 
                hue="Method",palette=clrs,alpha = 0.9)
    plt.xlabel('tumor type',fontsize = 14)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels)
    fnm = os.path.join(fig_dir,f'GDSC target count with low z-scored LN_IC50 among top {topN} ranked barplot.png')
    plt.savefig(fnm,dpi = 256)

    ### save dataframe
    fnm = os.path.join(mat_dir,f'compare_targets_in_GDSC_top_ranked_MVP_and_top_{topN}_altered.csv')
    summary_table.to_csv(fnm, index = False)
