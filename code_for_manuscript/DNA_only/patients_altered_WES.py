####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes find seed genes from WES data
# Input: TCGA tumor types of interest
####################################################################################

import os
import pandas as pd
import numpy as np
from scipy.stats import gmean

####################################################################################

def seed_genes(dir0,save_dir,t_types):
    
    """
    Compute average rank (based on MutSig2CV and altered rate) for mutation data
    Save the altered freq of seed genes (WES: Mutation, CNV and SV)
    """


    for t_type in t_types:
        if t_type!='SKCM':
            dir1 = os.path.join(dir0,'TCGA_'+t_type,'WES postprocessed data from cbioportal')
            file = os.path.join(dir1,'TCGA_'+t_type+'_WES_Primary_gene_summary.csv')  # CNA and mutation data
        else:
            dir1 = os.path.join(dir0,'TCGA_'+t_type,'WES postprocessed data from cbioportal','met')
            file = os.path.join(dir1,'TCGA_'+t_type+'_WES_met_gene_summary.csv')  # CNA and mutation data
        
        file_mutsig = os.path.join(dir1,'sig_genes.txt')
        sv = os.path.join(dir1,'Structural_Variant_Genes.txt')
        
        ### WES data patients with altered gene
        wes_df = pd.read_csv(file, sep = ',', index_col = 0)
        mut_sig = pd.read_csv(file_mutsig, sep = '\t')
        mut_sig.rename(columns={'gene':'Hugo_Symbol'},inplace = True)
        
        ### compute average rank for mutation data
        genes = set(wes_df.Hugo_Symbol.values).intersection(set(mut_sig.Hugo_Symbol.values))
        wes_df = wes_df.loc[wes_df.Hugo_Symbol.isin(genes),]
        mut_sig = mut_sig.loc[mut_sig.Hugo_Symbol.isin(genes),]

        wes_df.sort_values(by = 'Freq_MutatedSamples', ascending = False, inplace = True, ignore_index = True)
        wes_df['rank_mut_freq'] = wes_df.index+1

        mut_sig.sort_values(by = 'p', ascending = True, inplace = True, ignore_index = True)
        mut_sig['rank_p'] = mut_sig.index+1

        wes_mut_sig = pd.merge(wes_df,mut_sig,on = 'Hugo_Symbol')
        wes_mut_sig = wes_mut_sig[['Hugo_Symbol','Freq_MutatedSamples','rank_mut_freq','rank_p']]
        wes_mut_sig['avg_rank'] = wes_mut_sig[['rank_mut_freq','rank_p']].apply(gmean,axis = 1)
        wes_mut_sig.sort_values(by = 'avg_rank', ascending = True, inplace = True, ignore_index = True)
        
        ### save wes_mut_sig file
        save_dir1 = os.path.join(save_dir,t_type)
        os.makedirs(save_dir1, exist_ok=True)  # succeeds even if directory exists.

        wes_mutsig_file = os.path.join(save_dir1,'TCGA_'+t_type+'_wes_mut_sig_geo_avg_rank_genes.csv')
        wes_mut_sig.to_csv(wes_mutsig_file, index = False)
        
        ### compute rank for CNV data
        cnv_df = wes_df.copy()
        cnv_df = cnv_df[['Hugo_Symbol','Freq_CNV']]
        cnv_df.sort_values(by = 'Freq_CNV', ascending = False, inplace = True, ignore_index = True)

        ### load SV file
        sv_df = pd.read_csv(sv,sep = '\t')
        sv_df.sort_values(by = '#', ascending = False, inplace = True)
        tmp1 = np.array(sv_df['#'])/np.array(sv_df['Profiled Samples'])
        sv_df.Freq = tmp1
        
        # g1 = cnv_df.loc[cnv_df.Freq_CNV>0.02,'Hugo_Symbol'].values.tolist()
        g1 = cnv_df.Hugo_Symbol.values[:50]
        # g2 = wes_mut_sig.loc[wes_mut_sig.Freq_MutatedSamples>0.02,'Hugo_Symbol'].values.tolist()
        g2 = wes_mut_sig.Hugo_Symbol.values[:50]
        g3 = sv_df.loc[sv_df.Freq>0.02,'Gene'].values.tolist()
        g3 = g3[:50]

        ### make a single dataframe of wes, cnv and sv
        ### WES data patients with altered gene

        tmp1 = cnv_df.loc[cnv_df.Hugo_Symbol.isin(g1),['Hugo_Symbol','Freq_CNV']]
        tmp1 = tmp1.rename(columns={'Hugo_Symbol':'Gene','Freq_CNV':'Freq'})
        tmp2 = wes_mut_sig.loc[wes_mut_sig.Hugo_Symbol.isin(g2),['Hugo_Symbol','Freq_MutatedSamples']]
        tmp2 = tmp2.rename(columns={'Hugo_Symbol':'Gene','Freq_MutatedSamples':'Freq'})
        tmp3 = sv_df.loc[sv_df.Gene.isin(g3),['Gene','Freq']]

        tmp = pd.concat([tmp1,tmp2])
        tmp = pd.concat([tmp,tmp3])

        # check duplicates, mut
        dup_ind = tmp.duplicated(subset = ['Gene'],keep = False)

        # remove duplicates,keep the largest one
        wes_s = tmp.copy()
        wes_s.sort_values(by = 'Freq', ascending = False, inplace = True)
        wes_s = wes_s.drop_duplicates(subset='Gene', keep="first")
        
        # add CNV mutation and SV freq information
        cnv_f = []
        mut_f = []
        sv_f = []
        for gene in wes_s.Gene.values:
            if gene in cnv_df.Hugo_Symbol.values:
                cnv_f.append(cnv_df.loc[cnv_df.Hugo_Symbol==gene,'Freq_CNV'].values[0])
            else:
                cnv_f.append(0)
            if gene in sv_df.Gene.values:
                sv_f.append(sv_df.loc[sv_df.Gene==gene,'Freq'].values[0])
            else:
                sv_f.append(0)
            if gene in wes_mut_sig.Hugo_Symbol.values:
                mut_f.append(wes_mut_sig.loc[wes_mut_sig.Hugo_Symbol==gene,'Freq_MutatedSamples'].values[0])
            else:
                mut_f.append(0)
        wes_s['Freq_CNV'] = cnv_f
        wes_s['Freq_SV'] = sv_f
        wes_s['Freq_Mutation'] = mut_f

        fnm = os.path.join(save_dir1,'TCGA_'+t_type+'_wes_freq.csv')
        wes_s.to_csv(fnm,sep = ',',index = False)

        
def altered_rate_long(dir0,save_dir,t_types):
    
    """
    Save the altered freq of all available genes (WES: Mutation, CNV and SV)
    """

    for t_type in t_types:
        if t_type == 'SKCM':
            dir1 = os.path.join(dir0,'TCGA_'+t_type,'WES postprocessed data from cbioportal','met')
            file = os.path.join(dir1,'TCGA_'+t_type+'_WES_met_gene_summary.csv')  # CNA and mutation data
            sv = os.path.join(dir1,'Structural_Variant_Genes.txt')
        else:
            dir1 = os.path.join(dir0,'TCGA_'+t_type,'WES postprocessed data from cbioportal')
            file = os.path.join(dir1,'TCGA_'+t_type+'_WES_Primary_gene_summary.csv')  # CNA and mutation data
            sv = os.path.join(dir1,'Structural_Variant_Genes.txt')

        ### WES data patients with altered gene
        wes_df = pd.read_csv(file, sep = ',', index_col = 0)
        wes_df.sort_values(by = 'Freq_MutatedSamples', ascending = False, inplace = True, ignore_index = True)
        mut_df = wes_df[['Hugo_Symbol','Freq_MutatedSamples']]

        ### save wes_mut_sig file
        save_dir1 = os.path.join(save_dir,t_type)
        os.makedirs(save_dir1, exist_ok=True)  # succeeds even if directory exists.

        ### compute rank for CNV data
        cnv_df = wes_df.copy()
        cnv_df = cnv_df[['Hugo_Symbol','Freq_CNV']]
        cnv_df.sort_values(by = 'Freq_CNV', ascending = False, inplace = True, ignore_index = True)

        ### load SV file
        sv_df = pd.read_csv(sv,sep = '\t')
        sv_df.sort_values(by = '#', ascending = False, inplace = True)
        tmp1 = np.array(sv_df['#'])/np.array(sv_df['Profiled Samples'])
        sv_df.Freq = tmp1

        g1 = cnv_df['Hugo_Symbol'].values.tolist()
        g2 = mut_df['Hugo_Symbol'].values.tolist()
        g3 = sv_df['Gene'].values.tolist()

        ### make a single dataframe of wes, cnv and sv
        ### WES data patients with altered gene

        tmp1 = cnv_df.loc[cnv_df.Hugo_Symbol.isin(g1),['Hugo_Symbol','Freq_CNV']]
        tmp1 = tmp1.rename(columns={'Hugo_Symbol':'Gene','Freq_CNV':'Freq'})
        tmp2 = mut_df.loc[mut_df.Hugo_Symbol.isin(g2),['Hugo_Symbol','Freq_MutatedSamples']]
        tmp2 = tmp2.rename(columns={'Hugo_Symbol':'Gene','Freq_MutatedSamples':'Freq'})
        tmp3 = sv_df.loc[sv_df.Gene.isin(g3),['Gene','Freq']]

        tmp = pd.concat([tmp1,tmp2])
        tmp = pd.concat([tmp,tmp3])

        # check duplicates, mut
        dup_ind = tmp.duplicated(subset = ['Gene'],keep = False)

        # remove duplicates,keep the largest one
        wes_s = tmp.copy()
        wes_s.sort_values(by = 'Freq', ascending = False, inplace = True)
        wes_s = wes_s.drop_duplicates(subset='Gene', keep="first")

        # add CNV mutation and SV freq information
        cnv_f = []
        mut_f = []
        sv_f = []
        for gene in wes_s.Gene.values:
            if gene in cnv_df.Hugo_Symbol.values:
                cnv_f.append(cnv_df.loc[cnv_df.Hugo_Symbol==gene,'Freq_CNV'].values[0])
            else:
                cnv_f.append(0)
            if gene in sv_df.Gene.values:
                sv_f.append(sv_df.loc[sv_df.Gene==gene,'Freq'].values[0])
            else:
                sv_f.append(0)
            if gene in mut_df.Hugo_Symbol.values:
                mut_f.append(mut_df.loc[mut_df.Hugo_Symbol==gene,'Freq_MutatedSamples'].values[0])
            else:
                mut_f.append(0)
        wes_s['Freq_CNV'] = cnv_f
        wes_s['Freq_SV'] = sv_f
        wes_s['Freq_Mutation'] = mut_f

        fnm = os.path.join(save_dir1,'TCGA_'+t_type+'_wes_freq_long.csv')
        wes_s.to_csv(fnm,sep = ',',index = False)
