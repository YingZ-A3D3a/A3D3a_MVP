####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes perform DEGs for TCGA and GTEx
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import ranksums
from statsmodels.stats.multitest import fdrcorrection

####################################################################################
def compute_DEG(save_dir,t_types):

    for t_type in t_types:
        print(t_type)
        file_dir = os.path.join(save_dir,t_type)
        fnm = os.path.join(file_dir,'primary_normal_TCGA_GTEx_'+t_type+'_combatseq_protein_coding_log2.csv')
        df_a = pd.read_csv(fnm,index_col = 0)
        group = df_a.group
        
        # check how many samples have many 0 in gene expression, e.g. >70% of genes are of 0 expression
        tmp = df_a.iloc[:,:-1]
        ct = 0
        p_id = [] # sample id that has problem
        for i in range(tmp.shape[0]):
            if sum(tmp.iloc[i,:]!=0)/tmp.shape[1]<0.3:
                p_id.append(i)
                ct+=1
        print(ct)

        ## check the samples with problem
        for i in range(len(p_id)):
            print(tmp.iloc[p_id[i],].name)

        # remove the problematic sample
        for i in range(len(p_id)):
            df_a = df_a.drop([tmp.iloc[p_id[i],].name])
            
        ### select groups for comparions, combine GTEx_normal and TCGA_normal into the same normal group
        grp_cp = ['TCGA_normal','GTEx_normal']
        tmp = []
        for i in df_a.group.values:
            if i in grp_cp:
                tmp.append('normal')
            else:
                tmp.append('primary_tumor')
        df = df_a.copy()
        df.group = tmp
        
        p_values = []
        fold_change = []
        abs_lfc = [] # abolute log2 fold change
        mean_A = []
        mean_B = []

        grpA = df.loc[df.group=='primary_tumor',] # primary
        grpA = grpA.iloc[:,:-1]
        grpB = df.loc[df.group=='normal',] # normal
        grpB = grpB.iloc[:,:-1]

        for i in range(grpA.shape[1]):
            A = 2**(grpA.iloc[:,i].values)
            B = 2**(grpB.iloc[:,i].values)
            p_values.append(ranksums(A,B).pvalue)
            fold_change.append(np.mean(A)/np.mean(B))
            abs_lfc.append(abs(np.log2(np.mean(A)/np.mean(B))))
            mean_A.append(np.mean(A))
            mean_B.append(np.mean(B))

        p_df = pd.DataFrame({'p_value': p_values, 'mean_tumor':mean_A,'mean_normal':mean_B,'Fold_Change': fold_change,'abs_log2_fc': abs_lfc}, index = grpA.columns)
        p_df.sort_values(by = ['p_value','abs_log2_fc'],inplace = True, ascending = [True,False])
        p_df['gene'] = p_df.index

        ### remove inf cases
        p_df = p_df.loc[p_df.abs_log2_fc!=float('inf'),]
        
        fdrs = fdrcorrection(p_df.p_value.values)[1]
        p_df['FDR_BH'] = fdrs
        p_df = p_df[['gene','mean_tumor','mean_normal','Fold_Change','abs_log2_fc','p_value','FDR_BH']]
        
        ### save DEGs
        save_fnm= os.path.join(file_dir,'DEG_all_combatseq_TCGA_primary_vs_GTEx+NAT_022123.csv')
        p_df.to_csv(save_fnm,index = False)

        
    """
    Count the patients with altered gene
    If the fc>1 (UP), count the patients in the primary group that have expression above the 99 percentile of normal group
    If the fc<1 (DOWN), count the patients in the primary group that have expression below the 1 percentile of normal group
    """
    for t_type in t_types:
        ### rna, expression data
        file_dir = os.path.join(save_dir,t_type)
        rna_data = os.path.join(file_dir,'primary_normal_TCGA_GTEx_'+t_type+'_combatseq_protein_coding_log2.csv')
        ### rna, fold change between groups
        rna_deg = os.path.join(file_dir,'DEG_all_combatseq_TCGA_primary_vs_GTEx+NAT_022123.csv')
        grp_A = 'TCGA_primary_tumor'
        grp_B = ['TCGA_normal','GTEx_normal']

        r_df = pd.read_csv(rna_data,sep = ',',index_col = 0)
        group = pd.Series(r_df.group)
        r_df = r_df.iloc[:,:-1]

        deg_df = pd.read_csv(rna_deg,sep = ',',index_col = 0)

        genes = r_df.columns

        cts = []
        for i in range(len(genes)):
            genex = genes[i]
            A = r_df.loc[group.values==grp_A,genex]
            B = r_df.loc[group.isin(grp_B),genex]
            lth = np.percentile(B,99)
            ltl = np.percentile(B,1)
            if deg_df.loc[genex,'Fold_Change']>1:
                ct = (A>lth).sum()
            elif deg_df.loc[genex,'Fold_Change']<1:
                ct = (A<ltl).sum()
            # print(f"fold change: {g4.fc.values[i]}, altered samples: {ct}, percent altered: {ct/len(A)}")
            cts.append(ct)

        ### save dataframe of patients with altered gene
        Nt = (group.values==grp_A).sum()
        rna_alt = pd.DataFrame({'gene':genes,'count':cts,'percentage':[i/Nt for i in cts]})
        rna_alt.sort_values(by = 'percentage', ascending = False, inplace = True)

        rna_s = rna_alt.iloc[:,[0,2]]
        rna_s = rna_s.rename(columns = {'gene':'Gene','percentage':'Freq'})
        fnm = os.path.join(file_dir,'TCGA_'+t_type+'_rna_altered_freq_all_genes.csv')
        rna_s.to_csv(fnm,index = None)
