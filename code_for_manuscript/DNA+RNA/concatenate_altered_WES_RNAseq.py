####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes concatenate WES and RNAseq seed genes
####################################################################################

import os
import pandas as pd
import numpy as np

####################################################################################
def concatenate_wes_rna(save_dir,wes_dir,t_types):

    for t_type in t_types:
        print(t_type)
        fnm = os.path.join(wes_dir,t_type,'TCGA_'+t_type+'_wes_freq.csv')
        wes_s = pd.read_csv(fnm)

        ### rna
        rna_data = os.path.join(save_dir,t_type,'primary_normal_TCGA_GTEx_'+t_type+'_combatseq_protein_coding_log2.csv')
        rna_classifier = os.path.join(save_dir,t_type,'elastic net model','DEGs_with_elastic_info_GTEx+NAT_022123.csv')
        rna_df = pd.read_csv(rna_classifier,sep = ',')

        g4 = rna_df.loc[rna_df.elastic_net_select == 1,]
        deg = g4.gene.values.tolist()
        grp_A = 'TCGA_primary_tumor'
        grp_B = ['TCGA_normal','GTEx_normal']
        
        r_df = pd.read_csv(rna_data,sep = ',',index_col = 0)
        group = pd.Series(r_df.group)
        r_df = r_df.iloc[:,:-1]
        deg_rdf = r_df[deg]
        
        ### count the patients with altered gene
        # if the fc>1 (UP), count the patients in the primary group that have expression above 
        # the 99 percentile of normal group
        # if the fc<1 (DOWN), count the patients in the primary group that have expression below 
        # the 1 percentile of normal group
        cts = []
        for i in range(0,g4.shape[0]):
            genex = g4.gene.values[i]
            A = deg_rdf.loc[group.values==grp_A,genex]
            B = deg_rdf.loc[group.isin(grp_B),genex]
            lth = np.percentile(B,99)
            ltl = np.percentile(B,1)
            if g4.Fold_Change.values[i]>1:
                ct = (A>lth).sum()
            elif g4.Fold_Change.values[i]<1:
                ct = (A<ltl).sum()
            # print(f"fold change: {g4.fc.values[i]}, altered samples: {ct}, percent altered: {ct/len(A)}")
            cts.append(ct)

        ### save dataframe of patients with altered gene
        genes = g4.gene.values
        Nt = (group.values==grp_A).sum()
        rna_alt = pd.DataFrame({'gene':genes,'count':cts,'percentage':[i/Nt for i in cts]})
        rna_alt.sort_values(by = 'percentage', ascending = False, inplace = True)

        rna_s = rna_alt.iloc[:,[0,2]]
        rna_s = rna_s.rename(columns = {'gene':'Gene','percentage':'Freq'})

        # concatenate wes with rna, if a gene presents in both wes and rna, choose the larger freq
        conc = pd.concat([wes_s,rna_s])

        # check duplicates, mut
        dup_ind = conc.duplicated(subset = ['Gene'],keep = False)

        # if there are duplicates,remove duplicates,keep the largest one
        if dup_ind.sum()>0:
            conc.sort_values(by = 'Freq', ascending = False, inplace = True)
            conc = conc.drop_duplicates(subset='Gene', keep="first")

        # add wes or rna info
        in_wes = []
        in_rna = []
        for i in conc.Gene.values:
            if i in wes_s.Gene.values:
                in_wes.append(1)
            else:
                in_wes.append(0)
            if i in rna_s.Gene.values:
                in_rna.append(1)
            else:
                in_rna.append(0)
        conc['wes'] = in_wes
        conc['rna'] = in_rna
        conc.sort_values(by = 'Freq', ascending = False, inplace = True)
        conc = conc.reset_index(drop = True)
        
        file_dir = os.path.join(save_dir,t_type,'elastic net model')
        fnm = os.path.join(file_dir,'wes+rna_classifier_022123.csv')
        conc.to_csv(fnm,index = None)
