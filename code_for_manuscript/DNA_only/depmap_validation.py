####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes perform depmap validation
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import seaborn as sns

####################################################################################
### model name
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉") # function to convert to subscript
example_string = "A3D3a"
Model_nm = example_string.translate(SUB)+"'s MVP"
    
def depmap_validation(dir0,save_dir,t_types,mt_types,topN = 20):

    fig_dir = os.path.join(save_dir,'figures')
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
    plot_essential_and_common_essentials(t_types,ge_list_all_rm_ce,ge_list_all,num_cell_lines,ce_g,fig_dir,mat_dir)
    
    # including common essentials
    compare_top_predicted_and_top_altered(t_types,df_rnk,ge_list_all,num_cell_lines,overlap_ce,dir0,save_dir,fig_dir,mat_dir,ftype = 'including common essentials',topN = topN)
    
    # excluding common essentials
    compare_top_predicted_and_top_altered(t_types,df_rnk,ge_list_all_rm_ce,num_cell_lines,overlap_ce,dir0,save_dir,fig_dir,mat_dir,ftype = 'excluding common essentials',topN = topN)
    
def plot_essential_and_common_essentials(t_types,ge_list_all_rm_ce,ge_list_all,num_cell_lines,ce_g,fig_dir,mat_dir):
    ### distribution of essential genes, common essential genes removed
    ng_rm_ce = {} # number of essential genes for each cell line for every type of cancer, remove common essentials
    for t_type in t_types:
        ng_rm_ce[t_type] = [len(i) for i in ge_list_all_rm_ce[t_type]]
    ng = {} # number of essential genes for each cell line for every type of cancer
    for t_type in t_types:
        ng[t_type] = [len(i) for i in ge_list_all[t_type]]        
    labels, data1 = [*zip(*ng_rm_ce.items())]  # 'transpose' items to parallel key, value lists
    plt.rcParams["figure.figsize"] = (2,2)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    bpl = plt.boxplot(data1, positions=np.array(range(len(data1)))*2.0, sym='', widths=0.9)
    labels, data2 = [*zip(*ng.items())]  # 'transpose' items to parallel key, value lists
    bpr = plt.boxplot(data2, positions=np.array(range(len(data2)))*2.0, sym='', widths=0.9)
    colorl = '#E64B35FF'
    colorr = '#3C5488FF'
    set_box_color(bpl, colorl) # colors are from http://colorbrewer2.org/
    set_box_color(bpr, colorr)
    labels1 = [i+' ('+str(num_cell_lines[i])+')' for i in labels]
    plt.xticks(np.arange(0, 2*len(labels),2), labels1, rotation = 90, fontsize = 5)
    plt.yticks([0,500,1000,1500], fontsize = 5)
    plt.xlim([-1,2*len(labels)-1])
    plt.ylim([0,1500])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.plot([], c=colorl, label='without common essentials')
    plt.plot([], c=colorr, label='with common essentials')
    plt.legend(fontsize = 5)
    title = f'depmap essential genes'
    # plt.title(title, fontsize = 14)
    fnm = os.path.join(fig_dir,title+'.svg')
    plt.ylabel('Count',fontsize = 6)
    plt.tight_layout()
    plt.savefig(fnm,dpi = 300)

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
    x = np.arange(0,2*len(ng_pool),2)
    width = 0.40    
    # plot data in grouped manner of bar type
    plt.rcParams["figure.figsize"] = (2.5,2)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.bar(x, ng_rm_ce_pool, width,color = '#E64B35FF')
    plt.bar(x+0.5, ng_pool, width,color = '#3C5488FF')
    plt.xticks(x+0.25,labels1,fontsize = 5, rotation = 90)
    plt.yticks(fontsize = 5)
    title = f'depmap essential genes pooled all cell lines for each cancer type'
    # plt.title(title,fontsize = 5)
    plt.ylabel('Count',fontsize = 5)
    plt.ylim([0,2500])
    plt.legend(["without common essentials", "with common essentials"],fontsize = 5)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,title+'_common_essential_removed.svg')
    plt.savefig(fnm,dpi = 300)
    print(f"mean and std depmap essentials without common essentials: {np.mean(ng_rm_ce_pool)},{np.std(ng_rm_ce_pool)}")
    print(f"mean and std depmap essentials with common essentials: {np.mean(ng_pool)},{np.std(ng_pool)}")
    # save essential genes
    save_essential_genes(g_pool,ce_g,mat_dir)

def compare_top_predicted_and_top_altered(t_types,df_rnk,ge_list_all,num_cell_lines,overlap_ce,dir0,save_dir,fig_dir,mat_dir,ftype,topN):
    ### protein coding genes
    dir2 = 'databases/TCGA_GTEx'
    fnm = os.path.join(dir0,dir2,'TcgaTargetGTEX_gene_biotype_by_biomaRt.csv')
    fx = pd.read_csv(fnm)
    pcg = fx.loc[fx.gene_biotype=='protein_coding','hgnc_symbol'].values.tolist()
    print(f'number of protein coding genes {len(pcg)}')
    
    ### best Markov model, check depmap genes in the disease network
    df_rnk = df_rnk.sort_values(by = 'mean_rank')
    select_md = df_rnk.index[0]
    overlap_ct = {} # overlap percent between top-ranked and essential for every cancer type
    overlap_all = {} # overlapped genes between top-ranked and essential
    overlap_ct_fq = {} # overlap percent between top-altered and essential for every cancer type
    overlap_all_fq = {} # overlapped genes between top-altered and essential

    for t_type in t_types:
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        df = df.iloc[:topN,]
        g_d = df.genes.values.tolist()    # all disease genes
        
        fq_dir = os.path.join(save_dir, t_type, 'TCGA_'+t_type+'_wes_freq_long.csv')
        df_q = pd.read_csv(fq_dir)
        df_q = df_q.loc[df_q.Gene.isin(pcg),]
        df_q = df_q.sort_values(by = ['Freq'],ascending = False)
        df_q = df_q.reset_index(drop = True)
        top_f = df_q.iloc[:topN,0].values.tolist()
        
        ge_list = ge_list_all[t_type]
        overlap_num = []  # overlap of essential genes and top predicted genes
        overlap = []
        overlap_num_fq = []  # overlap of essential genes and top altered genes
        overlap_fq = []
        
        for l in ge_list:
            tmp = len(set(g_d).intersection(l))
            overlap_num.append(tmp)
            overlap.extend(set(g_d).intersection(l))
            
            tmp = len(set(top_f).intersection(l))
            overlap_num_fq.append(tmp)
            overlap_fq.extend(set(top_f).intersection(l))
            
        overlap_ct[t_type] = overlap_num
        overlap_all[t_type] = overlap
        overlap_ct_fq[t_type] = overlap_num_fq
        overlap_all_fq[t_type] = overlap_fq
    filenm = f'Depmap essential genes among top {topN} disease network genes {ftype} boxplot'
    multiple_boxplots(overlap_ct,overlap_ct_fq,num_cell_lines,filenm,fig_dir,ftype,topN = topN)
    
    ### save the number of depmap essentials to a csv
    tumor_type = []
    number_of_depmap_essentials = []
    for t_type in overlap_ct:
        tumor_type.extend([t_type for i in overlap_ct[t_type]])
        number_of_depmap_essentials.extend(overlap_ct[t_type])
    num_dep_ess = pd.DataFrame({'tumor_type':tumor_type,'number_of_depmap_essentials':number_of_depmap_essentials})
    fnm = os.path.join(mat_dir,f'number_of_depmap_essentials_top_{topN}_including_common_every_cell_line.csv')
    num_dep_ess.to_csv(fnm)
    
    ### save dataframe, topN ranked genes, whether depmap essential
    c_types = []
    top_genes = []
    if_essential = []
    if_common_essential = []
    top_index = []
    for t_type in t_types:
        c_types.extend([t_type for i in range(topN)])
        top_index.extend([i+1 for i in range(topN)])
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        df = df.iloc[:topN,]
        g_d = df.genes.values.tolist()    # all disease genes
        top_genes.extend(g_d)
        for item in g_d:
            if item in overlap_all[t_type]:
                if_essential.append(1)
            else:
                if_essential.append(0)
        for item in g_d:
            if item in overlap_ce:
                if_common_essential.append(1)
            else:
                if_common_essential.append(0)            
    top_ess_mat = pd.DataFrame({'cancer_type':c_types,'top_index':top_index,
                            'top_genes':top_genes,'depmap_essential':if_essential,
                            'common_essential':if_common_essential})
    print(top_ess_mat.shape)
    fnm = os.path.join(mat_dir,f'top_{topN}_genes_MVP_with_depmap_essential_info_{ftype}.csv')
    top_ess_mat.to_csv(fnm,index = False)
    
    ### check how many essential genes are first neighbors
    labels1 = [i+' ('+str(num_cell_lines[i])+')' for i in t_types]
    N_ess,N_ess_fn,N_ess_fq = check_first_neighors_are_essential(t_types,overlap_all,overlap_all_fq,labels1,select_md,topN,save_dir,fig_dir)
    statistical_test(N_ess,N_ess_fq,topN,fig_dir)
    if ftype=='including common essentials':
        save_table_essentials(t_types,N_ess,N_ess_fn,select_md,topN,overlap_all,overlap_all_fq,overlap_ce,pcg,save_dir,mat_dir)
    
def save_essential_genes(g_pool,ce_g,mat_dir):
    ### save the depmap essentials for every cancer type
    df_es = pd.DataFrame()
    for t_type in g_pool:
        tmp = g_pool[t_type]
        t1 = [t_type for i in range(len(tmp))]
        c1 = [int(i in ce_g) for i in tmp]
        tmp_df = pd.DataFrame({'tumor type':t1,'depmap essentials':tmp,'common essentials':c1})
        df_es = pd.concat([df_es, tmp_df], axis = 0)
    print(df_es.shape)
    df_es.head()
    fnm = os.path.join(mat_dir,'depmap essential genes for every tumor type.csv')
    df_es.to_csv(fnm, index = False)

def set_box_color(bp, color, linewidth0=0.75):
    ### set color of boxplot
    plt.setp(bp['boxes'], color=color,linewidth = linewidth0)
    plt.setp(bp['whiskers'], color=color,linewidth = linewidth0)
    plt.setp(bp['caps'], color=color,linewidth = linewidth0)
    plt.setp(bp['medians'], color='k',linewidth = linewidth0)

### make multiple boxplot comarisons between two groups
def multiple_boxplots(dict1, dict2, num_cell_lines, filenm, fig_dir, ftype, face_colors = ['#91D1C2FF', '#B09C85FF'], 
                      edge_colors = ['#00A087FF','#7E6148FF'],figsize = [6,3], topN = 20):
    """input data are two dictionaries, where keys can be cancer types"""
    labels, data = [*zip(*dict1.items())]  # 'transpose' items to parallel key, value lists
    plt.rcParams["figure.figsize"] = (figsize[0],figsize[1])
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()

    # box edge color
    color1,color2 = edge_colors[0],edge_colors[1]
    colors = face_colors
    # fill with colors
    boxpropsl = dict(linestyle='-', linewidth=1, facecolor=colors[0], edgecolor='k',alpha = 0.5)
    boxpropsr = dict(linestyle='-', linewidth=1, facecolor=colors[1], edgecolor='k',alpha = 0.5)

    bpl = plt.boxplot(data, positions=np.array(range(len(data)))*5.0-1, sym = '', widths=1.7, patch_artist=True,boxprops = boxpropsl)

    for i in range(len(data)):
        plt.scatter([i*5.0-1+np.random.normal(scale = 0.15) for k in data[i]],data[i],s = .02,c = 'k')

    labels, data_fq = [*zip(*dict2.items())]  # 'transpose' items to parallel key, value lists
    bpr = plt.boxplot(data_fq, positions=np.array(range(len(data_fq)))*5.0+1, sym = '', widths=1.7, patch_artist=True,boxprops = boxpropsr)

    for i in range(len(data_fq)):
        plt.scatter([i*5.0+1+np.random.normal(scale = 0.15) for k in data_fq[i]],data_fq[i],s = .02,c = 'k')
        tmp1 = max(data[i])+1
        tmp2 = max(data_fq[i])+1
        tmp = max(tmp1,tmp2)
        plt.vlines(x=i*5.0-1,ymin = tmp,ymax=tmp+topN/100)
        plt.vlines(x=i*5.0+1,ymin = tmp,ymax=tmp+topN/100)
        plt.hlines(y=tmp+topN/100,xmin = i*5.0-1,xmax = i*5.0+1)
        x = data[i]
        y = data_fq[i]
        z = ranksums(x,y)
        xlocation = i*5.0-1.2
        ylocation = tmp+(topN+5)/100
        if (z[1]<0.05)&(z[1]>=0.01):
            plt.text(x=xlocation,y=ylocation,s='*')
        elif (z[1]<0.01)&(z[1]>=0.001):
            plt.text(x=xlocation,y=ylocation,s='**')
        elif z[1]<0.001:
            plt.text(x=xlocation,y=ylocation,s='***')
        else:
            plt.text(x=xlocation,y=ylocation,s='n.s.')

    set_box_color(bpl, color1) # colors are from http://colorbrewer2.org/
    set_box_color(bpr, color2)

    labels1 = [i+' ('+str(num_cell_lines[i])+')' for i in labels]
    plt.xticks(np.arange(0, 5*len(labels),5), labels1, fontsize = 7, rotation = 90)
    plt.yticks(fontsize = 7)
    plt.xlim([-3,5*len(labels)-2])
    # plt.ylim([0,tmp])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)

    plt.plot([], c=color1, label=Model_nm)
    plt.plot([], c=color2, label='Top altered')
    plt.legend(fontsize = 7)
    plt.title(f'Top {topN}, {ftype}', fontsize = 5)
    plt.ylabel(f'Count of genetic dependencies')
    plt.tight_layout()
    fnm = os.path.join(fig_dir,filenm+f'_{ftype}.png')
    plt.savefig(fnm,dpi = 256)
    fnm = os.path.join(fig_dir,filenm+f'_{ftype}.svg')
    plt.savefig(fnm,dpi = 256)
    plt.show()
    
def check_first_neighors_are_essential(t_types,overlap_all,overlap_all_fq,labels1,select_md,topN,save_dir,fig_dir):
   ### check how many essential genes are first neighbors
    N_ess = []
    N_ess_fn = []
    N_ess_fq = []
    for t_type in t_types:
        genes = list(set(overlap_all[t_type]))
        N_ess.append(len(genes))        
        genes_fq = list(set(overlap_all_fq[t_type]))
        N_ess_fq.append(len(genes_fq))        
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        fn = df.loc[df.source=='first neighbor','genes'].values.tolist()
        N_ess_fn.append(len(set(genes).intersection(fn)))
    x = np.arange(0,2*len(N_ess),2)
    width = 0.4    
    # plot data in grouped manner of bar type
    plt.rcParams["figure.figsize"] = (3,2.2)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    plt.bar(x, N_ess, width)
    plt.bar(x, N_ess_fn, width)
    plt.bar(x+0.5, N_ess_fq, width,color = 'gray')
    plt.xticks(x+0.25,labels1,fontsize = 6,rotation = 90)
    plt.yticks(fontsize = 6)
    title = f'Count of genetic dependencies'
    # plt.title('Including common essentials')
    plt.ylabel(title,fontsize = 7)
    plt.ylim([0,topN/2+5])
    plt.legend(["Seed", "First neighbor",f'Top altered'],fontsize = 5)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,title+f' top {topN} including common essentials.svg')
    plt.savefig(fnm,dpi = 300, bbox_inches='tight')
    return N_ess,N_ess_fn,N_ess_fq

def save_table_essentials(t_types,N_ess,N_ess_fn,select_md,topN,overlap_all,overlap_all_fq,overlap_ce,pcg,save_dir,mat_dir):
    ### create a table for comparison depmap essentials between MC and top altered
    tumor_types = []
    num_depmap_essential = []
    types = []
    percent_depmap_essential = []
    depmap_essential_genes = []
    common_essential_genes = []
    
    for t_type in t_types:
        ### MC model ###
        md_dir = os.path.join(save_dir, t_type, 'markov chain models')
        s_md = 'markov_output_wes_'+select_md+'_022123.csv'
        df = pd.read_csv(os.path.join(md_dir,s_md))
        df = df.iloc[:topN,]
        g_d = df.genes.values.tolist()    # all disease genes        
        tmp_dep_ess = []
        tmp_com_ess = []        
        for item in g_d:
            if item in overlap_all[t_type]:
                tmp_dep_ess.append(item)
            if item in overlap_ce:
                tmp_com_ess.append(item)
        num_depmap_essential.append(len(tmp_dep_ess))
        tumor_types.append(t_type)
        types.append(f'MVP')
        depmap_essential_genes.append(','.join(sorted(tmp_dep_ess)))
        common_essential_genes.append(','.join(sorted(tmp_com_ess)))
        percent_depmap_essential.append(100*len(tmp_dep_ess)/topN)        
        ### top altered ###
        fq_dir = os.path.join(save_dir, t_type, 'TCGA_'+t_type+'_wes_freq_long.csv')
        df_q = pd.read_csv(fq_dir)
        df_q = df_q.loc[df_q.Gene.isin(pcg),]
        df_q = df_q.sort_values(by = 'Freq',ascending = False)
        df_q = df_q.reset_index(drop = True)
        top_f = df_q.iloc[:topN,0].values.tolist()        
        tmp_dep_ess = []
        tmp_com_ess = []        
        for item in top_f:
            if item in overlap_all_fq[t_type]:
                tmp_dep_ess.append(item)
            if item in overlap_ce:
                tmp_com_ess.append(item)
        num_depmap_essential.append(len(tmp_dep_ess))
        tumor_types.append(t_type)
        types.append(f'Top_altered')
        depmap_essential_genes.append(','.join(sorted(tmp_dep_ess)))
        common_essential_genes.append(','.join(sorted(tmp_com_ess)))
        percent_depmap_essential.append(100*len(tmp_dep_ess)/topN)    
    summary_table = pd.DataFrame({'Tumor type':tumor_types,'Method':types,
                                f'Count of genetic dependencies':num_depmap_essential,
                                f'Percent of genetic dependencies (%)':percent_depmap_essential,
                                f'Genetic dependencies':depmap_essential_genes,
                                f'Common essentials':common_essential_genes})
    # ### save dataframe
    # fnm = os.path.join(mat_dir,f'compare_genetic_dependencies_top_{topN}_ranked_and_top_altered.csv')
    # summary_table.to_csv(fnm, index = False)
    # ### save data frame of number of pooled genetic dependencies for every cancer type
    # num_pool_dep_ess = pd.DataFrame({'tumor_type':t_types,
    #                                 'number_of_genetic_dependencies':N_ess,
    #                                 'number_of_genetic_dependencies_that_are_first_neighbors':N_ess_fn})
    # fnm = os.path.join(mat_dir,f'number_of_pooled_genetic_dependencies_among_top_{topN}.csv')
    # num_pool_dep_ess.to_csv(fnm,index = False)

def statistical_test(N_ess,N_ess_fq,topN,fig_dir,text_location = 0.3):
    ### statistical test
    st = ranksums(N_ess,N_ess_fq)
    tmp = pd.DataFrame({Model_nm:N_ess,'top altered':N_ess_fq})
    plt.rcParams["figure.figsize"] = (1.5,1.7)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['left'].set_linewidth(0.75)
    colors = ['#91D1C2FF', '#B09C85FF']
    my_pal = {Model_nm: colors[0], "top altered":colors[1]}
    sns.boxplot(data = tmp, palette=my_pal, linewidth = 0.5, flierprops={"marker": "x", "markersize":0})
    plt.scatter([0+np.random.normal(scale = 0.01) for i in N_ess],N_ess,s = .5,c = 'k')
    plt.scatter([1+np.random.normal(scale = 0.01) for i in N_ess_fq],N_ess_fq,s = .5,c = 'k')
    
    tmp1 = 1.1*max(N_ess)
    tmp2 = 1.1*max(N_ess_fq)
    tmp0 = max(tmp1,tmp2)
    plt.vlines(x=0,ymin = tmp0,ymax=tmp0)
    plt.vlines(x=1,ymin = tmp0,ymax=tmp0)
    plt.hlines(y=tmp0, xmin = 0, xmax = 1)

    xlocation = text_location
    ylocation = tmp0
    z = st
    if (z[1]<0.05)&(z[1]>=0.01):
        plt.text(x=xlocation,y=ylocation,s='*')
    elif (z[1]<0.01)&(z[1]>=0.001):
        plt.text(x=xlocation,y=ylocation,s='**')
    elif z[1]<0.001:
        plt.text(x=xlocation,y=ylocation,s='***')
    else:
        plt.text(x=xlocation,y=ylocation,s='n.s.')
        
    plt.ylabel(f'Count of genetic dependencies',fontsize = 7)
    # plt.title(f'Including common essentials, p={np.round(st[1],4)}')
    plt.xticks(rotation = 90, fontsize = 6)
    plt.yticks(fontsize = 6)
    plt.tight_layout()
    fnm = os.path.join(fig_dir,f'compare_MC_and_top_{topN}_altered_essential_genes_including_common_essentials.svg')
    plt.savefig(fnm,dpi = 300)
