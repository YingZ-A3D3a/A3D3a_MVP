####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes perform data transformation of TCGA and GTEx, extract protein coding genes, and make UMAP
####################################################################################

library(dplyr)
library(stringr)
library(resample)
library(ggplot2)
library(umap)

####################################################################################

args = commandArgs(trailingOnly = TRUE)
dir0 = args[1]
sav_dir = args[2]
t_types = args[3]

dir1 = paste0(dir0,'/TCGA_GTEx')
dir_combat = paste0(dir0,'/TCGA_GTEx/Combatseq/All_Cancers')

# load gene_biotype
fnm = paste0(dir1,'/TcgaTargetGTEX_gene_biotype_by_biomaRt.csv')
fx = read.delim2(fnm,sep = ',')
fx_p = fx[fx$gene_biotype=='protein_coding',]
print(dim(fx_p))

### optional, check duplicated mapping
temp = length(unique(fx_p[,'ensembl_gene_id']))
temp = fx_p[,'ensembl_gene_id']
dup = temp[duplicated(temp)]
print(dup)
fx_p[fx_p$'ensembl_gene_id' %in% dup,]
fx_p = fx_p %>% distinct(ensembl_gene_id, .keep_all = TRUE) # keep the first row of values for duplicated ensembl id
print(dim(fx_p))

for (t_type in t_types){
    print(t_type)
    filenm = paste(dir_combat, paste0(t_type,'.profile.ComBatSeq.nor.txt'),sep = '/')
    dt = read.delim(filenm, sep = '\t')
    print(dim(dt))
    ## indx
    indx = c()
    x1 = row.names(dt)
    for (i in seq(1,length(x1))){
        if (x1[i] %in% fx_p$ensembl_gene_id){
            indx = c(indx,i)
        }
    }
    dt1 = dt[indx,]
    x1_p = x1[indx]
    dt1['sample_id'] = x1_p
    # map gene name
    temp = dt1[,'sample_id']
    row.names(fx_p) = fx_p$ensembl_gene_id
    g_nm = fx_p[temp,'hgnc_symbol']
    dt1['hgnc_symbol'] = g_nm
    colnm = colnames(dt1)
    new_colnm = str_replace_all(colnm,'\\.','-')
    colnames(dt1) = new_colnm

    ### data transformation
    df1 = t(dt1)
    n = dim(df1)[1]
    df1 = df1[1:(n-2),]
    print(dim(df1))
    colnames(df1) = dt1$hgnc_symbol
    head(df1, n = 3)
    df1 = as.data.frame(df1)
    df2 = data.frame(sapply(df1,function(x) as.numeric(x)))
    row.names(df2) = row.names(df1)

    ### log2(count+1) transformation
    magic_fun <- function(x){
        x = as.numeric(x)
        y = log2(x+1)
        return(y)
    }
    samples = colnames(df2)
    genes = row.names(df2)

    tmp = as.data.frame(df2)
    df3 = data.frame(lapply(tmp,magic_fun))
    colnames(df3) = samples
    row.names(df3) = genes
                            
    ### top 3000 genes with largest variances across samples
    save_dir = paste(sav_dir,t_type,sep = '/')
    dir.create(file.path(save_dir), showWarnings = FALSE)

    vars = colVars(as.matrix(df3[sapply(df3, is.numeric)]))
    tmp = sort(vars,decreasing = TRUE, index.return = TRUE)
    top3000_indx = tmp$ix[1:3000]
    top3000_genes_mat = df3[,top3000_indx]

    samples = row.names(df3)
    groups = c()
    for (i in seq(1,length(samples))){
        tmp = str_split(samples[i],'-')
        if (tmp[[1]][1]=='GTEX'){
            groups = c(groups,'GTEx_normal')
        }
        else if (str_sub(tmp[[1]][4],1,2)=='01'){
            groups = c(groups,'TCGA_primary_tumor')
        }
        else if (str_sub(tmp[[1]][4],1,2)=='11'){
            groups = c(groups,'TCGA_normal')
        }
        else if (str_sub(tmp[[1]][4],1,2)=='06'){
            groups = c(groups,'TCGA_metastatic_tumor')
        } 
        else if (str_sub(tmp[[1]][4],1,2)=='02'){
            groups = c(groups,'TCGA_recurrent_tumor')
        }    
        else {
            groups = c(groups,'other')
        }
    }

    groups_df = data.frame(groups)
    row.names(groups_df) = samples
    groups_df$sample_id = samples

    set.seed(1)
    df.umap = umap(top3000_genes_mat, n_neighbors = 15)  # default n_neighbors = 15
    umap_df = as.data.frame(df.umap$layout)
    colnames(umap_df) = c('UMAP_1','UMAP_2')
    umap_df$group = groups
    fnm = paste0(save_dir,'/umap.png')
    # png(fnm,width=6,height=6,units="in",res=1200)
    tmp_plot = ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = group))+geom_point(size = 2)+
        theme_classic()+theme(legend.position="bottom")+
        theme(axis.text = element_text(size = 10))
    ggsave(tmp_plot, file = fnm, width=10,height=10,units="in",dpi=1024)
                            
    ### select normal and primary tumor in TCGA
    is_select = groups_df[,'groups'] %in% c('TCGA_primary_tumor','TCGA_normal','GTEx_normal')
    tmp = groups_df[is_select,'sample_id']
    df_s = df3[tmp,]
    gps = groups_df$groups[is_select]
    df_s$group = gps
                            
    fnm = paste0(save_dir,
             '/primary_normal_TCGA_GTEx_',t_type,'_combatseq_protein_coding_log2.csv')
    write.table(df_s,fnm,sep = ',',col.names = NA)
}
