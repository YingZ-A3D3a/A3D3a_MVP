import pandas as pd
import PGM_PY.load_interactome as li

# check the edges between a selected gene and a list of genes

def edges(geneoe,g_list):
    df_int = li.load_interactome()
    index1 = (df_int.geneA==geneoe)&(df_int.geneB.isin(g_list))
    index2 = (df_int.geneB==geneoe)&(df_int.geneA.isin(g_list))
    df_select = df_int.loc[index1|index2,]
    return df_select