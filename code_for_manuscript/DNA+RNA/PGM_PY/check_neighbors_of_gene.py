import pandas as pd
import numpy as np
import PGM_PY.check_edge_type as CET

def check_neighbors(geneoe,glist):
    df_select = CET.edges(geneoe,glist)

    s1 = df_select[['geneA','geneB']]
    s1 = s1.drop_duplicates()
    s1 = s1.dropna()

    s1['geneA'] = s1['geneA'].astype(str)
    s1['geneB'] = s1['geneB'].astype(str)

    s3 = pd.DataFrame(np.sort(s1.values, axis=1), columns=s1.columns).drop_duplicates()
    # print(s3.shape)

    geneA = list(set(s3['geneA'].values))
    temp = geneA.copy()
    geneB = list(set(s3['geneB'].values))
    geneA.extend(geneB)
    # print(f"geneA:{len(temp)},geneB:{len(geneB)},total genes:{len(set(geneA))}")
    nb = [i for i in geneA if i!=geneoe]
    # print(f"neighbors of {geneoe}: {len(nb)}")
    return nb