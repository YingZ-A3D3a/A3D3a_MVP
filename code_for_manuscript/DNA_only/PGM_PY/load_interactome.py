### load the canSAR and KEGG interactome
import os
import pandas as pd

def load_interactome():
    dir0 = 'PGM_PY/other_files'
    fnm = os.path.join(dir0,'CanSAR_KEGG_Hi_union_combined_021524.csv')
    df = pd.read_csv(fnm, low_memory = False)
    return df