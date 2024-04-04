### load the gene name mapping
dir2 = '...'
filenm = os.path.join(dir2, 'ensembl_entrez_to_hgnc_human_102122.csv')
df0 = pd.read_csv(filenm)

dir0 = '...'
fnm = os.path.join(dir0,'HuRI.tsv')
huri = pd.read_csv(fnm, sep = '\t',header = None)

### map ensembl id to symbol
gA = []
gB = []
for i in range(huri.shape[0]):
    eA = huri.iloc[i,0]
    eB = huri.iloc[i,1]
    if (eA in df0.ensembl_id.values) and (eB in df0.ensembl_id.values):
        gA.append(df0.loc[df0.ensembl_id==eA,'hgnc_symbol'].values[0])
        gB.append(df0.loc[df0.ensembl_id==eB,'hgnc_symbol'].values[0])
    
huri_g = pd.DataFrame({'geneA':gA,'geneB':gB})
huri_g.head()

fnm = os.path.join(dir0,'HuRI_symbol.csv')
huri_g.to_csv(fnm, index = False)
