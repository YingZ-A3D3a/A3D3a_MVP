import biomart
# Set up connection to server                                               
server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')         
mart = server.datasets['hsapiens_gene_ensembl']
# List the types of data we want                                            
attributes = ['ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol']
response = mart.search({'attributes': attributes})

data = response.raw.data.decode('ascii')
ensembl_id = []
entrez_id = []
hgnc_symbol = []

# Store the data in a dataframe
for line in data.splitlines():                                              
    line = line.split('\t')                                                 
    # The entries are in the same order as in the `attributes` variable                                                 
    ensembl_id.append(line[0])
    entrez_id.append(line[1])
    hgnc_symbol.append(line[2])
ensembl2syb = pd.DataFrame({'ensembl_id':ensembl_id, 'entrez_id':entrez_id, 'hgnc_symbol':hgnc_symbol})
ensembl2syb.head()

### save the gene name mapping
dir2 = '...'
filenm = os.path.join(dir2, 'ensembl_entrez_to_hgnc_human_102122.csv')
ensembl2syb.to_csv(filenm, index = False)
