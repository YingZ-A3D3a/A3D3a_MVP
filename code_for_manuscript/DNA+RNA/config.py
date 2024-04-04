t_types = ['BLCA','BRCA',
          'COADREAD','GBM','KIRC',
          'LUAD','LUSC',
          'PRAD','STAD','THCA','UCEC'] # tumor types

to_remove = ['TTN','MUC16']
thre_FN_FDR = 0.05 # threshold for the FDR of the first neighbors
topN_tier1 = 20
topNs = [20,100]

# model configuration setting for number of first neighbors
fn_min = 50
fn_max = 550
fn_step = 100

# parameters for the MC
Wm_min = 0.3
Wm_max = 0.9
MC_alpha_min = 0
MC_alpha_max = 0.1
MC_alpha_step = 0.01

# parameters for the PR, and PPR
alpha_min = 0.05
alpha_max = 0.95
alpha_step = 0.05

# parameters for the diffusion
sigma2_min = 5
sigma2_max = 30
sigma2_step = 5

# depmap t_types
depmap_types = ['BLCA','BRCA',
        ['COAD','READ'],'GB','CCRCC',
        'LUAD','LUSC',
        'PRAD','STAD','THPA','UCEC'] # matched tumor types in depmap