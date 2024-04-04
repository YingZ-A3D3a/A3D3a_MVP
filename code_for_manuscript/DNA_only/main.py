import os
import subprocess
import logging
logging.basicConfig(level=logging.WARNING)

from config import t_types,to_remove,thre_FN_FDR
from config import topN_tier1,topNs
from config import fn_min,fn_max,fn_step,alpha_min,alpha_max,alpha_step
from config import depmap_types

root = '...'  # root directory
dir0 = os.path.join(root,'input')
# save directory
save_dir = os.path.join(root,'results','DNA_only')
os.makedirs(save_dir, exist_ok=True)  # succeeds even if directory exists.

""" Save the altered freq of seed genes (DNA: Mutation, CNV and SV) """
from patients_altered_WES import seed_genes as sg
sg(dir0,save_dir,t_types)

""" Save the altered freq of all available genes (DNA: Mutation, CNV and SV) """
from patients_altered_WES import altered_rate_long as arl
arl(dir0,save_dir,t_types)

""" Find first neighbors for each cancer type by permutation test """
from find_first_neighbors import find_fn
find_fn(save_dir,to_remove,t_types)

""" Tuning the parameters by grid search """
from tuning_PGM_parameters import tuning_pgm
tuning_pgm(save_dir,t_types,to_remove,thre_FN_FDR)

""" Optimize the MVP models by tier 1 genes """
from evaluation_of_models_by_tier1_genes import eval_mc_tier1 as emt
emt(dir0,save_dir,t_types,topN_tier1)

""" Run the PR, PPR models """
from prep_model_comparisons import prep_model_comparisons as pmc
pmc(save_dir,t_types,fn_min, fn_max, fn_step, alpha_min, alpha_max, alpha_step)

""" Run the diffusion models """
fnm = 'diffusion_models.R'
# Call the R script with input parameters using subprocess.Popen()
process = subprocess.Popen(['Rscript', fnm, save_dir]+t_types)
# Wait for the script to finish execution
process.wait()

""" Model tuning for PR PPR and difffusion, and model comparisons """
from model_comparisons import model_comparisons_pipeline as mcp
mcp(dir0,save_dir,t_types,topNs)

""" Compare study bias between models """
from study_bias_check import study_bias as sb
sb(save_dir,t_types,topNs)

""" These codes perform depmap validation """
from depmap_validation import depmap_validation as dv
dv(dir0,save_dir,t_types,depmap_types,topN = 100) # topN genes can be input as a parameter

""" These codes perform depmap saturation check """
from depmap_saturation_check import depmap_validation_detailed as dvd
dvd(dir0,save_dir,t_types,depmap_types)

""" These codes perform GDSC validation """
from GDSC_validation import GDSC_validation
GDSC_validation(dir0,save_dir,t_types,topN = 20)

""" These codes perform community detection """
from community_detection import community_detection_pipeline as cdp
cdp(dir0,save_dir,t_types,plot_topN = 60)

""" These codes set initial score of driver/onco genes to zero """
from set_initial_score_to_zero import set_initial_to_zero as sitz
sitz(dir0,save_dir)
