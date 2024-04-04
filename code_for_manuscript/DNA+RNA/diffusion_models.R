####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes run diffusion models
####################################################################################

library(igraph)
library(diffuStats)
library(ggplot2)
library(ggsci)
library(igraphdata)

####################################################################################
args = commandArgs(trailingOnly = TRUE)
save_dir = args[1]
t_types = args[2]

diff_models <- function(save_dir,t_types,sigma2_min = 5, sigma2_max = 30, 
                        sigma2_step = 5, fn_min = 50, fn_max = 550, fn_step = 100, 
                        n_perm = 1000){

    sigma_2s = seq(sigma2_min,sigma2_max,sigma2_step)
    fn_nums = seq(fn_min,fn_max,fn_step)

    for (t_type in t_types){
        om_dir = paste(save_dir, t_type, 'other models detailed',sep = '/')
        for (fn_num in fn_nums){
            edge_fnm = paste0(om_dir,'/fn_',as.character(fn_num),'_largest_connected_graph_edges.csv')
            node_fnm = paste0(om_dir,'/fn_',as.character(fn_num),'_largest_connected_graph_nodes.csv')
            edge_df = read.delim(edge_fnm, sep = ',')
            node_df = read.delim(node_fnm, sep = ',')

            g <- graph_from_data_frame(edge_df, directed=FALSE, vertices=node_df)
            # print(g, e=TRUE, v=TRUE)

            score0 = node_df$ini_score
            names(score0) = node_df$node

            ### save the output of the diffusion model
            ### comparing scores across methods, classical diffusion kernel
            for (sigma_2 in sigma_2s){
                K_d <- diffuStats::diffusionKernel(g, sigma2 = sigma_2)
                list_methods <- c("raw", "ber_s")
                df_diff <- diffuse_grid(
                    K = K_d,
                    scores = score0,
                    grid_param = expand.grid(method = list_methods),
                    n.perm = n_perm
                )

                raw = df_diff[df_diff$method=='raw',]
                ber_s = df_diff[df_diff$method=='ber_s',]

                raw = raw[order(raw$node_score, decreasing = TRUE),]
                ber_s = ber_s[order(ber_s$node_score, decreasing = TRUE),]

                rownames(raw) <- 1:nrow(raw)    # reset index after sorting
                rownames(ber_s) <- 1:nrow(ber_s)    # reset index after sorting

                raw$rank = row.names(raw)
                ber_s$rank = row.names(ber_s)

                model_tp = 'diffusion_model'
                method1 = list_methods[1] # raw model
                save_fnm = paste0(om_dir,'/',model_tp,'_wes+rna_fn_',as.character(fn_num),'_',method1,'_sigma2=',as.character(sigma_2),'.csv')
                write.table(raw,save_fnm,sep = ',',row.names = FALSE)

                method1 = list_methods[2] # ber_s model
                save_fnm = paste0(om_dir,'/',model_tp,'_wes+rna_fn_',as.character(fn_num),'_',method1,'_sigma2=',as.character(sigma_2),'.csv')
                write.table(ber_s,save_fnm,sep = ',',row.names = FALSE)
            }
        }
    }
}
diff_models(save_dir,t_types)
