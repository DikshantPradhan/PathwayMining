# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/data_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')

# lethal_pairs <- pairs[lethal_double_dels,]
# lethal_parts <- unique(lethal_pairs[,1], lethal_pairs[,2])
# lethal_part_g0_set_idxs <- c()
# for (part in lethal_parts){
#   lethal_part_g0_set_idxs <- c(get_set_idx(part, clean_mutans_g0_set), lethal_part_g0_set_idxs)
# }
# 
# for (set in lethal_part_set_idxs){
#   print(length(clean_mutans_g0_set[[set]]))
# }
# 
# 
# 
# lethal_part_g1_set_idxs <- c()
# for (part in lethal_parts){
#   lethal_part_g1_set_idxs <- c(get_set_idx(part, clean_mutans_g1_set), lethal_part_g1_set_idxs)
# }
# for (set in lethal_part_g1_set_idxs){
#   print(length(clean_mutans_g1_set[[set]]))
# }



# output <- get_path_between_reactions(mutans_falcon@S, 845, 846)