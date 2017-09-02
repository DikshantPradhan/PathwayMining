# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')

for (i in 1:length(og_set_list)){ # print sets joined by each deletion
  print(og_set_list[[i]])
  for (j in og_set_list[[i]]){
    print(j)
    print(test_composition_set[[get_rxn_idx(j)]])
  }
}
# 
# check_composing_sets <- function(idx){
#   print(composition_set[[idx]])
#   for (set in unlist(composition_set[[idx]])){
#     print(paste(set, ":"))
#     print(og_set_list[[set]])
#   }
#   print(set_lists[[idx]])
# }

# for (i in 1:length(og_set_list)){
#   print(og_set_list[[i]])
#   for (j in og_set_list[[i]]){
#     print(paste(j, "~~"))
#     for (k in og_set_list[[i]]){
#       # print(paste(k, ":"))
#       print(check_sets_for_containing(k, set_lists[[get_rxn_idx(j)]]))
#     }
#   }
# }

# for (i in unlist(set_lists[[17]])){
#   print(check_sets_for_containing(i, test_set_eno))
# }

# set_lists_2 <- c()
# for (i in 1:94){
#   set_lists_2[[i]] <- generate_set_lists(i)
# }
# 
# set_list_temp <- set_lists_2
# 
# for (i in 1:94){
#   set_lists_2[[i]] <- set_list_temp[[i]][[i]]
# }