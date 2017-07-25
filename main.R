source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')

for (i in 1:length(og_set_list)){
  print(og_set_list[[i]])
  for (j in og_set_list[[i]]){
    print(j)
    print(composition_set[[get_rxn_idx(j)]])
  }
}

check_composing_sets <- function(idx){
  print(composition_set[[idx]])
  for (set in unlist(composition_set[[idx]])){
    print(paste(set, ":"))
    print(og_set_list[[set]])
  }
  print(set_lists[[idx]])
}