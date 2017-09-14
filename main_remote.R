# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')

#source('~/GitHub/PathwayMining/load_mod.R')


#for (i in 1:length(ecoli_set_lists)){
#  for (j in 1:length(ecoli_set_lists[[i]])){
#    check <- check_sets_for_containing(ecoli_set_lists[[i]][[j]], set_lists[[i]])
#    if (!check){
#      print(c(i,j))
#      print(ecoli_set_lists[[i]][[j]])
#    }
#  }
#}

#for (i in 1:length(ecoli_r1_set)){
#  check <- check_sets_for_containing(union_set_list[[i]], ecoli_r1_set)
#  if (!check){print(c(i))}
#  else{print("good")}
#}

#for (i in 1:length(ecoli_set_lists_2)){
#  print(length(ecoli_set_lists[[i]]))
#}

for (i in 1:length(ecoli_pair_lists)){
  print(c(length(ecoli_pair_lists[[i]]), length(ecoli_pair_lists_2[[i]])))
  #for (j in 1:length(ecoli_pair_lists[[i]])){
  #  if (!all.equal(ecoli_pair_lists_2[[i]][j,],ecoli_pair_lists_2[[i]][j,])){print(c(i,j))}
  #}
}
