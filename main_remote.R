# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
# source('~/GitHub/PathwayMining/raptor_coupling.R')
# source('~/GitHub/PathwayMining/grb_tools.R')

# source('~/GitHub/PathwayMining/load_mod.R')

#rxns <- c()

# for (i in 1:length(ecoli_r1_set)){
#   # test_set <- get_list_of_sets(ecoli_pair_lists[[i]])
#
#   check <- check_sets_for_containing(union_set_list_2[[i]], ecoli_r1_set)
#   if (!check){
#     print(c(i))
#   }
#   # else{print("good")}
#
#   # for (j in 1:length(test_set)){
#   #   check <- check_sets_for_containing(test_set[[j]], set_lists[[i]])
#   #   if (!check){
#   #     print(c(i,j))
#   #     print(test_set[[j]])
#   #     # rxns <- c(rxns, unlist(test_set[[j]]))
#   #   }
#   # }
# }

#for (i in 1:length(ecoli_og_set_list)){
#
#  for (j in 1:length(ecoli_og_set_list[[i]])){
#    check <- check_sets_for_containing(ecoli_og_set_list[[i]][[j]], set_lists[[i]])
#    if (!check){
#      print(c(i,j))
#      print(ecoli_og_set_list[[i]][[j]])
#      # rxns <- c(rxns, unlist(test_set[[j]]))
#    }
#  }
#}

#for (i in 1:length(ecoli_r1_set_test)){
#  check <- check_sets_for_containing(ecoli_r1_set_test[[i]], union_set_list)
#  if (!check){print(c(i))}
#  #else{print("good")}
#}

#for (i in 1:length(ecoli_set_lists_2)){
#  print(length(ecoli_set_lists[[i]]))
#}

#for (i in 1:length(ecoli_pair_lists)){
  #print(c(length(ecoli_pair_lists[[i]]), length(ecoli_pair_lists_2[[i]])))
#  test_set <- get_list_of_sets(ecoli_pair_lists[[i]])
#  print(c(length(test_set), length(ecoli_set_lists[[i]])))
  #for (j in 1:length(test_set)){
  #  print(c(length(test_set[[j]]), length(ecoli_set_lists[[i]][[j]])))
  #  #for (k in 1:length(test_set[[i]])){
  #  #  if (!all.equal(test_set[[i]][j], ecoli_set_lists[[i]][j])){print(c(i,j))}
  #  #}
  #}
  #for (j in 1:length(ecoli_pair_lists[[i]])){
  #  if (!all.equal(ecoli_pair_lists_2[[i]][j,],ecoli_pair_lists_2[[i]][j,])){print(c(i,j))}
  #}
#}

for (i in 1:94){
  #print(c(length(ecoli_set_lists_1[[i]]), length(ecoli_set_lists_2[[i]])))
  for (j in 1:length(ecoli_set_lists_1[[i]])){
    check <- check_sets_for_containing(ecoli_set_lists_1[[i]][[j]], ecoli_set_lists_2[[i]])
    if (!check){print(c(i,j)); print(ecoli_set_lists_1[[i]][[j]])}
  }
}
