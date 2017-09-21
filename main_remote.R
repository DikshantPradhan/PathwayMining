# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')

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

#for (i in 1:94){
#  #print(c(length(ecoli_set_lists_1[[i]]), length(ecoli_set_lists_2[[i]])))
#  for (j in 1:length(ecoli_set_lists_1[[i]])){
#    check <- check_sets_for_containing(ecoli_set_lists_1[[i]][[j]], ecoli_set_lists_2[[i]])
#    if (!check){print(c(i,j)); print(ecoli_set_lists_1[[i]][[j]])}
#  }
#}

# test composition of set_lists (make sure that blocking any reaction in an og_set results in the same r1 set)
ecoli <- GRB_ecoli_model()
vars <- ecoli$get_names()$VarName

ecoli_og_set_list <- GRB_generate_set_list(ecoli)

ecoli <- GRB_ecoli_model()
ecoli_set_lists <- GRB_generate_set_lists(ecoli, 1:94)

ecoli_composition_set <- return_composition_sets(ecoli_og_set_list, ecoli_set_lists)

for (i in 1:length(ecoli_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_og_set_list[[i]])
  if (length(ecoli_og_set_list[[i]]) > 1){
  for (j in 1:(length(ecoli_og_set_list[[i]])-1)){
    rxn1 <- GRB_get_rxn_idx(ecoli, ecoli_og_set_list[[i]][[j]])
    rxn2 <- GRB_get_rxn_idx(ecoli, ecoli_og_set_list[[i]][[j+1]])
    #print(c(rxn1, rxn2))
    if (length(ecoli_composition_set[[rxn1]]) != length(ecoli_composition_set[[rxn2]])){
      print("composition error")
      print(c(ecoli_og_set_list[[i]][[j]], ";", ecoli_og_set_list[[i]][[j+1]]))
      print(c(ecoli_composition_set[[rxn1]], ";", ecoli_composition_set[[rxn2]]))
    }
    #print(j)
    #print(ecoli_composition_set[[GRB_get_rxn_idx(ecoli, j)]])
  }
  }
}

ecoli_deletion_list <- check_set_list_for_deletion(vars, ecoli_set_lists)

for (i in 1:length(ecoli_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_og_set_list[[i]])
  if (length(ecoli_og_set_list[[i]]) > 1){
  for (j in 1:(length(ecoli_og_set_list[[i]])-1)){
    rxn1 <- GRB_get_rxn_idx(ecoli, ecoli_og_set_list[[i]][[j]])
    rxn2 <- GRB_get_rxn_idx(ecoli, ecoli_og_set_list[[i]][[j+1]])
    #print(c(rxn1, rxn2))
    if (length(ecoli_deletion_list[[rxn1]]) != length(ecoli_deletion_list[[rxn2]])){
      print("deletion error")
      print(c(ecoli_og_set_list[[i]][[j]], ";", ecoli_og_set_list[[i]][[j+1]]))
      print(c(ecoli_deletion_list[[rxn1]], ";", ecoli_deletion_list[[rxn2]]))
    }
    #print(j)
    #print(ecoli_composition_set[[GRB_get_rxn_idx(ecoli, j)]])
  }
  }
}
