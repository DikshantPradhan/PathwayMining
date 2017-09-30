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
yeast <- GRB_yeast_model()
n <- yeast$get_sizes()$NumVars
vars <- yeast$get_names()$VarName

yeast_og_set_list <- GRB_generate_set_list(yeast)

yeast <- GRB_yeast_model()
yeast_set_lists <- GRB_generate_set_lists(yeast, 1:n)

yeast_composition_set_full <- return_composition_sets(yeast_og_set_list, yeast_set_lists)

yeast_composition_set <- yeast_composition_set_full$composition

for (i in 1:length(yeast_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_og_set_list[[i]])
  if (length(yeast_og_set_list[[i]]) > 1){
  for (j in 1:(length(yeast_og_set_list[[i]])-1)){
    rxn1 <- GRB_get_rxn_idx(yeast, yeast_og_set_list[[i]][[j]])
    rxn2 <- GRB_get_rxn_idx(yeast, yeast_og_set_list[[i]][[j+1]])
    #print(c(rxn1, rxn2))
    if (length(yeast_composition_set[[rxn1]]) != length(yeast_composition_set[[rxn2]])){
      print("composition error")
      print(c(yeast_og_set_list[[i]][[j]], ";", yeast_og_set_list[[i]][[j+1]]))
      print(c(yeast_composition_set[[rxn1]], ";", yeast_composition_set[[rxn2]]))
    }
    #print(j)
    #print(ecoli_composition_set[[GRB_get_rxn_idx(ecoli, j)]])
  }
  }
}

print(yeast_composition_set_full$error)

yeast_deletion_list <- check_set_list_for_deletion(vars, yeast_set_lists)

for (i in 1:length(yeast_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_og_set_list[[i]])
  if (length(yeast_og_set_list[[i]]) > 1){
  for (j in 1:(length(yeast_og_set_list[[i]])-1)){
    rxn1 <- GRB_get_rxn_idx(yeast, yeast_og_set_list[[i]][[j]])
    rxn2 <- GRB_get_rxn_idx(yeast, yeast_og_set_list[[i]][[j+1]])
    #print(c(rxn1, rxn2))
    if (length(yeast_deletion_list[[rxn1]]) != length(yeast_deletion_list[[rxn2]])){
      print("deletion error")
      print(c(yeast_og_set_list[[i]][[j]], ";", yeast_og_set_list[[i]][[j+1]]))
      print(c(yeast_deletion_list[[rxn1]], ";", yeast_deletion_list[[rxn2]]))
    }
    #print(j)
    #print(ecoli_composition_set[[GRB_get_rxn_idx(ecoli, j)]])
  }
  }
}
