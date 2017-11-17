# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/falcon_tools.R')

ptm <- proc.time()

# # test composition of set_lists (make sure that blocking any reaction in an og_set results in the same r1 set)

ecoli <- GRB_ecoli_falcon_model()
n <- ecoli$get_sizes()$NumVars
vars <- ecoli$get_names()$VarName

reaction_indexes <- c()
reaction_indexes <- grep('Ex_a', vars)

ecoli_falcon_og_set_list <- GRB_generate_set_list(ecoli, reaction_indexes = reaction_indexes)

ecoli <- GRB_ecoli_falcon_model()
ecoli_set_lists <- GRB_generate_set_lists(ecoli, ecoli_og_set_list, 1:n, reaction_indexes)

ecoli_composition_set_full <- return_composition_sets(ecoli_og_set_list, ecoli_set_lists, ecoli)

ecoli_composition_set <- ecoli_composition_set_full$composition

print('composition error check')

# check for errors in compositions of r1 sets (intermediate) (each blockage within the same set should have the same effects)
for (i in 1:length(ecoli_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_og_set_list[[i]])
  if (length(ecoli_og_set_list[[i]]) > 1){
    for (j in 1:(length(ecoli_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, ecoli_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, ecoli_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      #print(paste(ecoli_composition_set[[rxn1]], ecoli_composition_set[[rxn2]]))
      if (rxn1 > length(ecoli_composition_set) | rxn2 > length(ecoli_composition_set)){next}

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

print(ecoli_composition_set_full$error)

print('deletion error check')

ecoli_deletion_list <- check_set_list_for_deletion(vars, ecoli_set_lists)

# check for inconsistencies in deletions of sets (each blockage within the same set should have the same effects)
for (i in 1:length(ecoli_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_og_set_list[[i]])
  if (length(ecoli_og_set_list[[i]]) > 1){
    for (j in 1:(length(ecoli_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, ecoli_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, ecoli_og_set_list[[i]][[j+1]])
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

save(ecoli_og_set_list, ecoli_set_lists, ecoli_composition_set_full, ecoli_deletion_list, file = "ecoli_run_data.RData")

proc.time() - ptm

print('specific error')

ct <- 0
blocked <- c()
for (i in 1:length(ecoli_composition_set_full$error)){
  if (nchar(ecoli_composition_set_full$error[i]) > 8){
    ct <- ct + 1
    #blocked <- c(blocked, ecoli_vars[i])
    print(ecoli_composition_set_full$error[i])
    #print(ecoli_composition_set_full$error[i])
  }
}

print(ct)

ecoli_r0_pairs <- return_pairs_from_set_list(ecoli_og_set_list)
ecoli_new_r1_pairs <- new_pairs_from_composition(ecoli_og_set_list, ecoli_composition_set)
ecoli_r1_pairs <- append_pair_lists(ecoli_r0_pairs, ecoli_new_r1_pairs)

ecoli_r0_set_list <- ecoli_og_set_list
ecoli_r1_set_list <- get_list_of_sets(ecoli_r1_pairs)

print('FIN')
