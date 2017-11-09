# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
# source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')

# source('~/GitHub/PathwayMining/load_mod.R')

ptm <- proc.time()

# # test composition of set_lists (make sure that blocking any reaction in an og_set results in the same r1 set)

mutans <- GRB_mutans_model()
n <- mutans$get_sizes()$NumVars
vars <- mutans$get_names()$VarName

mutans_og_set_list <- GRB_generate_set_list(mutans)

mutans <- GRB_mutans_model()
mutans_set_lists <- GRB_generate_set_lists(mutans, mutans_og_set_list, 1:n)

mutans_composition_set_full <- return_composition_sets(mutans_og_set_list, mutans_set_lists, mutans)

mutans_composition_set <- mutans_composition_set_full$composition

# check for errors in compositions of r1 sets (intermediate) (each blockage within the same set should have the same effects)
for (i in 1:length(mutans_og_set_list)){ # print sets joined by each deletion
  #print(mutans_og_set_list[[i]])
  if (length(mutans_og_set_list[[i]]) > 1){
    for (j in 1:(length(mutans_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, mutans_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, mutans_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      #print(paste(mutans_composition_set[[rxn1]], mutans_composition_set[[rxn2]]))
      if (rxn1 > length(mutans_composition_set) | rxn2 > length(mutans_composition_set)){next}

      if (length(mutans_composition_set[[rxn1]]) != length(mutans_composition_set[[rxn2]])){
        print("composition error")
        print(c(mutans_og_set_list[[i]][[j]], ";", mutans_og_set_list[[i]][[j+1]]))
        print(c(mutans_composition_set[[rxn1]], ";", mutans_composition_set[[rxn2]]))
      }
      #print(j)
      #print(mutans_composition_set[[GRB_get_rxn_idx(mutans, j)]])
    }
  }
}

print(mutans_composition_set_full$error)

mutans_deletion_list <- check_set_list_for_deletion(vars, mutans_set_lists)

# check for inconsistencies in deletions of sets (each blockage within the same set should have the same effects)
for (i in 1:length(mutans_og_set_list)){ # print sets joined by each deletion
  #print(mutans_og_set_list[[i]])
  if (length(mutans_og_set_list[[i]]) > 1){
    for (j in 1:(length(mutans_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, mutans_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, mutans_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      if (length(mutans_deletion_list[[rxn1]]) != length(mutans_deletion_list[[rxn2]])){
        print("deletion error")
        print(c(mutans_og_set_list[[i]][[j]], ";", mutans_og_set_list[[i]][[j+1]]))
        print(c(mutans_deletion_list[[rxn1]], ";", mutans_deletion_list[[rxn2]]))
      }
      #print(j)
      #print(mutans_composition_set[[GRB_get_rxn_idx(mutans, j)]])
    }
  }
}

save(mutans_og_set_list, mutans_set_lists, mutans_composition_set_full, mutans_deletion_list, file = "mutans_run_data.RData")

proc.time() - ptm

ct <- 0
blocked <- c()
for (i in 1:length(mutans_composition_set_full$error)){
  if (nchar(mutans_composition_set_full$error[i]) > 8){
    ct <- ct + 1
    #blocked <- c(blocked, mutans_vars[i])
    print(mutans_composition_set_full$error[i])
    #print(mutans_composition_set_full$error[i])
  }
}

print(ct)

mutans_r0_pairs <- return_pairs_from_set_list(mutans_og_set_list)
mutans_new_r1_pairs <- new_pairs_from_composition(mutans_og_set_list, mutans_composition_set)
mutans_r1_pairs <- append_pair_lists(mutans_r0_pairs, mutans_new_r1_pairs)

mutans_r0_set_list <- mutans_og_set_list
mutans_r1_set_list <- get_list_of_sets(mutans_r1_pairs)

print('FIN')
