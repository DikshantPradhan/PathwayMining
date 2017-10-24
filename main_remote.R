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
yeast <- GRB_yeast_model()
n <- yeast$get_sizes()$NumVars
vars <- yeast$get_names()$VarName

yeast_og_set_list <- GRB_generate_set_list(yeast)

yeast <- GRB_yeast_model()
yeast_set_lists <- GRB_generate_set_lists(yeast, yeast_og_set_list, 1:n)

yeast_composition_set_full <- return_composition_sets(yeast_og_set_list, yeast_set_lists, yeast)
#
yeast_composition_set <- yeast_composition_set_full$composition
#
for (i in 1:length(yeast_og_set_list)){ # print sets joined by each deletion
  #print(yeast_og_set_list[[i]])
  if (length(yeast_og_set_list[[i]]) > 1){
  for (j in 1:(length(yeast_og_set_list[[i]])-1)){
    rxn1 <- get_rxn_idx(vars, yeast_og_set_list[[i]][[j]])
    rxn2 <- get_rxn_idx(vars, yeast_og_set_list[[i]][[j+1]])
    #print(c(rxn1, rxn2))
    #print(paste(yeast_composition_set[[rxn1]], yeast_composition_set[[rxn2]]))
    if (rxn1 > length(yeast_composition_set) | rxn2 > length(yeast_composition_set)){next}

    if (length(yeast_composition_set[[rxn1]]) != length(yeast_composition_set[[rxn2]])){
      print("composition error")
      print(c(yeast_og_set_list[[i]][[j]], ";", yeast_og_set_list[[i]][[j+1]]))
      print(c(yeast_composition_set[[rxn1]], ";", yeast_composition_set[[rxn2]]))
    }
    #print(j)
    #print(yeast_composition_set[[GRB_get_rxn_idx(yeast, j)]])
  }
  }
}

print(yeast_composition_set_full$error)

yeast_deletion_list <- check_set_list_for_deletion(vars, yeast_set_lists)

for (i in 1:length(yeast_og_set_list)){ # print sets joined by each deletion
  #print(yeast_og_set_list[[i]])
  if (length(yeast_og_set_list[[i]]) > 1){
  for (j in 1:(length(yeast_og_set_list[[i]])-1)){
    rxn1 <- get_rxn_idx(vars, yeast_og_set_list[[i]][[j]])
    rxn2 <- get_rxn_idx(vars, yeast_og_set_list[[i]][[j+1]])
    #print(c(rxn1, rxn2))
    if (length(yeast_deletion_list[[rxn1]]) != length(yeast_deletion_list[[rxn2]])){
      print("deletion error")
      print(c(yeast_og_set_list[[i]][[j]], ";", yeast_og_set_list[[i]][[j+1]]))
      print(c(yeast_deletion_list[[rxn1]], ";", yeast_deletion_list[[rxn2]]))
    }
    #print(j)
    #print(yeast_composition_set[[GRB_get_rxn_idx(yeast, j)]])
  }
  }
}

save(yeast_og_set_list, yeast_set_lists, yeast_composition_set_full, yeast_deletion_list, file = "yeast_run_data.RData")

proc.time() - ptm

ct <- 0
blocked <- c()
for (i in 1:length(yeast_composition_set_full$error)){
  if (nchar(yeast_composition_set_full$error[i]) > 8){
    ct <- ct + 1
    #blocked <- c(blocked, yeast_vars[i])
    print(yeast_composition_set_full$error[i])
    #print(yeast_composition_set_full$error[i])
  }
}

print(ct)
