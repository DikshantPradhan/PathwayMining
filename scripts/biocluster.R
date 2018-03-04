source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/data_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')

## script for biocluster

GRB_flux_coupling_raptor_wrapper <- function(i, vars, model_og, reaction_indexes = 1:length(vars), compare_mtx = FALSE, r0_coupling_mtx = c()){
  print(paste('suppression index:', i))

  #if (!(r0_coupling_mtx[i,i])){
  #  print(paste(vars[i], ' blocked'))
  #  next
  #}

  #prev_ub <- model$getattr("UB")[vars[i]]
  #prev_lb <- model$getattr("LB")[vars[i]]

  model <- model_og$copy() #GRB_ecoli_model()

  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))

  output <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_mtx, known_set_mtx = r0_coupling_mtx)$coupled

  return(output)
}

GRB_generate_set_lists_cluster <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
                                         compare_known_r0_sets = FALSE, optimize_suppr = FALSE){

  n <- model_og$get_sizes()$NumVars
  vars <- model_og$get_names()$VarName

  if (suppression_idxs == -1){
    if (length(reaction_indexes) > 0){
      suppression_idxs = reaction_indexes
    }
    else {
      suppression_idxs = 1:n
    }
  }

  # dim: rxns_row, rxns_col, deletions
  model <- model_og$copy()

  r0_coupling_mtx <- c()
  if (compare_known_r0_sets){
    r0_coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled
    r0_coupling_mtx <- fill_coupling_matrix(r0_coupling_mtx)
  }

  suppr_vector <- matrix(data = FALSE, nrow = 1, ncol = n)
  suppr_vector[suppression_idxs] <- TRUE
  if (optimize_suppr){
    i <- 1
    while (i <= n){
      if (suppr_vector[i]){ # if tagged to be suppressed
        set_idx <- which(r0_coupling_mtx[,i])[1] # which is first reaction (row) i is coupled to
        rxn_idxs <- which(r0_coupling_mtx[set_idx,]) # other reactions in set
        # only suppress first reaction in set since, theoretically, suppressing any should have the same effect
        suppr_vector[rxn_idxs] <- FALSE
        suppr_vector[rxn_idxs[1]] <- TRUE
      }
      i <- i+1
    }
  }

  # coupling_array <- array(data = FALSE, dim = c(n,n,n), dimnames = list(vars, vars, paste('del', vars, sep = "_")))
  print(paste("# of suppressions:", length(which(suppr_vector)), sep = " "))

  coupling <- lapply(which(suppr_vector), function(x) GRB_flux_coupling_raptor_wrapper(x, vars, model_og, reaction_indexes = reaction_indexes, compare_mtx = compare_known_r0_sets, r0_coupling_mtx = r0_coupling_mtx))

  return(coupling)
}

convert_coupling_list_to_array <- function(coupling_list){
  coupling_matrix <- coupling_list[[1]]
  vars <- rownames(coupling_matrix)
  coupling_array <- array(data = FALSE, dim = c(nrow(coupling_matrix),ncol(coupling_matrix),length(coupling_list)), dimnames = list(vars, vars, seq(from = 1, to = length(coupling_list))))
  for (i in 1:length(coupling_list)){
    coupling_array[,,i] <- coupling_list[[i]]
  }
  return(coupling_array)
}

#cluster_func <- function(){
#  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes,
#                                              compare_mtx = compare_known_r0_sets, known_set_mtx = r0_coupling_mtx)$coupled
#
#  return(coupling_mtx)
#}

model <- GRB_ecoli_model()
ecoli_coupling_list <- GRB_generate_set_lists_cluster(model, suppression_idxs = c(1,5,9), compare_known_r0_sets = TRUE, optimize_suppr = TRUE)
ecoli_coupling_array <- convert_coupling_list_to_array(ecoli_coupling_list)

ecoli_r1_matrix <- coupling_matrix_from_array(ecoli_coupling_array)
ecoli_r1_matrix <- (ecoli_r1_matrix > 0)
ecoli_r1_sets_cluster <- list(get_list_of_sets(return_couples(ecoli_r1_matrix)))

ecoli <- GRB_ecoli_model()
n <- ecoli$get_sizes()$NumVars
vars <- ecoli$get_names()$VarName
ecoli_coupling_array <- GRB_generate_set_lists_array(ecoli, suppression_idxs = c(1,5,9), compare_known_r0_sets = TRUE, optimize_suppr=TRUE)
ecoli_r1_matrix <- coupling_matrix_from_array(ecoli_coupling_array)
ecoli_r1_matrix <- (ecoli_r1_matrix > 0)
ecoli_r1_sets <- list(get_list_of_sets(return_couples(ecoli_r1_matrix)))

print(compare_sets(ecoli_r1_sets, ecoli_r1_sets_cluster))
