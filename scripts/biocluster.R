source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
#source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
# source('~/GitHub/PathwayMining/logic_tools.R')
# source('~/GitHub/PathwayMining/falcon_tools.R')
#source('~/GitHub/PathwayMining/load_mod.R')
#source('~/GitHub/PathwayMining/data_tools.R')
#source('~/GitHub/PathwayMining/gene_tools.R')

## script for biocluster

GRB_flux_coupling_raptor_wrapper <- function(i, vars, model_og, reaction_indexes = 1:length(vars), compare_mtx = FALSE, r0_coupling_mtx = c()){
  print(paste('suppression index:', i))

  #if (!(r0_coupling_mtx[i,i])){
  #  print(paste(vars[i], ' blocked'))
  #  next
  #}

  #prev_ub <- model$getattr("UB")[vars[i]]
  #prev_lb <- model$getattr("LB")[vars[i]]

  model <- model_og$copy()

  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))

  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_mtx, known_set_mtx = r0_coupling_mtx)$coupled
  coupling_list <- Matrix(data = FALSE, nrow = 1, ncol = length(coupling_mtx))
  coupling_list[which(coupling_mtx)] <- TRUE

  #output <- list(which(coupling_mtx)) #list(get_list_of_sets(return_couples(coupling_mtx)))
  return(coupling_list)
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
  suppr_vector <- Matrix(data = FALSE, nrow = 1, ncol = n, sparse = TRUE)
  suppr_vector[suppression_idxs] <- TRUE
  if (optimize_suppr){
    i <- 1
    while (i <= n){
      if (suppr_vector[i]){ # if tagged to be suppressed
        set_idx <- which(r0_coupling_mtx[,i])[1] # which is first reaction (row) i is coupled to
        if (!is.na(set_idx)){
          rxn_idxs <- which(r0_coupling_mtx[set_idx,]) # other reactions in set
          # only suppress first reaction in set since, theoretically, suppressing any should have the same effect
          suppr_vector[rxn_idxs] <- FALSE
          suppr_vector[rxn_idxs[1]] <- TRUE
        }
        else {
          suppr_vector[i] <- FALSE
        }

      }
      i <- i+1
    }
  }

  print(paste("# of suppressions:", length(which(suppr_vector)), sep = " "))
  coupling <- lapply(which(suppr_vector), function(x) GRB_flux_coupling_raptor_wrapper(x, vars, model_og, reaction_indexes = reaction_indexes, compare_mtx = compare_known_r0_sets, r0_coupling_mtx = r0_coupling_mtx))

  return(coupling)
}


#ecoli <- GRB_ecoli_model()
#n <- ecoli$get_sizes()$NumVars
#vars <- ecoli$get_names()$VarName

#ecoli_r0_set_list <- GRB_generate_set_list(ecoli)

#proc.time() - ptm # timing end

#ecoli <- GRB_ecoli_model()
#ecoli_set_lists <- GRB_generate_set_lists(ecoli, 1:n, compare_known_r0_sets = TRUE, optimize_suppr=TRUE) #GRB_generate_set_lists(ecoli, ecoli_r0_set_list, 1:n)


ptm <- proc.time() # timing start
#load('maranas_falcon_model.RData')
#sybil_falcon_model <- yeast_falcon_model # load sybil falcon model
model <- GRB_yeast_falcon_model() #GRB_generate_falcon_model(sybil_falcon_model, falcon_model = TRUE)
n <- model$get_sizes()$NumVars
vars <- model$get_names()$VarName

load('~/GitHub/PathwayMining/data/yeast_model/Maranas_model/maranas_model_lipid_exch.RData')

sybil_model <- yeast_model
non_gene_assc_rxns <- which(sybil_model@genes == "")

gene_indexes <- c()
gene_indexes <- grep('Ex_a', vars)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)
#suppr_indexes <- suppr_indexes[1:150]

#reaction_indexes <- c()
#reaction_indexes <- grep('Ex_a', vars)

yeast_r0_coupling_mtx <- flux_coupling_raptor(model$copy(), reaction_indexes = gene_indexes)$coupled

save(yeast_r0_coupling_mtx, file = 'yeast_r0_coupling_mtx.RData')

set_lists <- GRB_generate_set_lists_cluster(model, suppression_idxs = suppr_indexes, reaction_indexes = suppr_indexes, compare_known_r0_sets = TRUE, optimize_suppr=TRUE)
proc.time() - ptm # timing end
#print(compare_sets(ecoli_r1_sets, ecoli_r1_sets_cluster))
