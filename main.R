# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')



#model <- model_og$copy()
r0_coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled

GRB_set_list_mtx <- function(){
  
  # load model, r0_mtx as rdata file, make copies
  # params: suppression indices, reaction test indices
  
  print(paste('suppression index: ', i))
  if (!(r0_coupling_mtx[i,i])){
    print(paste(vars[i], ' blocked'))
    next
  }
  
  #prev_ub <- model$getattr("UB")[vars[i]]
  #prev_lb <- model$getattr("LB")[vars[i]]
  
  model <- model_og$copy() #GRB_ecoli_model()
  
  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))
  
  coupling_array[,,i] <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled
  
  # unfix i
  #model$setattr("UB", prev_ub)
  #model$setattr("LB", prev_lb)
  
}

lapply(suppr_idxs,
       FUN = GRB_set_list_mtx
)