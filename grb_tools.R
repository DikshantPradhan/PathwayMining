

GRB_generate_set_lists <- function(model, suppression_idxs){
  set_lists <- c()
  
  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  
  for (i in suppression_idxs){
    
    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]
    
    # block i
    model$setattr("UB", setNames(0, vars[i]))
    model$setattr("LB", setNames(0, vars[i]))
    
    set_lists[i] <- list(get_list_of_sets(return_couples(flux_coupling_raptor(model)$coupled)))
    
    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)
    
  }
  
  return(set_lists)
}