GRB_ecoli_model <- function(){
  data(Ec_core);
  model=Ec_core;
  model <- changeBounds(model, 11, lb = 0) # this and next lead to no PGI in GND coset
  # model <- changeBounds(model, 13, lb = 0, ub = 0)
  model <- rmReact(model = model, react = 13)
  for (i in findExchReact(model)@react_pos){
    model <- changeBounds(model, i, lb = -1000, ub = 1000)
    # if (model@lowbnd[i] == 0){
    #   model <- changeBounds(model, i, lb = -1000)
    # }
  }

  ecoli <- as_GRBmodel(model)
  ecoli$show_output(FALSE)
  return(ecoli)
}

GRB_yeast_model <- function(){

  setwd("~/GitHub/PathwayMining/data/yeast_model")
  yeast_model <- readTSVmod(reactList = "Y7_test_react.tsv", metList = "Y7_met.tsv")
  yeast_model <- rmReact(model = yeast_model, react = 1606)
  yeast_model <- rmReact(model = yeast_model, react = 1590)

  setwd("~/GitHub/PathwayMining/")

  yeast <- as_GRBmodel(yeast_model)
  yeast$show_output(FALSE)

  return(yeast)
}

GRB_get_rxn_idx <- function(model, rxn){
  vars <- model$get_names()$VarName
  return(return(which(vars == rxn)))
}

GRB_generate_pair_list <- function(model){
  return(return_couples(raptor::flux_coupling(model)$coupled))
}

GRB_generate_set_list <- function(model){
  return(get_list_of_sets(return_couples(flux_coupling_raptor(model, min_fva_cor=0.99)$coupled)))
}

GRB_generate_pair_lists <- function(model, suppression_idxs){
  pair_lists <- c()

  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName

  for (i in suppression_idxs){

    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]

    #model <- GRB_ecoli_model()

    # block i
    model$setattr("UB", setNames(0, vars[i]))
    model$setattr("LB", setNames(0, vars[i]))

    pair_lists[i] <- list(return_couples(flux_coupling_raptor(model)$coupled))
    #pair_lists[i] <- list(return_couples(raptor::flux_coupling(model)$coupled))

    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)
  }
  return(pair_lists)
}

GRB_generate_set_lists <- function(model, suppression_idxs){
  set_lists <- c()

  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName

  for (i in suppression_idxs){

    #prev_ub <- model$getattr("UB")[vars[i]]
    #prev_lb <- model$getattr("LB")[vars[i]]

    model <- GRB_ecoli_model()

    # block i
    model$setattr("UB", setNames(0, vars[i]))
    model$setattr("LB", setNames(0, vars[i]))

    set_lists[i] <- list(get_list_of_sets(return_couples(flux_coupling_raptor(model, min_fva_cor=0.99)$coupled)))

    # unfix i
    #model$setattr("UB", prev_ub)
    #model$setattr("LB", prev_lb)

  }

  return(set_lists)
}

GRB_get_union_set_from_degen_pairs <- function(model, pair_lists){

  set_list = model$get_names()$VarName

  for (i in 1:length(pair_lists)){
    set_list <- get_list_of_sets(pair_lists[[i]], rxns_list = set_list)
    # print(set_list)
    print(paste(i, ": ", length(set_list)))
  }

  return(set_list)
}

GRB_r1_set <- function(model){
  return(GRB_get_union_set_from_degen_pairs(model, GRB_generate_pair_lists(model, 1:model$get_sizes()$NumVars)))
}
