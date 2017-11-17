source('falcon_tools.R')

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
  #yeast_model <- readTSVmod(reactList = "Y7_test_react.tsv", metList = "Y7_met.tsv")
  #yeast_model <- rmReact(model = yeast_model, react = 1606)
  #yeast_model <- rmReact(model = yeast_model, react = 1590)

  #setwd("~/GitHub/PathwayMining/")

  #yeast_model <- get_yeast_model()

  load('yeast_open_mod.RData')

  for (i in 1:length(yeast_open_mod@react_id)){
    yeast_open_mod@lowbnd[i] <- -1000
    yeast_open_mod@uppbnd[i] <- 1000
  }

  #yeast_model@lowbnd[223] <- 0

  setwd("~/GitHub/PathwayMining/")

  yeast <- as_GRBmodel(yeast_open_mod)
  yeast$show_output(FALSE)

  return(yeast)
}

GRB_mutans_model <- function(){

  #setwd("~/GitHub/PathwayMining/data/mutans_model")

  #setwd("~/GitHub/PathwayMining/")

  #load('mutans_model.RData')
  load('data/mutans_model/mutans_model.RData')

  for (i in 1:length(mutans@react_id)){
    mutans@lowbnd[i] <- -1000
    mutans@uppbnd[i] <- 1000
  }

  setwd("~/GitHub/PathwayMining/")

  mutans <- as_GRBmodel(mutans)
  mutans$show_output(FALSE)

  return(mutans)
}


GRB_generate_falcon_model <- function(sybil_model, r0_gene_set = c(), r0_rxn_set_list = c()){
  sybil_falcon_model <- generate_falcon_model(sybil_model, r0_gene_set, r0_rxn_set_list)

  grb_falcon_model <- as_GRBmodel(sybil_falcon_model)
  grb_falcon_model$show_output(FALSE)

  ## ADD NECESSARY CONSTRAINTS TO MODEL
  vars <- grb_falcon_model$get_names()$VarName

  split_fwd_rxns <- vars[grep('fwd', vars)]
  split_rev_rxns <- vars[grep('rev', vars)]
  #print(split_rxns)
  #conv_num <- sapply(split_rxns, function(x) strsplit(x, ' |_')[[1]][5])
  #dir <- sapply(split_rxns, function(x) strsplit(x, ' |_')[[1]][3])
  #split_rxns <- sapply(split_rxns, function(x) strsplit(x, ' |_')[[1]][2])
  #split_rxns <- unique(split_rxns)

  print(split_fwd_rxns)
  print(split_rev_rxns)
  #print(conv_num)
  #print(dir)

  #print(grb_falcon_model)
  if (length(split_fwd_rxns) != length(split_rev_rxns)){
    print('fed rev matchup error')
    return()
  }

  # add bounds on split conversion reactions (gene -> activity_[rxn])
  for (i in 1:length(split_fed_rxns)){
    fwd <- split_fwd_rxns[i]
    rev <- split_rev_rxns[i]

    # get bounds (care about fwd_ub & rev_lb)
    fwd_ub <- grb_falcon_model$getattr("UB")[fwd]
    fwd_lb <- grb_falcon_model$getattr("LB")[fwd]
    rev_ub <- grb_falcon_model$getattr("UB")[rev]
    rev_lb <- grb_falcon_model$getattr("LB")[rev]

    a_rxn <- strsplit(fwd, ' ')[[1]][1]
    I <- paste('I', a_rxn, sep = '_')

    #grb_falcon_model$addvar(name = I, vtype = 'B') #??? not sure how to add binary constraint
    #grb_falcon_model$addconstr(paste(fwd, '*', I, sep = ''),
        #sense="<=", rhs= fwd_ub, name = paste(a_rxn, 'fwd', sep = '_')) # bound on fwd conversion
    #grb_falcon_model$addconstr(paste(rev, '*(1 - ', I, ')',  sep = ''),
        #sense=">=", rhs= rev_lb, name = paste(a_rxn, 'rev', sep = '_')) # bound on rev conversion
  }

  # for each reaction w fwd and rev components:
  #   add constraint so that only one can run at a time
  #   binary vars I_fwd, I_rev <- {0, 1} and I_fwd + I_rev = 1

  return(grb_falcon_model)
}

GRB_get_rxn_idx <- function(model, rxn){
  vars <- model$get_names()$VarName
  return(return(which(vars == rxn)))
}

GRB_generate_pair_list <- function(model_og){
  model <- model_og$copy()
  return(return_couples(flux_coupling_raptor(model)$coupled))
}

GRB_generate_set_list <- function(model_og, reaction_indexes = c()){
  model <- model_og$copy()
  return(get_list_of_sets(return_couples(flux_coupling_raptor(model,
                                reaction_indexes = reaction_indexes)$coupled)))
}

GRB_generate_pair_lists <- function(model_og, suppression_idxs){
  pair_lists <- c()

  n <- model_og$get_sizes()$NumVars
  vars <- model_og$get_names()$VarName

  for (i in suppression_idxs){

    #prev_ub <- model$getattr("UB")[vars[i]]
    #prev_lb <- model$getattr("LB")[vars[i]]

    model <- model_og$copy() #GRB_ecoli_model()

    # block i
    model$setattr("UB", setNames(0, vars[i]))
    model$setattr("LB", setNames(0, vars[i]))

    pair_lists[i] <- list(return_couples(flux_coupling_raptor(model)$coupled))
    #pair_lists[i] <- list(return_couples(raptor::flux_coupling(model)$coupled))

    # unfix i
    #model$setattr("UB", prev_ub)
    #model$setattr("LB", prev_lb)
  }
  return(pair_lists)
}

GRB_generate_set_lists <- function(model_og, og_set_list, suppression_idxs, reaction_indexes = c()){
  set_lists <- c()
  unblocked_rxns <- unlist(og_set_list)

  n <- model_og$get_sizes()$NumVars
  vars <- model_og$get_names()$VarName

  for (i in suppression_idxs){
    print(i)
    if (!(vars[i] %in% unblocked_rxns)){
      print(paste(vars[i], ' blocked'))
      next
    }

    #prev_ub <- model$getattr("UB")[vars[i]]
    #prev_lb <- model$getattr("LB")[vars[i]]

    model <- model_og$copy() #GRB_ecoli_model()

    # block i
    model$setattr("UB", setNames(0, vars[i]))
    model$setattr("LB", setNames(0, vars[i]))

    set_lists[i] <- list(get_list_of_sets(return_couples(flux_coupling_raptor(model,
                            reaction_indexes = reaction_indexes)$coupled)))

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

GRB_get_blocked <- function(model){
  lb <- which(model$getattr('LB') > -1000)
  ub <- which(model$getattr('UB') < 1000)

  potential_blocked <- intersect(lb, ub)
  blocked <- c()
  for (i in potential_blocked){
    if (model$getattr('LB')[i] == 0 & model$getattr('UB')[i] == 0){
      blocked <- c(blocked, i)
    }
  }

  return(model$get_names()$VarName[blocked])
}
