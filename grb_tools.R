source('~/GitHub/PathwayMining/falcon_tools.R')

## MODELS

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

GRB_ecoli_falcon_model <- function(){
  sybil_ecoli <- get_ecoli_model()
  ecoli_falcon_model <- GRB_generate_falcon_model(sybil_ecoli)
  return(ecoli_falcon_model)
}

GRB_yeast_model <- function(){
  #maranas_model_exch_add_biom_rm/maranas_model_lipid_exch
  load('~/GitHub/PathwayMining/data/yeast_model/Maranas_model/yeast_model.RData')

  # setwd("~/GitHub/PathwayMining/")

  yeast <- as_GRBmodel(yeast_model)
  yeast$show_output(FALSE)

  return(yeast)
}

GRB_yeast_falcon_model <- function(){

  #load('~/GitHub/PathwayMining/data/yeast_model/Price Models/yeast_open_mod.RData')
  load('~/GitHub/PathwayMining/data/yeast_model/Maranas_model/yeast_model.RData')

  sybil_yeast <- yeast_model
  yeast_falcon_model <- GRB_generate_falcon_model(sybil_yeast)
  return(yeast_falcon_model)
}

GRB_mutans_model <- function(){

  #setwd("~/GitHub/PathwayMining/data/mutans_model")

  #setwd("~/GitHub/PathwayMining/")

  #load('mutans_model.RData')
  load('~/GitHub/PathwayMining/data/mutans_model/mutans_model.RData')
  # for (i in findExchReact(mutans)@react_pos){
  #   mutans <- changeBounds(mutans, i, lb = -1000, ub = 1000)
  # }
  #for (i in 1:length(mutans@react_id)){
  #  mutans@lowbnd[i] <- -1000
  #  mutans@uppbnd[i] <- 1000
  #}

  # setwd("~/GitHub/PathwayMining/")

  mutans <- as_GRBmodel(mutans)
  mutans$show_output(FALSE)

  return(mutans)
}

GRB_mutans_falcon_model <- function(){
  load('~/GitHub/PathwayMining/data/mutans_model/mutans_model.RData')

  sybil_mutans <- mutans
  mutans_falcon_model <- GRB_generate_falcon_model(sybil_mutans)
  return(mutans_falcon_model)
}

GRB_pao_model <- function(){
  load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')

  pao <- as_GRBmodel(pao_model)
  pao$show_output(FALSE)

  return(pao)
}

GRB_pao_falcon_model <- function(){
  load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')

  pao_falcon <- GRB_generate_falcon_model(pao_model)
  pao_falcon$show_output(FALSE)

  return(pao_falcon)
}

## FUNCTIONS

GRB_generate_falcon_model <- function(sybil_model, falcon_model = FALSE, r0_gene_set = c(), r0_rxn_set_list = c()){

  sybil_falcon_model <- sybil_model

  if (!falcon_model){ # if user has not passed in a falcon model already
    sybil_falcon_model <- generate_falcon_model(sybil_model, r0_gene_set, r0_rxn_set_list)
  }

  grb_falcon_model <- as_GRBmodel(sybil_falcon_model)
  grb_falcon_model$show_output(FALSE)

  ## ADD NECESSARY CONSTRAINTS TO MODEL
  vars <- grb_falcon_model$get_names()$VarName

  split_fwd_rxns <- vars[grep('fwd', vars)]
  split_rev_rxns <- vars[grep('rev', vars)]

  split_rxns <- sapply(split_fwd_rxns, function(x) strsplit(x, ' ')[[1]][1])
  split_rxns <- unique(split_rxns)

  if (length(split_fwd_rxns) != length(split_rev_rxns)){
    print('fed rev matchup error')
    return()
  }

  for (rxn in split_rxns){
    #a_rxn <- paste('a', rxn, sep = '_')
    I <- paste('I', rxn, sep = '_')
    .Call("GRB_addvar", grb_falcon_model$exptr, 0L, integer(0), numeric(0), 1.0, 0.0, 1.0, 'B', I)
  }

  .Call("GRB_updatemodel", grb_falcon_model$exptr)
  vars <- grb_falcon_model$get_names()$VarName
  n <- grb_falcon_model$get_sizes()$NumVars

  # add bounds on split conversion reactions (gene -> activity_[rxn])
  for (i in 1:length(split_fwd_rxns)){
    fwd <- split_fwd_rxns[i]
    rev <- split_rev_rxns[i]

    #print('fwd, rev')
    #print(fwd)
    #print(rev)

    # get bounds (care about fwd_ub & rev_lb)
    fwd_ub <- grb_falcon_model$getattr("UB")[[fwd]]
    fwd_lb <- grb_falcon_model$getattr("LB")[[fwd]]
    rev_ub <- grb_falcon_model$getattr("UB")[[rev]]
    rev_lb <- grb_falcon_model$getattr("LB")[[rev]]

    a_rxn <- strsplit(fwd, ' ')[[1]][1]
    I <- paste('I', a_rxn, sep = '_')
    i_idx <- which(vars == I)
    fwd_idx <- which(vars == fwd)
    rev_idx <- which(vars == rev)

    #print(paste(fwd, fwd_ub, fwd_lb, rev, rev_ub, rev_lb, I))
    fwd_vec <- c(fwd_ub, -1)
    rev_vec <- c(-1*rev_ub, -1)
    fwd_idxs <- c(i_idx, which(vars == fwd))
    rev_idxs <- c(i_idx, which(vars == rev))

    fwd_name <- paste(fwd, 'I', sep = ' ')
    rev_name <- paste(rev, 'I', sep = ' ')

    #.Call("GRB_updatemodel", grb_falcon_model$exptr) ##
    .Call("GRB_addconstr", grb_falcon_model$exptr, 2L, as.integer(fwd_idxs-1), fwd_vec, ">=", 0.0,
      fwd_name)
    #.Call("GRB_addconstr", grb_falcon_model$exptr, 2L, as.integer(rev_idxs-1), rev_vec, "<=", (-1*rev_lb),
    #  rev_name)
    .Call("GRB_addconstr", grb_falcon_model$exptr, 2L, as.integer(rev_idxs-1), rev_vec, ">=", (-1*rev_ub),
      rev_name)
    .Call("GRB_updatemodel", grb_falcon_model$exptr)
    #grb_falcon_model$addconstr(paste(fwd, '*', I, sep = ''), sense="<=", rhs= fwd_ub, name = paste(a_rxn, 'fwd', sep = '_')) # bound on fwd conversion
    #grb_falcon_model$addconstr(paste(rev, '*(1 - ', I, ')',  sep = ''), sense=">=", rhs= rev_lb, name = paste(a_rxn, 'rev', sep = '_')) # bound on rev conversion
    #.Call("GRB_updatemodel", grb_falcon_model$exptr)
  }

  # for each reaction w fwd and rev components:
  #   add constraint so that only one can run at a time
  #   binary vars I_fwd, I_rev <- {0, 1} and I_fwd + I_rev = 1

  .Call("GRB_updatemodel", grb_falcon_model$exptr)
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
  return(get_list_of_sets(return_couples(flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled)))
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
    print(paste('suppression index: ', i))
    if (!(vars[i] %in% unblocked_rxns)){
      print(paste(vars[i], ' blocked'))
      next
    }

    #prev_ub <- model$getattr("UB")[vars[i]]
    #prev_lb <- model$getattr("LB")[vars[i]]

    model <- model_og$copy()

    # block i
    model$setattr("UB", setNames(0, vars[i]))
    model$setattr("LB", setNames(0, vars[i]))

    set_lists[i] <- list(get_list_of_sets(return_couples(flux_coupling_raptor(model,reaction_indexes = reaction_indexes)$coupled)))

    # unfix i
    #model$setattr("UB", prev_ub)
    #model$setattr("LB", prev_lb)

  }

  return(set_lists)
}

GRB_generate_set_lists <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
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

  set_lists <- c()
  print(length(suppr_vector))
  print(paste("# of suppressions:", length(which(suppr_vector)), sep = " "))
  for (i in which(suppr_vector)){
    print(paste('suppression index:', i))
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

    coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_known_r0_sets, known_set_mtx = r0_coupling_mtx)$coupled
    set_lists[i] <- list(get_list_of_sets(return_couples(coupling_mtx)))
    # unfix i
    #model$setattr("UB", prev_ub)
    #model$setattr("LB", prev_lb)
  }

  #output <- list(r0_mtx <- r0_coupling_mtx, r1_coupling_array = coupling_array)

  #return(coupling_array)
  return(set_lists)
}

GRB_generate_set_lists_array <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
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

  coupling_array <- array(data = FALSE, dim = c(n,n,n), dimnames = list(vars, vars, paste('del', vars, sep = "_")))
  print(paste("# of suppressions:", length(which(suppr_vector)), sep = " "))
  for (i in which(suppr_vector)){
    print(paste('suppression index:', i))
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

    coupling_array[,,i] <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_known_r0_sets, known_set_mtx = r0_coupling_mtx)$coupled

    # unfix i
    #model$setattr("UB", prev_ub)
    #model$setattr("LB", prev_lb)

  }

  output <- list(r0_mtx <- r0_coupling_mtx, r1_coupling_array = coupling_array)

  #return(coupling_array)
  return(coupling_array)
}

coupling_matrix_from_array <- function(coupling_array){
  coupling_matrix <- matrix(data = 0, nrow = dim(coupling_array)[1], ncol = dim(coupling_array)[2],
                            dimnames = dimnames(coupling_array)[1:2])

  for (i in 1:nrow(coupling_matrix)){
    for (j in i:ncol(coupling_matrix)){
      #print(paste(i,j))
      coupling_matrix[i,j] <- length(which(coupling_array[i,j,]))
      #print(length(which(coupling_array[i,j,])))
    }
  }

  #coupling_matrix <- apply(coupling_array, 3, sum)
  return(coupling_matrix)
}

GRB_flux_coupling_raptor_wrapper <- function(i, vars, model_og, reaction_indexes = 1:length(vars), compare_mtx = FALSE, init_coupling_mtx = c(), file_output = NULL){
  print(paste('suppression index:', i))

  #if (!(init_coupling_mtx[i,i])){
  #  print(paste(vars[i], ' blocked'))
  #  next
  #}

  #prev_ub <- model$getattr("UB")[vars[i]]
  #prev_lb <- model$getattr("LB")[vars[i]]

  model <- model_og$copy()

  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))

  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_mtx, known_set_mtx = init_coupling_mtx, stored_obs = 4000)$coupled
  # coupling_list <- Matrix(data = FALSE, nrow = 1, ncol = length(coupling_mtx))
  # coupling_list[which(coupling_mtx)] <- TRUE

  #output <- list(which(coupling_mtx)) #list(get_list_of_sets(return_couples(coupling_mtx)))
  # return(coupling_list)
  coupling_idxs <- which(coupling_mtx)
  if (!is.null(file_output)){
    write(paste(c(i,coupling_idxs), collapse = ','),file = file_output,append=TRUE)
  }
  return(coupling_idxs)
}

GRB_generate_set_lists_cluster <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
                                         compare_known_init_sets = FALSE, optimize_suppr = FALSE, optimize_rxns = FALSE, cores = 1, avoid_idxs = c(), file_output = NULL){

  n <- model_og$get_sizes()$NumVars
  vars <- model_og$get_names()$VarName

  if (suppression_idxs[1] == -1){
    if (length(reaction_indexes) > 0){
      suppression_idxs = reaction_indexes
    }
    else {
      suppression_idxs = 1:n
    }
  }

  # dim: rxns_row, rxns_col, deletions
  model <- model_og$copy()

  #init_coupling_mtx <- c()
  if (compare_known_init_sets){
    init_coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled
    init_coupling_mtx <- fill_coupling_matrix(init_coupling_mtx)
  }
  suppr_vector <- Matrix(data = FALSE, nrow = 1, ncol = n, sparse = TRUE)
  suppr_vector[suppression_idxs] <- TRUE
  if (compare_known_init_sets & optimize_suppr){
    i <- 1
    while (i <= n){
      if (suppr_vector[i]){ # if tagged to be suppressed
        set_idx <- which(init_coupling_mtx[,i])[1] # which is first reaction (row) i is coupled to
        if (!is.na(set_idx)){
          rxn_idxs <- which(init_coupling_mtx[set_idx,]) # other reactions in set
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

    if (optimize_rxns){
      reaction_indexes <- which(suppr_vector)
    }
  }

  if (length(avoid_idxs) > 0){
    suppr_vector[avoid_idxs] <- FALSE
  }

  print(paste("# of suppressions:", length(which(suppr_vector)), sep = " "))
  coupling <- mclapply(which(suppr_vector), function(x) GRB_flux_coupling_raptor_wrapper(x, vars, model_og, reaction_indexes = reaction_indexes, compare_mtx = compare_known_init_sets, init_coupling_mtx = init_coupling_mtx, file_output = file_output), mc.cores = cores)

  #coupling_idxs <- c()
  #for (i in which(suppr_vector)){
  #  intermediate_coupling_idxs <- GRB_flux_coupling_raptor_wrapper(i, vars, model_og, reaction_indexes = reaction_indexes, compare_mtx = compare_known_init_sets, init_coupling_mtx = init_coupling_mtx)
  #  #write(paste(),file="myfile",append=TRUE)
  #  coupling_idxs <- c(coupling_idxs, intermediate_coupling_idxs)
  #  coupling_idxs <- unique(coupling_idxs)
  #}

  return(coupling)
}

coupling_matrix_from_coupling_vector_list <- function(coupling_list, n_react, vars, init_sets = NULL){

  #len <- length(coupling_list[[1]])
  #len <- n_react*n_react
  #coupling_vector <- Matrix(data = FALSE, nrow = 1, ncol = len)

  #for (i in 1:length(coupling_list)){
  #  coupling_vector[coupling_list[[i]]] <- TRUE
  #}

  #matrix_dim_size <- sqrt(len)
  coupling_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(coupling_matrix) <- vars
  colnames(coupling_matrix) <- vars
  if (!is.null(init_sets)){
    coupling_matrix <- fill_coupling_matrix_from_sets(coupling_matrix, init_sets)
  }
  
  for (i in 1:length(coupling_list)){
    if (is.null(coupling_list[[i]])){next}
    coupling_matrix[coupling_list[[i]]] <- TRUE
  }
  
  coupling_matrix <- fill_coupling_matrix(coupling_matrix)
  #coupling_matrix[which(coupling_vector)] <- TRUE

  return(coupling_matrix)
}

full_ish_coupling_matrix_from_coupling_vector_list <- function(coupling_list, n_react, vars, init_sets){
  
  #matrix_dim_size <- sqrt(len)
  coupling_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(coupling_matrix) <- vars
  colnames(coupling_matrix) <- vars
  for (i in 1:length(coupling_list)){
    if (is.null(coupling_list[[i]])){next}
    print(i)
    intermediate_mtx <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
    rownames(intermediate_mtx) <- vars
    colnames(intermediate_mtx) <- vars
    intermediate_mtx[coupling_list[[i]]] <- TRUE
    intermediate_mtx <- fill_coupling_matrix_from_sets(intermediate_mtx, init_sets)
    coupling_matrix[which(intermediate_mtx)] <- TRUE
  }
  #coupling_matrix[which(coupling_vector)] <- TRUE
  
  return(coupling_matrix)
}

identify_intermediate_uncoupled <- function(full_coupling_mtx, full_ish_coupling_mtx, n_react){
  uncoupled_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(uncoupled_matrix) <- vars
  colnames(uncoupled_matrix) <- vars
  
  uncoupled <- which(full_coupling_mtx & !full_ish_coupling_mtx)
  uncoupled_matrix[uncoupled] <- TRUE
  
  return(uncoupled_matrix)
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

GRB_maximize <- function(model_og, obj, suppress = c(), max = TRUE){ # suppress is characters
  model <- model_og$copy()

  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  # clear obj
  model$setattr("Obj", setNames(rep(0.0, times = n), vars))

  # set suppressions
  if (length(suppress) > 0){
    #print(suppress)
    suppr_idxs <- which(vars %in% suppress)
    model$setattr("UB", setNames(rep(0.0, times = length(suppr_idxs)), vars[suppr_idxs]))
    model$setattr("LB", setNames(rep(0.0, times = length(suppr_idxs)), vars[suppr_idxs]))
  }

  # set obj
  model$setattr("Obj", setNames(1.0, vars[obj]))
  model$set_model_sense(maximize=TRUE)
  if (!max){
    model$set_model_sense(minimize=TRUE)
  }
  model$optimize()
  sol <- model$get_solution()
  obj_max <- sol$ObjVal
  return(obj_max)
}
