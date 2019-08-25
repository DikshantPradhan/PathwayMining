library(gurobi)
library(sybil)
library(dplyr)


insert_into_list <- function(list, insert, idx){
  segment_1 <- list[1:idx]
  segment_2 <- list[(idx+1):length(list)]
  
  new_list <- c(segment_1, insert, segment_2)
  return(new_list)
}

# compare_lp_files <- 

rewrite_lp <- function(filename = 'temp.lp', model){
  lines <- readLines(filename)
  
  for (i in rev(1:length(model@react_id))){
    rxn <- model@react_id[i]
    rxn <- gsub(' ', replacement = '_', rxn)
    lp_rxn_id <- paste('x_', i, sep = '')
    # print(paste(lp_rxn_id, rxn, sep = '; '))
    lines <- gsub(pattern = lp_rxn_id, replacement = rxn, lines)
  }
  
  for (i in rev(1:length(model@met_id))){
    met <- model@met_id[i]
    met <- gsub(' ', replacement = '_', met)
    lp_met_id <- paste('r_', i, sep = '')
    # print(paste(lp_met_id, met, sep = '; '))
    lines <- gsub(lp_met_id, met, lines)
  }
  
  return(lines)
}

convert_sybil_to_lp <- function(sybil, output = 'temp.lp'){
  # optimizeProb(sybil, poCmd = list(c("writeProb", "LP_PROB", output_file, "'lp'")))
  # new_lines <- rewrite_lp(output, sybil)
  # write(new_lines, file = output)
  prob <- sysBiolAlg(sybil, useNames=TRUE)
  writeProb(problem(prob), fname = output)
}

convert_sybil_to_gurobi <- function(sybil, output = 'temp.lp'){
  convert_sybil_to_lp(sybil, output = output)
  # prob <- sysBiolAlg(sybil, useNames=TRUE)
  # writeProb(problem(prob), fname = output)
  model <- gurobi_read(output)
  return(model)
}

GRB_generate_falcon_model <- function(sybil_model, falcon_model = FALSE, r0_gene_set = c(), r0_rxn_set_list = c()){
  
  sybil_falcon_model <- sybil_model
  
  if (!falcon_model){ # if user has not passed in a falcon model already
    sybil_falcon_model <- generate_falcon_model(sybil_model, r0_gene_set, r0_rxn_set_list)
  }
  
  # grb_falcon_model <- as_GRBmodel(sybil_falcon_model)
  # grb_falcon_model$show_output(FALSE)
  
  ## ADD NECESSARY CONSTRAINTS TO MODEL
  vars <- sybil_falcon_model@react_id #grb_falcon_model$get_names()$VarName
  
  split_fwd_rxns <- vars[grep('fwd', vars)]
  split_rev_rxns <- vars[grep('rev', vars)]
  
  split_rxns <- sapply(split_fwd_rxns, function(x) strsplit(x, ' ')[[1]][1])
  split_rxns <- unique(split_rxns)
  
  if (length(split_fwd_rxns) != length(split_rev_rxns)){
    print('fed rev matchup error')
    return()
  }
  
  Binaries <- c('', 'Binaries')
  
  for (rxn in split_rxns){
    I <- paste('I', rxn, sep = '_')
    Binaries <- c(Binaries, I)
  }
  
  # vars <- grb_falcon_model$get_names()$VarName
  # n <- grb_falcon_model$get_sizes()$NumVars
  
  new_constr <- c()
  
  # add bounds on split conversion reactions (gene -> activity_[rxn])
  for (i in 1:length(split_fwd_rxns)){
    fwd <- split_fwd_rxns[i]
    fwd <- gsub(pattern = ' ', replacement = '_', fwd)
    rev <- split_rev_rxns[i]
    rev <- gsub(pattern = ' ', replacement = '_', rev)
    
    # fwd_ub <- grb_falcon_model$getattr("UB")[[fwd]]
    # fwd_lb <- grb_falcon_model$getattr("LB")[[fwd]]
    # rev_ub <- grb_falcon_model$getattr("UB")[[rev]]
    # rev_lb <- grb_falcon_model$getattr("LB")[[rev]]
    
    a_rxn <- strsplit(fwd, ' ')[[1]][1]
    I <- paste('I', a_rxn, sep = '_')
    # i_idx <- which(vars == I)
    # fwd_idx <- which(vars == fwd)
    # rev_idx <- which(vars == rev)
    
    #print(paste(fwd, fwd_ub, fwd_lb, rev, rev_ub, rev_lb, I))
    # fwd_vec <- c(fwd_ub, -1)
    # rev_vec <- c(-1*rev_ub, -1)
    # fwd_idxs <- c(i_idx, which(vars == fwd))
    # rev_idxs <- c(i_idx, which(vars == rev))
    
    fwd_name <- paste(fwd, 'I', sep = '_')
    # fwd_name <- gsub(' ', '_', fwd_name)
    rev_name <- paste(rev, 'I', sep = '_')
    # rev_name <- gsub(' ', '_', rev_name)
    
    fwd_line <- paste(' ', fwd_name, ": - ", fwd, ' + 1000 ', I, ' >= 0', sep = '')
    rev_line <- paste(' ',rev_name, ": - ", rev, ' - 1000 ', I, ' >= -1000', sep = '')
    
    new_constr <- c(new_constr, fwd_line, rev_line)
  }
  
  # write and read model
  convert_sybil_to_lp(sybil_falcon_model)
  test_model <- gurobi_read('temp.lp')
  print('tested')
  lp_lines <- rewrite_lp(model = sybil_falcon_model)
  
  constr_idx <- grep('Bounds', lp_lines)
  lp_lines <- insert_into_list(lp_lines, new_constr, constr_idx-2)
  end_idx <- grep('End', lp_lines)
  lp_lines <- lp_lines <- insert_into_list(lp_lines, Binaries, (end_idx-2))
  write(lp_lines, file = 'temp.lp')
  model <- gurobi_read('temp.lp')
  
  return(model)
}

solve_model <- function(model, obj_idx, sense = 'max'){
  model$obj <- rep(0, length(model$obj))
  model$obj[obj_idx] <- 1
  model$modelsense <- sense
  sol <- gurobi(model, list(OutputFlag = 0))
  return(sol)
}

flux_coupling_raptor <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01,
                                 bnd_tol = 0.1, stored_obs = 4000, cor_iter = 5, cor_check = TRUE,
                                 reaction_indexes = c(), compare_mtx = FALSE, 
                                 known_set_mtx = Matrix(data = FALSE, nrow = 1, ncol = 1, sparse = TRUE),
                                 directional = FALSE) {
  
  # min_fva_cor is minimum correlation between fluxes
  # bnd_tol is allowed error in comparing max & min flux
  # fix_frac is const value used in fixing flux at non-zero value
  # fix_tol_frac is error allowed in determining whether flux is fixed
  # stored_obs is # not flux values to be stored
  # cor_iter is number of iterations after which correlation is considered in checking coupling
  
  if (!cor_check){
    stored_obs = 1
  }
  
  vars <- model$varnames
  n <- length(vars)
  
  if (is.null(known_set_mtx)){
    compare_mtx <- FALSE
  }
  else {
    if ((nrow(known_set_mtx) < n) | (ncol(known_set_mtx) < n)){compare_mtx <- FALSE}
  }
  
  # if empty set, then assume all reactions are to be inspected
  if (length(reaction_indexes) == 0){
    reaction_indexes <- c(1:n)
  }
  
  prev_obj <- model$obj
  model$obj <- rep(0,n) # clear the objective
  prev_sense <- model$modelsense
  
  original_ub <- model$ub
  original_lb <- model$lb
  
  global_max <- rep(0, n)
  global_min <- rep(0, n)
  
  coupled <- Matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars), sparse = TRUE)
  
  blocked <- rep(FALSE, n)
  active <- rep(TRUE, n)
  
  not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }
  
  rxn_fix <- function(max_, min_){
    
    avg <- min_ + fix_frac*(max_ - min_)
    ct = 0
    while (near(avg, 0) & ct < 5){
      print(c(avg, max_, min_))
      avg <- avg  + fix_frac*(max_ - min_)
      ct = ct + 1
    }
    
    return(avg)
  }
  
  correlation_check <- function(flux, i, j){ # return true if correlation is high, or no correlation --> continue comparing flux
    # if false, skip in depth comparison
    
    n_entries <- length(which(!near(flux[,i], 0, tol = bnd_tol) | !near(flux[,j], 0, tol = bnd_tol)))
    # n_entries <- length(which(flux[,i] != 0 | flux[,j] != 0))
    if (n_entries > cor_iter){
      C <- cor(flux[,i], flux[,j])
      
      if ((is.na(C)) | (abs(C) < min_fva_cor)){
        return(FALSE)
      }
    }
    return(TRUE)
  }
  
  flux <- matrix(c(0), nrow = stored_obs, ncol = n)
  lp_calls <- 0
  
  update_flux <- function(flux_, idx, sol){
    if (stored_obs > 0){
      flux_[idx,] <- sol
    }
    return(flux_)
  }
  
  # if j is coupled to i, couple the rxns known to be coupled to j to i as well
  known_set_coupling <- function(i, j, coupled, active){
    #sets <- which(known_set_mtx[,j])
    set <- which(known_set_mtx[j,])
    
    #for (k in set){
    #  coupled[i, k] <- TRUE
    #  coupled[k, k] <- TRUE
    #}
    
    coupled[i, set] <- TRUE
    coupled[set, set] <- TRUE
    active[set] <- FALSE
    #for (set in sets){
    #known_coupled_rxns <- which(known_set_mtx[set,])
    
    #coupled[i, known_coupled_rxns] <- TRUE
    #active[known_coupled_rxns] <- FALSE
    #}
    
    list(coupled = coupled, active = active)
  }
  
  check_directional_coupling <- function(i, j, fixed = TRUE){
    
    info <- list(lp_calls = 0,
                 coupled = FALSE)
    
    prev_ub_i <- model$ub[i]
    prev_lb_i <- model$lb[i]
    
    if (!fixed){
      fixed_val <- 0.01
      if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
        fixed_val <- rxn_fix(global_max[i], global_min[i])
      }
      else {
        if (!near(model$ub[i], 0)){ #model$getattr("UB")[vars[i]] > tol
          
          model$obj[i] <- 1
          model$modelsense <- 'max'
          sol <- gurobi(model, list(OutputFlag = 0))
          info$lp_calls <- info$lp_calls + 1
          global_max <- pmax(global_max, sol$x)
          global_min <- pmin(global_min, sol$x)
          #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
          # flux[lp_calls%%stored_obs,] <- sol$x
          
          if (!near(global_max[i], 0) | !near(global_min[i], 0)){
            fixed_val <- rxn_fix(global_max[i], global_min[i])
          }
          
          model$obj <- rep(0, n)
        }
        if (!near(model$lb[i], 0)){ #model$getattr("LB")[vars[i]] < (-1*tol_)
          
          model$obj[i] <- 1
          model$modelsense <- 'min'
          sol <- gurobi(model, list(OutputFlag = 0))
          info$lp_calls <- info$lp_calls + 1
          #print(is.null(flux))
          global_max <- pmax(global_max, sol$x)
          global_min <- pmin(global_min, sol$x)
          #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
          # flux[lp_calls%%stored_obs,] <- sol$x
          
          if (!near(global_max[i], 0) | !near(global_min[i], 0)){
            fixed_val <- rxn_fix(global_max[i], global_min[i])
          }
          
          model$obj <- rep(0, n)
        }
        
        model$ub[i] <- fixed_val
        model$lb[i] <- fixed_val
      }
    }
    
    prev_ub_j <- model$ub[j]
    prev_lb_j <- model$lb[j]
    
    # zero out target bounds
    model$ub[j] <- 0
    model$lb[j] <- 0
    
    # check feasibility
    sol <- gurobi(model, list(OutputFlag = 0))
    
    info$lp_calls <- info$lp_calls + 1
    
    feasible <- !(length(sol$x) == 0)
    info$coupled <- feasible
    
    model$ub[j] <- prev_ub_j
    model$lb[j] <- prev_lb_j
    
    model$ub[i] <- prev_ub_i
    model$lb[i] <- prev_lb_i
    
    
    return(info)
  }
  
  
  for (idx in 1:(length(reaction_indexes))) { # (i in 1:(n-1))
    # iterate over passed in idxs instead (idx in 1:length(reaction_indexes)); i <-  reaction_indexes[idx]
    i <-  reaction_indexes[idx]
    
    if (!active[i] | blocked[i]) next
    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)
    prev_ub <- model$ub[i]
    prev_lb <- model$lb[i]
    
    if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
      fixed_val <- rxn_fix(global_max[i], global_min[i])
    }
    else {
      if (!near(model$ub[i], 0)){ #model$getattr("UB")[vars[i]] > tol
        
        model$obj[i] <- 1
        model$modelsense <- 'max'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        #print(is.null(flux))
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        model$obj <- rep(0, n)
      }
      if (!near(model$lb[i], 0)){ #model$getattr("LB")[vars[i]] < (-1*tol_)
        
        model$obj[i] <- 1
        model$modelsense <- 'min'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        #print(is.null(flux))
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        model$obj <- rep(0, n)
      }
      if (near(global_max[i], 0) & near(global_min[i], 0)){ #(abs(global_max[i]) < tol) & (abs(global_min[i]) < tol)
        #print(paste('blocked:', i))
        blocked[i] <- TRUE
        active[i] <- FALSE
        next
      }
    }
    
    # set new bounds for selected rxn (temporarily)
    model$ub[i] <- fixed_val #setattr("UB", setNames(fixed_val + 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    model$lb[i] <- fixed_val #setattr("LB", setNames(fixed_val - 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    
    # couple reaction to itself if not blocked
    if (!blocked[i]){
      coupled[i,i] <- TRUE
      active[i] <- FALSE
    }
    
    if (idx == length(reaction_indexes)){break}
    
    for (idx2 in (idx+1):length(reaction_indexes)) { # (j in (i+1):n)
      # also keep this in passed in idxs (idx2 in (idx+1):length(reaction_indexes)); j <-  reaction_indexes[idx2]
      j <-  reaction_indexes[idx2]
      # check for fixed or blocked
      if (!active[j] | blocked[j]){next}
      if (not_fixed(sub_max[j], sub_min[j])){next}
      #print(paste(sub_max[j], sub_min[j]))
      
      # check for uncoupled via correlation
      if (cor_check){
        if (!correlation_check(flux, i, j)){next}
        # C <- cor(flux[,i], flux[,j])
        # if (!is.na(C) & (abs(C) < min_fva_cor)){next}
      }
      
      skip <- FALSE
      
      max <- 0
      min <- 0
      
      # model$setattr("Obj", setNames(1.0, vars[j]))
      model$obj[j] <- 1
      if (!skip) {
        model$modelsense <- 'max'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        #print(is.null(flux))
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        sub_max <- pmax(sub_max, sol$x)
        sub_min <- pmin(sub_min, sol$x)
        
        #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        max <- sol$x[j]
        
        skip <- not_fixed(sub_max[j], sub_min[j])
      }
      
      if (!skip) {
        model$modelsense <- 'min'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        #print(is.null(flux))
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        sub_max <- pmax(sub_max, sol$x)
        sub_min <- pmin(sub_min, sol$x)
        
        #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        max <- sol$x[j]
        
        skip <- not_fixed(sub_max[j], sub_min[j])
      }
      
      if (near(max, 0) & near(min, 0)){skip = TRUE}
      
      if (!skip) { # finally label as coupled
        coupled[i,j] <- TRUE
        coupled[j,j] <- TRUE # make sure the reaction doesn't look blocked
        active[j] <- FALSE
        
        if (compare_mtx){
          output <- known_set_coupling(i, j, coupled, active)
          coupled <- output$coupled
          active <- output$active
        }
      }
      
      if (skip & directional){
        dir_i_j <- check_directional_coupling(i, j, fixed = TRUE)
        lp_calls <- lp_calls + dir_i_j$lp_calls
        
        if (dir_i_j$coupled){
          coupled[i,j] <- TRUE
          next
        }
        
        dir_j_i <- check_directional_coupling(j, i, fixed = FALSE)
        lp_calls <- lp_calls + dir_j_i$lp_calls
        
        if (dir_j_i$coupled){
          coupled[j,i] <- TRUE
        }
        
        if (dir_i_j$coupled | dir_j_i$coupled){
          coupled[i,j] <- TRUE
          coupled[j,i] <- TRUE
          active[j] <- FALSE
        }
      }
      
      # if (skip){
      #   
      # }
      
      # print(lp_calls)
      
      # model$setattr("Obj", setNames(0.0, vars[j]))
      model$obj <- rep(0, n)
    }
    
    
    # unfix i
    model$ub <- original_ub #setattr("UB", prev_ub)
    model$lb <- original_lb #setattr("LB", prev_lb)
  }
  
  #i <- reaction_indexes[length(reaction_indexes)]
  #if (!blocked[i] & active[i]){
  #  coupled[i,i] <- TRUE
  #  active[i] <- FALSE
  #}
  
  model$obj <- prev_obj#setattr("Obj", prev_obj)
  model$modelsense <- prev_sense #setattr("ModelSense", prev_sense)
  
  print(lp_calls)
  
  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}


partial_flux_coupling_raptor <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01,
                                 bnd_tol = 0.1, stored_obs = 4000, cor_iter = 5, cor_check = TRUE,
                                 reaction_indexes = c(), compare_mtx = FALSE, 
                                 known_set_mtx = Matrix(data = FALSE, nrow = 1, ncol = 1, sparse = TRUE)) {
  
  # min_fva_cor is minimum correlation between fluxes
  # bnd_tol is allowed error in comparing max & min flux
  # fix_frac is const value used in fixing flux at non-zero value
  # fix_tol_frac is error allowed in determining whether flux is fixed
  # stored_obs is # not flux values to be stored
  # cor_iter is number of iterations after which correlation is considered in checking coupling
  
  if (!cor_check){
    stored_obs = 1
  }
  
  vars <- model$varnames
  n <- length(vars)
  
  if (is.null(known_set_mtx)){
    compare_mtx <- FALSE
  }
  else {
    if ((nrow(known_set_mtx) < n) | (ncol(known_set_mtx) < n)){compare_mtx <- FALSE}
  }
  
  # if empty set, then assume all reactions are to be inspected
  if (length(reaction_indexes) == 0){
    reaction_indexes <- c(1:n)
  }
  
  prev_obj <- model$obj
  model$obj <- rep(0,n) # clear the objective
  prev_sense <- model$modelsense
  
  original_ub <- model$ub
  original_lb <- model$lb
  
  global_max <- rep(0, n)
  global_min <- rep(0, n)
  
  coupled <- Matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars), sparse = TRUE)
  
  blocked <- rep(FALSE, n)
  active <- rep(TRUE, n)
  
  not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }
  
  rxn_fix <- function(max_, min_){
    
    avg <- min_ + fix_frac*(max_ - min_)
    ct = 0
    while (near(avg, 0) & ct < 5){
      print(c(avg, max_, min_))
      avg <- avg  + fix_frac*(max_ - min_)
      ct = ct + 1
    }
    
    return(avg)
  }
  
  correlation_check <- function(flux, i, j){ # return true if correlation is high, or no correlation --> continue comparing flux
    # if false, skip in depth comparison
    
    n_entries <- length(which(!near(flux[,i], 0, tol = bnd_tol) | !near(flux[,j], 0, tol = bnd_tol)))
    # n_entries <- length(which(flux[,i] != 0 | flux[,j] != 0))
    if (n_entries > cor_iter){
      C <- cor(flux[,i], flux[,j])
      
      if ((is.na(C)) | (abs(C) < min_fva_cor)){
        return(FALSE)
      }
    }
    return(TRUE)
  }
  
  flux <- matrix(c(0), nrow = stored_obs, ncol = n)
  lp_calls <- 0
  
  update_flux <- function(flux_, idx, sol){
    if (stored_obs > 0){
      flux_[idx,] <- sol
    }
    return(flux_)
  }
  
  # if j is coupled to i, couple the rxns known to be coupled to j to i as well
  known_set_coupling <- function(i, j, coupled, active){
    #sets <- which(known_set_mtx[,j])
    set <- which(known_set_mtx[j,])
    
    #for (k in set){
    #  coupled[i, k] <- TRUE
    #  coupled[k, k] <- TRUE
    #}
    
    coupled[i, set] <- TRUE
    coupled[set, set] <- TRUE
    active[set] <- FALSE
    #for (set in sets){
    #known_coupled_rxns <- which(known_set_mtx[set,])
    
    #coupled[i, known_coupled_rxns] <- TRUE
    #active[known_coupled_rxns] <- FALSE
    #}
    
    list(coupled = coupled, active = active)
  }
  
  fix_flux_val <- function(i){ # idx of reaction to fix
    
    fixed_val <- 0
    
    if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
      fixed_val <- rxn_fix(global_max[i], global_min[i])
    }
    else {
      if (!near(model$ub[i], 0)){
        
        model$obj[i] <- 1
        model$modelsense <- 'max'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        #print(is.null(flux))
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        model$obj <- rep(0, n)
      }
      if (!near(model$lb[i], 0)){
        
        model$obj[i] <- 1
        model$modelsense <- 'min'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        #print(is.null(flux))
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        model$obj <- rep(0, n)
      }
      if (near(global_max[i], 0) & near(global_min[i], 0)){ #(abs(global_max[i]) < tol) & (abs(global_min[i]) < tol)
        #print(paste('blocked:', i))
        blocked[i] <- TRUE
        active[i] <- FALSE
        #next
        return(0)
      }
    }
    
    # set new bounds for selected rxn (temporarily)
    # model$ub[i] <- fixed_val #setattr("UB", setNames(fixed_val + 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    # model$lb[i] <- fixed_val #setattr("LB", setNames(fixed_val - 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    return(fixed_val)
  }
  
  check_directional_coupling <- function(i, j){
    coupled <- FALSE
    ret <- list(coupled = coupled,
                flux = flux,
                sub_max = sub_max,
                sub_min = sub_min,
                lp_calls = lp_calls)
    
    # check stored values
    idxs <- which(flux[,i] == 0)
    flag <- any(flux[idxs,j] != 0)
    if (flag){
      return(ret)
    }
    
    
    # fix i to zero
    model$ub[i] <- 0
    model$lb[i] <- 0
    
    # observe values for j
    skip <- FALSE
    
    max <- Inf #0
    min <- -Inf #0
    
    # model$setattr("Obj", setNames(1.0, vars[j]))
    model$obj[j] <- 1
    if (!skip) {
      model$modelsense <- 'max'
      sol <- gurobi(model, list(OutputFlag = 0))
      lp_calls <- lp_calls + 1
      if (length(sol$x) == 0){return(ret)}
      #print(is.null(flux))
      # global_max <- pmax(global_max, sol$x)
      # global_min <- pmin(global_min, sol$x)
      sub_max <- pmax(sub_max, sol$x)
      sub_min <- pmin(sub_min, sol$x)
      
      #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
      flux[lp_calls%%stored_obs,] <- sol$x
      
      max <- sol$x[j]
      
      #skip <- not_fixed(sub_max[j], sub_min[j])
    }
    
    if (!skip) {
      model$modelsense <- 'min'
      sol <- gurobi(model, list(OutputFlag = 0))
      lp_calls <- lp_calls + 1
      if (length(sol$x) == 0){return(ret)}
      #print(is.null(flux))
      # global_max <- pmax(global_max, sol$x)
      # global_min <- pmin(global_min, sol$x)
      sub_max <- pmax(sub_max, sol$x)
      sub_min <- pmin(sub_min, sol$x)
      
      #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
      flux[lp_calls%%stored_obs,] <- sol$x
      
      min <- sol$x[j]
      
      #skip <- not_fixed(sub_max[j], sub_min[j])
    }
    # return flux, sub_max, sub_min, lp
    
    #print(paste(min, max))
    if (length(min) < 1 | length(max) < 1){
      ret$coupled <- FALSE
      return(ret)
      # return(FALSE)
      } 
    if (!near(max, 0) | !near(min, 0)){
      ret$coupled <- FALSE
      return(ret)
      # return(FALSE)
      }#{skip = TRUE}
    
    ## fix i to non-zero value
    
    fixed_val <- fix_flux_val(i)
    if (fixed_val == 0){
      ret$coupled <- FALSE
      return(ret)
      # return(FALSE)
      }
    
    model$ub[i] <- fixed_val
    model$lb[i] <- fixed_val
    
    ## Check j
    
    model$modelsense <- 'max'
    sol <- gurobi(model, list(OutputFlag = 0))
    lp_calls <- lp_calls + 1
    if (length(sol$x) == 0){return(ret)}
    #print(is.null(flux))
    # global_max <- pmax(global_max, sol$x)
    # global_min <- pmin(global_min, sol$x)
    sub_max <- pmax(sub_max, sol$x)
    sub_min <- pmin(sub_min, sol$x)
    
    #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
    flux[lp_calls%%stored_obs,] <- sol$x
    
    max <- sol$x[j]
    if (!near(max, 0)){
      ret$coupled <- TRUE
      ret$flux <- flux
      ret$sub_max <- sub_max
      ret$sub_min <- sub_min
      ret$lp_calls <- lp_calls
      return(ret)
      # return(TRUE)
      }
    
    model$modelsense <- 'min'
    sol <- gurobi(model, list(OutputFlag = 0))
    lp_calls <- lp_calls + 1
    #print(is.null(flux))
    # global_max <- pmax(global_max, sol$x)
    # global_min <- pmin(global_min, sol$x)
    sub_max <- pmax(sub_max, sol$x)
    sub_min <- pmin(sub_min, sol$x)
    
    min <- sol$x[j]
    if (!near(min, 0)){
      ret$coupled <- TRUE
      ret$flux <- flux
      ret$sub_max <- sub_max
      ret$sub_min <- sub_min
      ret$lp_calls <- lp_calls
      return(ret)
      # return(TRUE)
      }
    
    ret$coupled <- FALSE
    ret$flux <- flux
    ret$sub_max <- sub_max
    ret$sub_min <- sub_min
    ret$lp_calls <- lp_calls
    return(ret)
    
    # return(FALSE)
  }
  
  for (idx in 1:(length(reaction_indexes))) { # (i in 1:(n-1))
    # iterate over passed in idxs instead (idx in 1:length(reaction_indexes)); i <-  reaction_indexes[idx]
    i <-  reaction_indexes[idx]
    
    if (!active[i] | blocked[i]) next
    
    # set i to 0, check for variable or zero flux
    sub_max <- rep(0, n)
    sub_min <- rep(0, n)
    prev_ub <- model$ub[i]
    prev_lb <- model$lb[i]
    
    if (!blocked[i]){
      coupled[i,i] <- TRUE
      active[i] <- FALSE
    }
    
    if (idx == length(reaction_indexes)){break}
    
    for (idx2 in (idx+1):length(reaction_indexes)) { # (j in (i+1):n)
      # also keep this in passed in idxs (idx2 in (idx+1):length(reaction_indexes)); j <-  reaction_indexes[idx2]
      j <-  reaction_indexes[idx2]
      # check for fixed or blocked
      if (!active[j] | blocked[j]){next}
      
      # CHECK I -> J
      if ((sub_max[j] == 0) & (sub_min[j] == 0)){
        ret <- check_directional_coupling(i,j)
        
        dir_i_j <- ret$coupled
        flux <- ret$flux
        sub_max <- ret$sub_max
        sub_min <- ret$sub_min
        lp_calls <- ret$lp_calls
        
      }
      # CHECK J -> I
      
      ret <- check_directional_coupling(j,i)
      
      dir_j_i <- ret$coupled
      flux <- ret$flux
      sub_max <- ret$sub_max
      sub_min <- ret$sub_min
      lp_calls <- ret$lp_calls
        
      # dir_i_j <- check_directional_coupling(i,j)
      # dir_j_i <- check_directional_coupling(j,i)
      
      #coupled[i,j] <- dir_i_j
      #coupled[j,j] <- dir_j_i
      if (dir_i_j & dir_j_i){coupled[i,j] <- TRUE; coupled[j,i] <- TRUE}
      
      # model$setattr("Obj", setNames(0.0, vars[j]))
      model$obj <- rep(0, n)
    }
    
    
    # unfix i
    model$ub <- original_ub #setattr("UB", prev_ub)
    model$lb <- original_lb #setattr("LB", prev_lb)
  }
  
  
  model$obj <- prev_obj#setattr("Obj", prev_obj)
  model$modelsense <- prev_sense #setattr("ModelSense", prev_sense)
  
  print(lp_calls)
  
  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}


GRB_get_rxn_idx <- function(model, rxn){
  vars <- model$varnames
  return(return(which(vars == rxn)))
}

GRB_generate_pair_list <- function(model_og){
  model <- model_og
  return(return_couples(flux_coupling_raptor(model)$coupled))
}

GRB_generate_set_list <- function(model_og, reaction_indexes = c()){
  model <- model_og
  return(get_list_of_sets(return_couples(flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled)))
}

GRB_flux_coupling_raptor_wrapper <- function(i, vars, model_og, reaction_indexes = 1:length(vars), compare_mtx = FALSE, init_coupling_mtx = c(), file_output = NULL){
  print(paste('suppression index:', i))
  
  model <- model_og
  
  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))
  
  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_mtx, known_set_mtx = init_coupling_mtx, stored_obs = 4000)$coupled
  
  coupling_idxs <- which(coupling_mtx)
  if (!is.null(file_output)){
    write(paste(c(i,coupling_idxs), collapse = ','),file = file_output,append=TRUE)
  }
  return(coupling_idxs)
}

GRB_generate_set_lists_cluster <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
                                           compare_known_init_sets = FALSE, optimize_suppr = FALSE, optimize_rxns = FALSE, cores = 1, avoid_idxs = c(), file_output = NULL){
  
  vars <- model_og$varnames
  n <- length(vars)
  
  if (suppression_idxs[1] == -1){
    if (length(reaction_indexes) > 0){
      suppression_idxs = reaction_indexes
    }
    else {
      suppression_idxs = 1:n
    }
  }
  
  # dim: rxns_row, rxns_col, deletions
  model <- model_og
  
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
  
  return(coupling)
}

GRB_maximize <- function(model_og, obj, suppress = c(), max = TRUE){ # suppress is characters
  model <- model_og
  
  vars <- model$varnames
  n <- length(vars)
  
  # clear obj
  model$obj <- rep(0.0, times = n) # setattr("Obj", setNames(rep(0.0, times = n), vars))
  
  # set suppressions
  if (length(suppress) > 0){
    suppr_idxs <- which(vars %in% suppress)
    model$ub[suppr_idxs] <- 0 #setattr("UB", setNames(rep(0.0, times = length(suppr_idxs)), vars[suppr_idxs]))
    model$lb[suppr_idxs] <- 0 #setattr("LB", setNames(rep(0.0, times = length(suppr_idxs)), vars[suppr_idxs]))
  }
  
  # set obj
  sense <- 'max'
  if (!max){
    sense <- 'min' # model$set_model_sense(minimize=TRUE)
  }
  sol <- solve_model(model = model, obj_idx = obj, sense = sense)
  obj_max <- sol$objval
  return(obj_max)
}

## extra functions

fill_coupling_matrix <- function(coupled){
  rows <- nrow(coupled)
  
  for (i in 1:nrow(coupled)){
    #identify set
    # if (!coupled[i,i]){next}
    set <- which(coupled[i,]) # true values in row
    if (length(set) < 1){next}
    set <- unique(c(i, set))
    coupled[set,set] <- TRUE
    
  }
  
  coupled[lower.tri(coupled)] <- FALSE
  return(coupled)
}

fill_coupling_matrix_from_sets <- function(mtx, sets){
  for (set in sets){
    mtx[set,set] <- TRUE
  }
  
  mtx <- fill_coupling_matrix(mtx)
  
  return(mtx)
}

set_vector <- function(coupled){
  active <- matrix(data = TRUE, nrow = 1, ncol = ncol(coupled))
  set_num <- matrix(data = 0, nrow = 1, ncol = ncol(coupled))
  
  set_iter <- 1
  for (i in nrow(coupled)){
    if (!active[i]){next} # skip if already in a set
    set <- which(coupled[i,])
    if (length(set) < 1){next} # skip if blocked reaction (will not be coupled to itself)
    active[set] <- FALSE
    set_num[set] <- set_iter # enter same set # into all indexes in the same set
    set_iter <- set_iter + 1
  }
  
  return(set_num)
}


# data("Ec_core")
#load('GitHub/PathwayMining_Package/data/pao_model.RData')
# model <- mutans_model
# optimizeProb(Ec_core, poCmd = list(c("writeProb", "LP_PROB", "'temp2.lp'", "'lp'")))
# test_model <- gurobi_read('temp2.lp')
# model <- convert_sybil_to_gurobi(pao_model)
# model <- gurobi_read('~/GitHub/PathwayMining/new_ecoli.lp')
#model <- convert_sybil_to_gurobi(Ec_core)
# falcon_model <- GRB_generate_falcon_model(model)

#r_mtx <- flux_coupling_raptor(model, fix_tol_frac = 0.00001)$coupled
#r_sets <- get_list_of_sets_from_mtx(r_mtx)
