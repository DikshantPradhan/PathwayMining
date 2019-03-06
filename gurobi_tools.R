library(gurobi)
library(sybil)
library(dplyr)

insert_into_list <- function(list, insert, idx){
  segment_1 <- list[1:idx]
  segment_2 <- list[(idx+1):length(list)]
  
  new_list <- c(segment_1, insert, segment_2)
  return(new_list)
}

rewrite_lp <- function(filename = 'temp.lp', model){
  lines <- readLines(filename)
  
  for (i in rev(1:length(model@react_id))){
    rxn <- model@react_id[i]
    lp_rxn_id <- paste('x_', i, sep = '')
    # print(paste(lp_rxn_id, rxn, sep = '; '))
    lines <- gsub(pattern = lp_rxn_id, replacement = rxn, lines)
  }
  
  for (i in rev(1:length(model@met_id))){
    met <- model@met_id[i]
    lp_met_id <- paste('r_', i, sep = '')
    # print(paste(lp_met_id, met, sep = '; '))
    lines <- gsub(lp_met_id, met, lines)
  }
  
  return(lines)
}

convert_sybil_to_lp <- function(sybil, output_file = "'temp.lp'", output = 'temp.lp'){
  optimizeProb(sybil, poCmd = list(c("writeProb", "LP_PROB", output_file, "'lp'")))
  new_lines <- rewrite_lp(output, sybil)
  write(new_lines, file = output)
}

convert_sybil_to_gurobi <- function(sybil, output_file = "'temp.lp'", output = 'temp.lp'){
  convert_sybil_to_lp(sybil, output_file = output_file, output = output)
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
    rev <- split_rev_rxns[i]
    
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
  lp_lines <- rewrite_lp(model = sybil_falcon_model)
  
  # model <- convert_sybil_to_gurobi(sybil_falcon_model)
  
  constr_idx <- grep('Subject To', lp_lines)
  lp_lines <- insert_into_list(lp_lines, new_constr, constr_idx)
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
  
  n <- length(model$varnames) #model$get_sizes()$NumVars
  
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
  
  vars <- model$varnames
  prev_obj <- model$obj
  model$obj <- rep(0, length(model$obj)) # clear the objective
  prev_sense <- model$modelsense
  
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
    set <- which(known_set_mtx[j,])
    
    coupled[i, set] <- TRUE
    coupled[set, set] <- TRUE
    active[set] <- FALSE
    
    list(coupled = coupled, active = active)
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
        # maximize
        sol <- solve_model(model, i, sense = 'max')
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        # model$setattr("Obj", setNames(0.0, vars[i]))
      }
      if (!near(model$lb[i], 0)){ #model$getattr("LB")[vars[i]] < (-1*tol_)
        # minimize
        sol <- solve_model(model, i, sense = 'min')
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_min[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        # model$setattr("Obj", setNames(0.0, vars[i]))
      }
      if (near(global_max[i], 0) & near(global_min[i], 0)){ #(abs(global_max[i]) < tol) & (abs(global_min[i]) < tol)
        blocked[i] <- TRUE
        active[i] <- FALSE
        next
      }
    }
    
    # set new bounds for selected rxn (temporarily)
    model$ub[i] <- fixed_val + 0.0*fix_tol_frac*abs(fixed_val)  #setattr("UB", setNames(fixed_val + 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    model$lb[i] <- fixed_val - 0.0*fix_tol_frac*abs(fixed_val)#setattr("LB", setNames(fixed_val - 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    
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
      
      # check for uncoupled via correlation
      if (cor_check){
        if (!correlation_check(flux, i, j)){next}
        
      }
      
      skip <- FALSE
      
      max <- 0
      min <- 0
      
      # model$setattr("Obj", setNames(1.0, vars[j]))
      if (!skip) {
        # maximize
        sol <- solve_model(model, j, sense = 'max')
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        sub_max <- pmax(sub_max, sol$x)
        sub_min <- pmin(sub_min, sol$x)
        
        flux[lp_calls%%stored_obs,] <- sol$x
        
        max <- sol$x[j]
        
        skip <- not_fixed(sub_max[j], sub_min[j])
      }
      
      if (!skip) {
        # minimize
        sol <- solve_model(model, j, sense = 'min')
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        sub_max <- pmax(sub_max, sol$x)
        sub_min <- pmin(sub_min, sol$x)
        
        flux[lp_calls%%stored_obs,] <- sol$x
        
        min <- sol$x[j]
        
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
      
      # model$setattr("Obj", setNames(0.0, vars[j]))
    }
    
    
    # unfix i
    model$ub[i] <- prev_ub
    model$lb[i] <- prev_lb
    # model$setattr("UB", prev_ub)
    # model$setattr("LB", prev_lb)
  }
  
  model$obj <- prev_obj
  model$modelsense <- prev_sense
  # model$setattr("Obj", prev_obj)
  # model$setattr("ModelSense", prev_sense)
  
  print(lp_calls)
  
  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}

data("Ec_core")
model <- convert_sybil_to_gurobi(Ec_core)
# falcon_model <- GRB_generate_falcon_model(Ec_core)

r_mtx <- flux_coupling_raptor(model)
r_sets <- get_list_of_sets_from_mtx(r_mtx)