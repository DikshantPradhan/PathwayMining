# new coupling protocol adapted from jensenlab/raptor

# SUB-ROUTINES

update_flux <- function(global, newflux, max){
  if (max == TRUE){
    global[i] = pmax(global, newflux)
  }
  if (max == FALSE){
    global[i] = pmin(global, newflux)
  }

  #for (i in 1:length(global)){
  #  if (max == TRUE){
  #    global[i] = max(global[i], newflux[i])
  #  }
  #  if (max == FALSE){
  #    global[i] = min(global[i], newflux[i])
  #  }
  #}

  return(global)
}

set_model_bounds <- function(model, rxn, UB, LB){
  model$setattr("UB", setNames(UB, rxn))
  model$setattr("LB", setNames(LB, rxn))

  return(model)
}

optimize_rxn <- function(model, rxn, max){
  model$setattr("Obj", setNames(1.0, rxn))
  if (max){
    model$set_model_sense(maximize=TRUE)
  }
  else {
    model$set_model_sense(minimize=TRUE)
  }
  model$optimize()
  sol <- model$get_solution()

  return(sol)
}

rxn_fix <- function(max, min){
  #max <- optimize_rxn(model, rxn, max = TRUE)$X[rxn]
  #min <- optimize_rxn(model, rxn, max = FALSE)$X[rxn]
  if (is.infinite(max)){
    max = 1000
  }
  if (is.infinite(min)){
    min = -1000
  }
  avg <- mean(c(max, min))
  if (avg == 0){
    avg <- mean(c(avg, max))
  }

  return(avg)
}

# MAIN FUNCTION
flux_coupling_raptor_1 <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01) {
  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  prev_obj <- model$getattr("Obj")
  model$setattr("Obj", setNames(numeric(n), vars)) # clear the objective
  # prev_sense <- model$getattr("ModelSense")

  not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }

  set_list <- vars

  global_max <- rep(-Inf, n)
  global_min <- rep(Inf, n)
  coupled <- rep(FALSE, n)

  ## test this section (initialization) ##
  start_idx <- 1
  start = FALSE
  while(!start){
    fixed_val <- rxn_fix(model, vars[start_idx])
    if (fixed_val == 0){
      start_idx = start_idx + 1
    }
    else {
      start = TRUE
    }
  }

  global_max[start_idx] <- optimize_rxn(model, vars[start_idx], max = TRUE)$X[vars[start_idx]]
  global_min[start_idx] <- optimize_rxn(model, vars[start_idx], max = FALSE)$X[vars[start_idx]]
  ##

  for (i in start_idx:(n-1)){
    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)

    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]

    # set new bounds for selected rxn (temporarily)
    model$setattr("UB", setNames(fixed_val + 0.5*fix_tol_frac*abs(fixed_val), vars[i]))
    model$setattr("LB", setNames(fixed_val - 0.5*fix_tol_frac*abs(fixed_val), vars[i]))

    # fixed_val <- rxn_fix(model, vars[i])
    # if (fixed_val == 0) next

    fixed_val <- mean(c(global_max[i], global_min[i]))
    if (fixed_val == 0){
      fixed_val <- mean(c(fixed_val, global_max[i]))
    }

    #model_i <- set_model_bounds(model, vars[i], UB = fixed_val + 0.5*fix_tol_frac*abs(fixed_val),
    #  LB = fixed_val - 0.5*fix_tol_frac*abs(fixed_val))

    for (j in (i+1):n){
      if (coupled[j] == TRUE) next

      skip <- FALSE
      model$setattr("Obj", setNames(1.0, vars[j]))
      if (!skip) {
        model$set_model_sense(maximize=TRUE)
        model$optimize()

        sol <- model$get_solution()
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) {
        model$set_model_sense(minimize=TRUE)
        model$optimize()

        sol <- model$get_solution()
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) {
        #coupled[i,j] <- TRUE
        coupled[j] <- TRUE
        set_list <- group_sets(set_list, vars[i], vars[j])
        #active[j] <- FALSE
      }

      # max_sol <- optimize_rxn(model, vars[j], max = TRUE)$X
      # min_sol <- optimize_rxn(model, vars[j], max = FALSE)$X

      global_max <- pmax(global_max, sub_max)
      global_min <- pmin(global_min, sub_min)

      #if (!not_fixed(max_sol[vars[j]], min_sol[vars[j]])){
      #  print(c(vars[i], vars[j]))
      #  set_list <- group_sets(set_list, vars[i], vars[j])
      #  coupled[i] = TRUE
      #  coupled[j] = TRUE

      #  if (j < n){ # epilogue
      #    for (k in (j+1):n){
      #      if (!not_fixed(max_sol[vars[k]], min_sol[vars[k]])){
      #        ## ???

      #        max_sol <- optimize_rxn(model, vars[k], max = TRUE)$X
      #        min_sol <- optimize_rxn(model, vars[k], max = FALSE)$X

      #        if (!not_fixed(max_sol[vars[k]], min_sol[vars[k]])){
      #          print(c(vars[i], vars[k]))
      #          set_list <- group_sets(set_list, vars[i], vars[k])
      #          coupled[k] = TRUE

      #          global_max <- pmax(global_max, max_sol)
      #          global_min <- pmin(global_min, min_sol)
      #        }
      #      }
      #    }
      #  }


      #}

    }

    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)
  }

  print(global_max)
  print(global_min)
  return(set_list)
}

flux_coupling_raptor_2 <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01) {
  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  prev_obj <- model$getattr("Obj")
  model$setattr("Obj", setNames(numeric(n), vars)) # clear the objective
  # prev_sense <- model$getattr("ModelSense")

  global_max <- rep(-Inf, n)
  global_min <- rep(Inf, n)

  blocked <- rep(FALSE, n)
  active <- rep(TRUE, n)

  set_list <- vars

  not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }

  flux <- data.frame()
  flux <- rbind(flux, rep(0, n))
  flux_idx <- 1

  for (i in 1:(n-1)){

    if (blocked[i] | !active[i]){print("bad rxn"); next}

    if (global_max[i] != 0 | global_min[i] != 0){
      fixed_val <- rxn_fix(global_max[i], global_min[i])
    }
    else {
      if (model$getattr("UB")[vars[i]] > 0){

        sol <- optimize_rxn(model, vars[i], max = TRUE)
        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)

        flux <- rbind(flux, sol$X)
        flux_idx <- flux_idx +1
      }
      if (model$getattr("LB")[vars[i]] < 0){

        sol <- optimize_rxn(model, vars[i], max = FALSE)
        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)

        flux <- rbind(flux, sol$X)
        flux_idx <- flux_idx +1
      }
      if (global_max[i] == 0 & global_min[i] == 0){
        blocked[i] <- TRUE
        next
      }
    }

    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]

    fixed_val <- rxn_fix(global_max[i], global_min[i])
    #print(fixed_val)
    model$setattr("UB", setNames(fixed_val + 0.5*fix_tol_frac*abs(fixed_val), vars[i]))
    model$setattr("LB", setNames(fixed_val - 0.5*fix_tol_frac*abs(fixed_val), vars[i]))

    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)

    for (j in (i+1):n){
      if (blocked[j] | !active[j]){print("bad rxn"); next}
      if (not_fixed(sub_min[j], sub_max[j])) next
      C <- cor(flux[,i], flux[,j])
      #print(C)
      if (!is.na(C) & C < min_fva_cor){next}
      skip <- FALSE
      model$setattr("Obj", setNames(1.0, vars[j]))
      if (!skip) {
        model$set_model_sense(maximize=TRUE)
        model$optimize()
        sol <- model$get_solution()

        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)

        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)

        skip <- not_fixed(sub_max[j], sub_min[j])

        flux <- rbind(flux, sol$X)
        flux_idx <- flux_idx +1
      }

      if (!skip) {
        model$set_model_sense(minimize=TRUE)
        model$optimize()
        sol <- model$get_solution()

        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)

        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])

        flux <- rbind(flux, sol$X)
        flux_idx <- flux_idx +1
      }

      if (!skip) {
        # coupled[i,j] <- TRUE
        print(c(i, j))
        #print(vars[i], vars[j])
        set_list <- group_sets(set_list, vars[i], vars[j])
        active[j] <- FALSE
      }
    }

    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)

  }

  #print(vars)
  #print(flux_idx)
  #print(flux)
  return(set_list)
}

ecoli <- as_GRBmodel(model)
ecoli$show_output(FALSE)
print(flux_coupling_raptor_2(ecoli))

flux_coupling_raptor_test <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01) {
  n <- model$get_sizes()$NumVars
  vars <- model$get_names()$VarName
  prev_obj <- model$getattr("Obj")
  model$setattr("Obj", setNames(numeric(n), vars)) # clear the objective
  prev_sense <- model$getattr("ModelSense")

  fva <- flux_variability(model, obj_frac=NA, return_fluxes=TRUE) # NEED TO REPLACE THIS
  use_min_fva_cor <- min_fva_cor != 0.0
  if (use_min_fva_cor) {
    fva_cor <- cor(fva$fluxes)
  }

  coupled <- matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars))
  fixed <- near(fva$minflux, fva$maxflux)
  blocked <- near(fva$minflux, 0) & near(fva$maxflux, 0)
  active <- !(fixed | blocked)

  not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }

  lp_calls <- 0
  for (i in 1:(n-1)) {
    if (!active[i]) next
    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)
    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]
    fixed_val <- fva$minflux[i] + fix_frac*(fva$maxflux[i] - fva$minflux[i])

    # set new bounds for selected rxn (temporarily)
    model$setattr("UB", setNames(fixed_val + 0.5*fix_tol_frac*abs(fixed_val), vars[i]))
    model$setattr("LB", setNames(fixed_val - 0.5*fix_tol_frac*abs(fixed_val), vars[i]))

    for (j in (i+1):n) {
      # check for fixed or blocked
      if (!active[j]) next
      # check for uncoupled
      if (use_min_fva_cor && !is.na(fva_cor[i,j]) && abs(fva_cor[i,j]) < min_fva_cor) next
      if (not_fixed(sub_max[j], sub_min[j])) next

      skip <- FALSE
      model$setattr("Obj", setNames(1.0, vars[j]))
      if (!skip) {
        model$set_model_sense(maximize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) {
        model$set_model_sense(minimize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)
        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) { # finally label as coupled
        coupled[i,j] <- TRUE
        active[j] <- FALSE
      }

      model$setattr("Obj", setNames(0.0, vars[j]))
    }
    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)
  }

  model$setattr("Obj", prev_obj)
  model$setattr("ModelSense", prev_sense)

  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}
