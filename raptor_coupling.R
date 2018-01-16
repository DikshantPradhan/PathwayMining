# new coupling protocol adapted from jensenlab/raptor

library('raptor')

# SUB-ROUTINES

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

# MAIN FUNCTION

flux_coupling_raptor <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01,
      bnd_tol = 0.1, stored_obs = 4000, cor_iter = 3, reaction_indexes = c()) {

  # min_fva_cor is minimum correlation between fluxes
  # bnd_tol is allowed error in comparing max & min flux
  # fix_frac is const value used in fixing flux at non-zero value
  # fix_tol_frac is error allowed in determining whether flux is fixed
  # stored_obs is # not flux values to be stored
  # cor_iter is number of iterations after which correlation is considered in checking coupling
  
  n <- model$get_sizes()$NumVars

  # if empty set, then assume all reactions are to be inspected
  if (length(reaction_indexes) == 0){
    reaction_indexes <- c(1:n)
  }

  vars <- model$get_names()$VarName
  prev_obj <- model$getattr("Obj")
  model$setattr("Obj", setNames(numeric(n), vars)) # clear the objective
  prev_sense <- model$getattr("ModelSense")

  global_max <- rep(0, n)
  global_min <- rep(0, n)

  coupled <- matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars))

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

  flux <- matrix(c(0), nrow = stored_obs, ncol = n)
  lp_calls <- 0

  update_flux <- function(flux_, idx, sol){
    if (stored_obs > 0){
      flux_[idx,] <- sol
    }
    return(flux_)
  }

  for (idx in 1:(length(reaction_indexes)-1)) { # (i in 1:(n-1)) 
    # iterate over passed in idxs instead (idx in 1:length(reaction_indexes)); i <-  reaction_indexes[idx]
    i <-  reaction_indexes[idx]

    if (!active[i] | blocked[i]) next
    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)
    prev_ub <- model$getattr("UB")[vars[i]]
    prev_lb <- model$getattr("LB")[vars[i]]

    if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
      fixed_val <- rxn_fix(global_max[i], global_min[i])
    }
    else {
      if (!near(model$getattr("UB")[vars[i]], 0)){ #model$getattr("UB")[vars[i]] > tol

        model$setattr("Obj", setNames(1.0, vars[i]))
        model$set_model_sense(maximize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()

        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)

        flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)

        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }

        model$setattr("Obj", setNames(0.0, vars[i]))
      }
      if (!near(model$getattr("LB")[vars[i]], 0)){ #model$getattr("LB")[vars[i]] < (-1*tol_)

        model$setattr("Obj", setNames(1.0, vars[i]))
        model$set_model_sense(minimize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()

        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)

        flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)

        if (!near(global_min[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }

        model$setattr("Obj", setNames(0.0, vars[i]))
      }
      if (near(global_max[i], 0) & near(global_min[i], 0)){ #(abs(global_max[i]) < tol) & (abs(global_min[i]) < tol)
        print(paste('blocked:', i))
        blocked[i] <- TRUE
	      active[i] <- FALSE
        next
      }
    }

    # set new bounds for selected rxn (temporarily)
    model$setattr("UB", setNames(fixed_val + 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
    model$setattr("LB", setNames(fixed_val - 0.0*fix_tol_frac*abs(fixed_val), vars[i]))

    # couple reaction to itself if not blocked
    if (!blocked[i]){
      coupled[i,i] <- TRUE
      active[i] <- FALSE
    }

    for (idx2 in (idx+1):length(reaction_indexes)) { # (j in (i+1):n) 
      # also keep this in passed in idxs (idx2 in (idx+1):length(reaction_indexes)); j <-  reaction_indexes[idx2]
      j <-  reaction_indexes[idx2]
      #if (j == 137){print('137')}
      # check for fixed or blocked
      if (!active[j] | blocked[j]) next
      # check for uncoupled
      if ((stored_obs > 0) & (i > cor_iter)){
        C <- cor(flux[,i], flux[,j])
        if (is.na(C) | (abs(C) < min_fva_cor)){next}
      }

      if (not_fixed(sub_max[j], sub_min[j])) next

      skip <- FALSE

      max <- 0
      min <- 0

      model$setattr("Obj", setNames(1.0, vars[j]))
      if (!skip) {
        model$set_model_sense(maximize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()

        flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)

        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)

        max <- sol$X[j]

        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (!skip) {
        model$set_model_sense(minimize=TRUE)
        model$optimize()
        lp_calls <- lp_calls + 1
        sol <- model$get_solution()

        flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)

        global_max <- pmax(global_max, sol$X)
        global_min <- pmin(global_min, sol$X)
        sub_max <- pmax(sub_max, sol$X)
        sub_min <- pmin(sub_min, sol$X)

        min <- sol$X[j]

        skip <- not_fixed(sub_max[j], sub_min[j])
      }

      if (near(max, 0) & near(min, 0)){skip = TRUE}

      if (!skip) { # finally label as coupled
        coupled[i,j] <- TRUE
        active[j] <- FALSE
      }

      # print(lp_calls)

      model$setattr("Obj", setNames(0.0, vars[j]))
    }


    # unfix i
    model$setattr("UB", prev_ub)
    model$setattr("LB", prev_lb)
  }

  model$setattr("Obj", prev_obj)
  model$setattr("ModelSense", prev_sense)

  print(lp_calls)

  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}

#ecoli <- as_GRBmodel(model)
#ecoli$show_output(FALSE)
#print(flux_coupling_raptor(ecoli))
#print(get_list_of_sets(return_couples(flux_coupling_raptor(ecoli)$coupled)))
