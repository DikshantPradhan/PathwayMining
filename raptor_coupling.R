# new coupling protocol adapted from jensenlab/raptor

# SUB-ROUTINES

update_flux <- function(global, newflux, max){
  for (i in 1:length(global)){
    if (max == TRUE){
      global[i] = max(global[i], newflux[i])
    }
    if (max == FALSE){
      global[i] = min(global[i], newflux[i])
    }
  }

  return(global)
}



# MAIN FUNCTION
flux_coupling <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01) {
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
