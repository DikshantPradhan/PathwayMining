library(sybilcycleFreeFlux)
data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=500
nPnts=1000

lm_fitting <- function(model, rxn_idx){
  sample_df = sampler(model)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample_df2 = sampler(suppressed_model(model, rxn_idx)) # phosphate exchange
  
  sample_df <- rbind.data.frame(sample_df, sample_df2)
  
  #output <- logical(length = nrow(sample_df))
  #output <- !output
  fit <- lm(sample_df[,rxn_idx] ~ ., data = sample_df)
  return(fit)
}

lm_fitting_probit <- function(model, rxn_idx){
  sample_df = sampler(model)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample_df2 = sampler(suppressed_model(model, rxn_idx)) # phosphate exchange
  
  sample_df[,rxn_idx] = rep(1, nrow(sample_df))
  
  sample_df <- rbind.data.frame(sample_df, sample_df2)
  
  #output <- logical(length = nrow(sample_df))
  #output <- !output
  fit <- glm(sample_df[,rxn_idx] ~ ., family = binomial(link = "probit"), data = sample_df, maxit = 100)
  return(fit)
}

lm_fitting_logit <- function(model, rxn_idx){
  sample_df = sampler(model)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample_df2 = sampler(suppressed_model(model, rxn_idx)) # phosphate exchange
  
  sample_df[,rxn_idx] = rep(1, nrow(sample_df))
  
  sample_df <- rbind.data.frame(sample_df, sample_df2)
  
  #output <- logical(length = nrow(sample_df))
  #output <- !output
  fit <- glm(sample_df[,rxn_idx] ~ ., family = "binomial", data = sample_df, maxit = 100)
  return(fit)
}

flux_comparison <- function(model, rxn_idx){
  sample_df = sampler(model)
  sample_df2 = sampler(suppressed_model(model, rxn_idx))
  
  sample_diff = c()
  
  for (i in 1:ncol(sample_df)){
    diff = mean(sample_df2[,i] - sample_df[,1])
    sample_diff = c(sample_diff, diff)
  }
  
  #colnames(sample_diff) <- model@react_id
  
  return(sample_diff)
}

flux_coupling <- function(sample){
  rxn_ct = ncol(samples)
  
  coupling = array(0, dim = c(2, rxn_ct, rxn_ct)) # layer 1 is coupling ratio, layer 2 is r-squared
  
  for(i in 1:rxn_ct){
    for(j in i:rxn_ct){
      fit <- lm(sample[,i] ~ sample[,j], data = sample)
      coupling[1,i,j] <- fit$coefficients[[2]]
      coupling[2,i,j] <- summary(fit)[[8]]
    }
  }
  
  return(coupling)
}

flux_coupling_fcf <- function(sample){
  rxn_ct = ncol(sample)
  obs_ct = nrow(sample)
  
  coupling = array(0, dim = c(rxn_ct, rxn_ct)) # layer 1 is coupling ratio, layer 2 is r-squared
  
  for(i in 1:rxn_ct){
    for(j in i:rxn_ct){
      R <- fcf_R(sample, i, j)
      #print(length(R))
      #print(length(coupling[ ,i,j]))
      coupling[i,j] <- fcf_analysis(R)
    }
  }
  
  return(coupling)
}

fcf_analysis <- function(R){
  if (is.nan(min(R)) | is.nan(max(R))){
    return(-1)
  }
  if (min(R) == 0 & max(R) < Inf){
    return(1) # directional
  }
  if (min(R) < Inf & max(R) < Inf){
    return(2) # partial
  }
  if (min(R) == max(R)){
    return(3) # fully
  }
  if (max(R) == Inf){
    if (min(R) == 0){
      return(4) # uncoupled
    }
    if (min(R) > 0 & min(R) < Inf){
      return(1) # directional
    }
  }

  return(-1)
}

fcf_R <- function(sample, v1, v2){
  sample_num = nrow(sample)
  R <- rep(0, sample_num)
  
  for (i in 1:sample_num){
    R[i] = sample[i,v1]/sample[i,v2]
  }
  
  return(R)
}

sampler <- function(model){
  sample = ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample = t(sample$Points)
  colnames(sample) <- model@react_id
  sample_df <- as.data.frame(sample)
  return(sample_df)
}

suppressed_model <- function(model, rxn_idx){
  model@lowbnd[rxn_idx] <- 0
  model@uppbnd[rxn_idx] <- 0
  return(model)
}

main <- function(){
  # fit <- lm_fitting(model)
  # 
  # model2 <- suppressed_model(model = model, 54)
  # 
  # fit2 <- lm_fitting(model2)
  
  sample <- sampler(model)
  coupling <- flux_coupling(sample)
}

#main()

# sample <- sampler(model)
# coupling <- flux_coupling_fcf(sample)

# fit <- lm_fitting(model, 37)
# sample_diff <- flux_comparison(model, 37)