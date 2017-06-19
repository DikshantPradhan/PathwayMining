library(sybilcycleFreeFlux)
data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=500
nPnts=5000

lm_fitting <- function(model, rxn_idx){
  sample_df = sampler(model)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample_df2 = sampler(suppressed_model(model, rxn_idx)) # phosphate exchange
  
  sample_df <- rbind.data.frame(sample_df, sample_df2)
  
  #output <- logical(length = nrow(sample_df))
  #output <- !output
  fit <- lm(sample_df[,rxn_idx] ~ ., data = sample_df)
  return(fit)
}

lm_fitting_probit <- function(sample_df, rxn_idx){
  fit <- glm(sample_df[,rxn_idx] ~ ., family = binomial(link = "probit"), data = sample_df, maxit = 100)
  return(fit)
}

lm_fitting_logit <- function(sample_df, rxn_idx){
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

flux_coupling <- function(sample, binary = TRUE){
  rxn_ct = ncol(sample)
  
  coupling = array(0, dim = c(2, rxn_ct, rxn_ct), 
                   dimnames = list(c("coupling ratio", "r-squared"), colnames(sample), colnames(sample))) # layer 1 is coupling ratio, layer 2 is r-squared
  
  for(i in 1:rxn_ct){
    for(j in i:rxn_ct){
      if (binary){
        fit <- glm(sample[,i] ~ sample[,j], family = "binomial", data = sample_df, maxit = 100)
      }
      else{
        fit <- lm(sample[,i] ~ sample[,j], data = sample)
      }
      sfit <- summary(fit)
      coupling[1,i,j] <- sfit$coefficients[2] #fit$coefficients[[2]]
      coupling[2,i,j] <- sfit$coefficients[4] #summary(fit)[[8]]
      # fit <- lm(sample[,i] ~ sample[,j], data = sample)
      # coupling[1,i,j] <- fit$coefficients[[2]]
      # coupling[2,i,j] <- summary(fit)[[8]]
    }
  }
  #colnames(coupling) <- colnames(sample)
  #rownames(coupling) <- colnames(sample)
  coupling[1, ,] <- coupling_generalize(coupling[1, ,])
  return(coupling)
}

flux_coupling_specific <- function(sample, rxn_idx, binary = TRUE){ # samples, # of rxn to compare with all others
  rxn_ct = ncol(sample)
  
  coupling = array(0, dim = c(2, rxn_ct)) # row 1 is coupling ratio, row 2 is standard error
  
  for(i in 1:rxn_ct){
    #fit <- lm(sample[,rxn_idx] ~ sample[,i], data = sample)
    if (binary){
      fit <- glm(sample[,rxn_idx] ~ sample[,i], family = "binomial", data = sample_df, maxit = 100)
    }
    else{
      fit <- lm(sample[,rxn_idx] ~ sample[,i], data = sample)
    }
    sfit <- summary(fit)
    coupling[1,i] <- sfit$coefficients[2] #fit$coefficients[[2]]
    coupling[2,i] <- sfit$coefficients[4] #summary(fit)[[8]]
  }
  colnames(coupling) <- colnames(sample)
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
  
  colnames(coupling) <- colnames(sample)
  
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

coupling_generalize <- function(array){
  
  for (i in 1:nrow(array)){
    for (j in 1:ncol(array)){
      if (array[i,j] > 0.7){
        array[i,j] = 1
      }
      else if (array[i,j] > 0.2){
        array[i,j] = 0.5
      }
      else if (array[i,j] < -0.2){
        array[i,j] = -0.5
      }
      else if (array[i,j] < -0.7){
        array[i,j] = -1
      }
      else {
        array[i,j] = 0
      }
    }
  }
  
  return(array)
}

#coupling_change_analysis <- function()

sampler <- function(model){
  sample = ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample = t(sample$Points)
  colnames(sample) <- model@react_id
  sample_df <- as.data.frame(sample)
  return(sample_df)
}

sampler_lm_fitting <- function(model, rxn_idx, binary = FALSE){
  sample_df = sampler(model)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample_df2 = sampler(suppressed_model(model, rxn_idx)) # phosphate exchange
  
  if (binary){
    sample_df[,rxn_idx] = rep(1, nrow(sample_df)) 
  }
  
  sample_df <- rbind.data.frame(sample_df, sample_df2)
  
  return(sample_df["Biomass_Ecoli_core_w_GAM" > 0])
}

rescale_sample <- function(sample, rxn_idx = 0){
  rxn_list = c(1:ncol(sample))
  if (rxn_idx != 0){
    rxn_list <- rxn_list[-rxn_idx]
  }
  
  #print(rxn_list)
  
  for(i in rxn_list){
    sample[,i] <- rescale(sample[,i], c(-1,1))
  }
  
  return(sample)
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

# sample <- sampler_lm_fitting(model, 37, binary = TRUE)
# coupling1 <- flux_coupling_specific(rescale_sample(sample,37), 37, binary = TRUE)

sample_og <- sampler(model)
#sample_og[,37] = rep(1, nrow(sample_og))
sample_suppr <- sampler(suppressed_model(model, 37))

coupling_og <- flux_coupling(rescale_sample(sample_og), binary = FALSE)

coupling_suppr <- flux_coupling(rescale_sample(sample_suppr), binary = FALSE)
#coupling_suppr <- flux_coupling(rescale_sample(sample_suppr, 37))