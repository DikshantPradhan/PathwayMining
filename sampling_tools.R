library(sybilcycleFreeFlux)
library(plotrix)

lm_fitting <- function(model, rxn_idx){
  sample_df = sampler(model)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  sample_df2 = sampler(suppressed_model(model, rxn_idx))
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
  
  return(sample_diff)
}

flux_coupling <- function(sample, binary = TRUE){
  rxn_ct = ncol(sample)
  coupling = array(0, dim = c(rxn_ct, rxn_ct), 
                   dimnames = list(colnames(sample), colnames(sample))) # layer 1 is coupling ratio, layer 2 is r-squared
  
  for(i in 1:rxn_ct){
    for(j in i:rxn_ct){
      if (binary){
        fit <- glm(sample[,i] ~ sample[,j], family = "binomial", data = sample, maxit = 100)
      }
      else{
        fit <- lm(sample[,i] ~ sample[,j], data = sample)
      }
      sfit <- summary(fit)
      coupling[i,j] <- sfit$coefficients[2]
    }
  }
  #coupling <- coupling_generalize(coupling)
  return(coupling)
}

flux_coupling_cor <- function(sample){
  rxn_ct = ncol(sample)
  coupling = array(0, dim = c(rxn_ct, rxn_ct), 
                   dimnames = list(colnames(sample), colnames(sample)))
  
  for(i in 1:rxn_ct){
    for(j in i:rxn_ct){
      cor <- cor(sample[,i], sample[,j])
      if (is.na(cor)){
        cor <- 0
      }
      coupling[i,j] <- cor
    }
  }
  #coupling <- coupling_generalize(coupling)
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
    coupling[1,i] <- sfit$coefficients[2]
    coupling[2,i] <- sfit$coefficients[4]
  }
  colnames(coupling) <- colnames(sample)
  return(coupling)
}

flux_coupling_fcf <- function(sample){
  rxn_ct = ncol(sample)
  obs_ct = nrow(sample)
  coupling = array(0, dim = c(rxn_ct, rxn_ct))
  
  for(i in 1:rxn_ct){
    for(j in i:rxn_ct){
      R <- fcf_R(sample, i, j)
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
      if (is.na(array[i,j])){
        # do nothing
      }
      else if (array[i,j] > 0.7){
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

return_couples <- function(array){ # correlation array (output from flux_coupling_???)
  
  #couple_list <- c()
  row <- dimnames(array)[[1]]
  col <- dimnames(array)[[2]]
  
  # rxn_col1 <- c("Biomass_Ecoli_core_w_GAM")
  # rxn_col2 <- c("Biomass_Ecoli_core_w_GAM")
  
  rxn_col1 <- c()
  rxn_col2 <- c()
  
  # print(dim(array))
  
  for (i in 1:dim(array)[1]){
    for (j in 1:dim(array)[2]){
      # print(array[i,j])
      if (is.na(array[i,j])){
        # do nothing
      }
      else if (abs(array[i,j]) > 0.99){
        # print("ok")
        #couple_list <- c(couple_list, paste(row[i],col[j], sep = "__"))
        rxn_col1 <- c(rxn_col1, row[i])
        rxn_col2 <- c(rxn_col2, col[j])
      }
    }
  }
  
  # if (length(rxn_col1) == 0){
  #   rxn_col1 <- c(rxn_col1, "Biomass_Ecoli_core_w_GAM")
  #   rxn_col2 <- c(rxn_col2, "Biomass_Ecoli_core_w_GAM")
  # }
  
  rxns <- cbind(rxn_col1, rxn_col2)
  colnames(rxns) <- c("rxn1", "rxn2")
  return(rxns)
}

return_coupling_change <- function(og_array, suppr_array, couples_list){
  
  row <- dimnames(og_array)[[1]]
  col <- dimnames(og_array)[[2]]
  #changes <- c()
  
  rxn_col1 <- c()
  rxn_col2 <- c()
  old_cor <- c()
  new_cor <- c()
  
  for (i in 1:nrow(couples_list)){
    rxns = couples_list[i,] #strsplit(i, split = "__")[[1]]
    idx_i = which(row == rxns[1])
    idx_j = which(col == rxns[2])
    change_string <- paste(rxns[1], " & ", rxns[2], ": ", og_array[idx_i, idx_j], " --> ", suppr_array[idx_i, idx_j])
    print(change_string)
    
    rxn_col1 <- c(rxn_col1, rxns[1])
    rxn_col2 <- c(rxn_col2, rxns[2])
    old_cor <- c(old_cor, as.numeric(og_array[idx_i, idx_j]))
    new_cor <- c(new_cor, as.numeric(suppr_array[idx_i, idx_j]))
    #changes <- c(changes, change_string)
  }
  
  old_cor <- as.numeric(old_cor)
  new_cor <- as.numeric(new_cor)
  
  rxns <- cbind(rxn_col1, rxn_col2, old_cor, new_cor)
  colnames(rxns) <- c("rxn1", "rxn2", "old_cor", "new_cor")
  return(rxns)
}

coupling_change <- function(og_array, suppr_array){
  
  couples <- return_couples(suppr_array - og_array)
  d_couples <- return_coupling_change(og_array, suppr_array, couples)
  
  return(d_couples)
}

find_coupling_change <- function(sample_og, sample_suppr){
  coupling_og <- flux_coupling_cor(sample_og)
  coupling_suppr <- flux_coupling_cor(sample_suppr)
  
  return(coupling_change(coupling_og, coupling_suppr))
}

sampler <- function(model, W=200, nPnts=5000, steps=1, Biomass = FALSE, Floor = TRUE){
  sample = ACHR(model,W,nPoints=nPnts,stepsPerPoint=steps)
  sample = t(sample$Points)
  colnames(sample) <- model@react_id
  sample_df <- as.data.frame(sample)
  if (Floor){
    for (i in 1:nrow(sample_df)){
      for (j in 1:ncol(sample_df)){
        if (abs(sample_df[i,j]) < 1.0e-9){
          sample_df[i,j] <- 0
        }
      }
    }
  }
  if (Biomass){
    # sample_df <- sample_df[which(sample_df[,13] > 0),]
  }
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

rescale_sample <- function(sample, model, rxn_idx = 0){ # rxn_idx is idx of blocked reaction, 0 if nothing is blocked
  
  opt <- fluxVar(model, percentage = 99)
  model_fva <- opt@lp_obj
  fva_min <- model_fva[1:95]
  fva_max <- model_fva[96:190]
  
  rxn_list = c(1:ncol(sample))
  if (rxn_idx != 0){
    rxn_list <- rxn_list[-rxn_idx]
  }
  
  for(i in rxn_list){ # rescale to flux variability range
    rescaled <- rescale(c(fva_min[i], sample[,i], fva_max[i]), c(-1,1))
    sample[,i] <- rescaled[-c(1, length(rescaled))]
  }
  
  return(sample)
}

suppressed_model <- function(model, rxn_idx){
  # model@lowbnd[rxn_idx] <- 0
  # model@uppbnd[rxn_idx] <- 0
  if (rxn_idx == 0){
    return(model)
  }
  
  for (i in rxn_idx){
    model <- changeBounds(model, rxn_idx, lb = 0, ub = 0)
  }
  
  return(model)
}

maxDiff_dist <- function(sample, model){
  
  rescaled <- rescale_sample(sample, model = model)
  maxDiff <- array(0, dim = c(nrow(sample)-1, ncol(sample)))
  
  for (i in 1:ncol(rescaled)){
    sorted <- sort(rescaled[,i])
    max = 0
    for (j in 1:(length(sorted)-1)){
      temp_max <- sorted[j+1] - sorted[j]
      maxDiff[j,i] <- temp_max
    }
    # maxDiff <- c(maxDiff, max)
  }
  
  colnames(maxDiff) <- colnames(sample)
  return(maxDiff)
}
