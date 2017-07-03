library(tree)

get_set_idx <- function(rxn, rxns_list){
  idx <- grep(core_rxn_id(rxn), rxns_list)
  
  for (j in idx){
    if (rxn %in% rxns_list[[j]]){
      # idx <- c()
      # idx <- j
      return(j)
    }
  }
  return(integer(0))
}

check_sets_for_containing <- function(rxns, set_list){
  
  for (i in 1:length(set_list)){
    if (all(rxns %in% set_list[[i]])){
      return(TRUE)
    }
  }
  
  return(FALSE)
}

check_set_list_for_containing <- function(rxns, set_lists){
  lists <- c()
  
  for (i in 1:length(set_lists)){
    if (check_sets_for_containing(rxns, set_lists[[i]])){
      lists <- c(lists, i)
    }
  }
  
  return(lists)
}

find_all_sets_for_rxn <- function(rxn_id, set_lists){
  for (i in 1:95){
    print(paste(i, ":"))
    print(set_lists[[i]][get_set_idx(rxn_id, set_lists[[i]])])
  }
}

correlation_predictor <- function(sample_df, rxn1, rxn2){
  cor_list <- rep(FALSE, 95)
  cor_list[check_set_list_for_containing(c(rxn1, rxn2), set_lists)] <- TRUE
  
  correlation_list <- c()
  
  
  for (j in which(sample_df["block"] == 0)){
    correlation_list[j] <- check_sets_for_containing(c(rxn1, rxn2), og_set_list)
  }

  for (i in 1:95){
    for (j in which(sample_df["block"] == i)){
      correlation_list[j] <- check_sets_for_containing(c(rxn1, rxn2), set_lists[[i]]) #cor_list[i]
    }
  }
  
  correlation_list <- as.factor(correlation_list)
  
  # sample_df <- data.frame(sample_df)
  sample <- sample_df[-c(96)]
  fit <- tree(correlation_list ~ ., data = sample)
  
  plot(fit)
  text(fit)
  return(fit)
}

predictor_matrix_generator <- function(sample_df){ # need to complete this function
  
  predictors <- array()
  
  for (i in 1:95){
    for (j in i:95){
      predictors[i,j] <- correlation_predictor()
    }
  }
  
  return(predictors)
}

pair_lists_from_predictors <- function(predictors){
  
  pairs <- c()
  
  for (i in 1:95){
    for (j in i:95){
      prediction <- c()
      if (prediction){
        # add pair to list
      }
    }
  }
  
  return(pairs)
}

sample_df_generator <- function(nPnts = 500, steps = 10){
  sample_df = sampler(model, nPnts = nPnts, steps = steps)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)

  sample_df["block"] <- rep(0, nrow(sample_df))

  for (i in 1:95){
    sample_df2 = sampler(suppressed_model(model, i), nPnts = nPnts, steps = steps)
    sample_df2["block"] <- rep(i, nrow(sample_df2))
    sample_df <- rbind.data.frame(sample_df, sample_df2)
  }
  sample_df <- data.frame(sample_df)
  return(sample_df)
}

test_fit <- function(fit, rxn1, rxn2){
  report <- c()
  
  # data <- data.frame(sampler(suppressed_model(model, i), nPnts = 10))
  # 
  # df <- rep(0, 95)
  # for (i in 1:95){
  #   df[i] <- mean(data[,i])
  # }
  # names(df) <- colnames(data)
  
  # print(df)
  
  df <- data <- data.frame(sampler(suppressed_model(model, i), nPnts = 1))
  
  for (i in 1:95){
    data <- data.frame(sampler(suppressed_model(model, i), nPnts = 100))
    
    # df <- rep(0, 95)
    for (j in 1:95){
      df[j] <- mean(data[,j])
    }
    # names(df) <- colnames(data)
    
    result <- predict(fit, newdata = data.frame(df))
    prediction <- colnames(result)[which.max(result)]
    known <- check_sets_for_containing(c(rxn1, rxn2), set_lists[[i]])
    
    # report <- c(report, paste(i, ":", prediction, known))
    
    if (prediction == "TRUE" & !known){
      report <- c(report, paste(i, ":", prediction, known))
    }
    if (prediction == "FALSE" & known){
      report <- c(report, paste(i, ":", prediction, known))
    }
  
  }
  return(report)
}