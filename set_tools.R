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

correlation_predictor <- function(sample_df){ # need to complete this function
  #
  cor_list <- rep(FALSE, 95)
  cor_list[check_set_list_for_containing(c("CYTBD", "NADH16"), set_lists)] <- TRUE
  
  correlation_list <- c()
  
  sample_df = sampler(model, nPnts = 500)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
  for (i in 1:nrow(sample_df)){
    if (check_sets_for_containing(c("CYTBD", "NADH16"), og_set_list)){
      correlation_list <- c(correlation_list, TRUE)
    }
    else {
      correlation_list <- c(correlation_list, FALSE)
    }
  }
  
  
  report <- paste(i, "::", "sample_df:", nrow(sample_df), ";", "corr:", length(correlation_list))
  
  for (i in 1:95){
    sample_df2 = sampler(suppressed_model(model, i), nPnts = 500)
    sample_df <- rbind.data.frame(sample_df, sample_df2)
    for (j in 1:(nrow(sample_df) - length(correlation_list))){
      correlation_list <- c(correlation_list, cor_list[i])
    }
    report <- c(report, paste(i, "::", "sample_df:", nrow(sample_df), ";", "corr:", length(correlation_list)))
  }
  sample_df <- data.frame(sample_df)
  fit <- tree(correlation_list ~ ., data = sample_df)
  
  plot(fit)
  text(fit)
  return(fit)
}

predictor_matrix_generator <- function(){ # need to complete this function
  
  sample_df <- c() # fill this in
  predictors <- c()
  
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

cor_list <- rep(FALSE, 95)
cor_list[check_set_list_for_containing(c("CYTBD", "NADH16"), set_lists)] <- TRUE

correlation_list <- c()

sample_df = sampler(model, nPnts = 500)# ACHR(model,W,nPoints=nPnts,stepsPerPoint=10)
for (i in 1:nrow(sample_df)){
  if (check_sets_for_containing(c("CYTBD", "NADH16"), og_set_list)){
    correlation_list <- c(correlation_list, TRUE)
  }
  else {
    correlation_list <- c(correlation_list, FALSE)
  }
}

sample_df["block"] <- rep(0, nrow(sample_df))

report <- paste(i, "::", "sample_df:", nrow(sample_df), ";", "corr:", length(correlation_list))

for (i in 1:95){
  sample_df2 = sampler(suppressed_model(model, i), nPnts = 500)
  sample_df2["block"] <- rep(i, nrow(sample_df2))
  sample_df <- rbind.data.frame(sample_df, sample_df2)
  for (j in 1:(nrow(sample_df) - length(correlation_list))){
    correlation_list <- c(correlation_list, cor_list[i])
  }
  report <- c(report, paste(i, "::", "sample_df:", nrow(sample_df), ";", "corr:", length(correlation_list)))
}
sample_df <- data.frame(sample_df)
fit1 <- tree(correlation_list ~ ., data = sample_df)

plot(fit1)
text(fit1)