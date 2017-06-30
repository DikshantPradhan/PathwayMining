data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=200
nPnts=500
steps=5

# opt <- fluxVar(model, percentage = 99)
# model_fva <- opt@lp_obj
# fva_min <- model_fva[1:95]
# fva_max <- model_fva[96:190]

S <- model@S

# sample_og <- sampler(model)
#rescaled_sample <- rescale_sample(sample)

sample_suppr_generator <- function(){
  sample_suppr <- array(0, dim = c(20, nPnts, 95))
  
  for (i in 20:39){
    sample <- sampler(suppressed_model(model, i)) # exchange reactions: 20 - 39
    for (j in 1:nPnts){
      for (k in 1:95){
        sample_suppr[i - 19, j, k] <- sample[j,k]
      }
    }
  }
  
  return(sample_suppr)
}

sample_suppr_full <- array(0, dim = c(20, nPnts, 95))
sample_comparison <- function(rxn_idx){
  sample <- sampler(suppressed_model(model, rxn_idx)) # exchange reactions: 20 - 39
  # print(dim(sample_compar))
  d_coupling <- find_coupling_change(sample_og, sample)
  
  for (j in 1:nPnts){
    for (k in 1:95){
      sample_suppr_full[i - 19, j, k] <- sample[j,k]
    }
  }
  
  return(d_coupling)
}

media_cond_generator <- function(){
  media_cond <- data.frame("rxn_idx" <- numeric(20))
  names(media_cond) <- c("rxn_idx")
  
  for (i in 20:39){
    coupling <- sample_comparison(i)
    j <- i - 19
    
    media_cond$rxn_idx[j] <- i
    media_cond$rxn_name[j] <- get_rxn_name_from_idx(i)
    
    if ("Biomass_Ecoli_core_w_GAM" %in% coupling[,1] | "Biomass_Ecoli_core_w_GAM" %in% coupling[,2]){
      media_cond$essential[j] <- TRUE
    }
    else {
      media_cond$essential[j] <- FALSE
    }
    
    media_cond$gained_pairs[j] <- list(find_gained_pairs(coupling))
    media_cond$lost_pairs[j] <- list(find_lost_pairs(coupling))
  }
  
  return(media_cond)
}

# sample_suppr <- sample_suppr_generator()
# media_cond <- media_cond_generator()

# compare_correlation_sets(12)
