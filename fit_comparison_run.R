data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=200
nPnts=5000
steps=5

opt <- fluxVar(model, percentage = 99)
model_fva <- opt@lp_obj
fva_min <- model_fva[1:95]
fva_max <- model_fva[96:190]

S <- model@S

sample_og <- sampler(model)
#rescaled_sample <- rescale_sample(sample)

sample_comparison <- function(rxn_idx){
  sample_suppr <- sampler(suppressed_model(model, rxn_idx)) # exchange reactions: 20 - 39
  d_coupling <-find_coupling_change(sample_og, sample_suppr)
  return(d_coupling)
}

media_cond <- data.frame("rxn_idx" <- numeric(20))
names(media_cond) <- c("rxn_idx")

for (i in 20:30){
  coupling <- sample_comparison(i)
  j <- i- 19
  
  media_cond$rxn_idx[j] <- i
  media_cond$rxn_name[j] <- get_rxn_name_from_idx(i)
  
  if ("Biomass_Ecoli_core_w_GAM" %in% coupling[,1] | "Biomass_Ecoli_core_w_GAM" %in% coupling[,2]){
    media_cond$essential[j] <- TRUE
  }
  else {
    media_cond$essential[j] <- FALSE
  }
  
  media_cond$gained_pairs[j] <- find_gained_pairs(coupling)
  media_cond$lost_pairs[j] <- find_lost_pairs(coupling)
}

# find_coupling_connections(d_coupling)

# sample <- sampler_lm_fitting(model, 37, binary = TRUE)
# coupling1 <- flux_coupling_specific(rescale_sample(sample,37), 37, binary = TRUE)

# plotMaxDiff <- function(maxDiff, idx){
#   hist(maxDiff[,idx], breaks = seq(0, max(maxDiff[,idx])+0.005, 0.0001))
# }
# 
# hist(rescaled_sample_[,40], breaks = seq(-1, 1, 0.001))
