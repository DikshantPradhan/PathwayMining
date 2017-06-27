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

compare_correlation_sets <- function(comparison_num){
  og_pairs <- return_couples(flux_coupling_cor(sample_og))
  og_rxn_set <- get_list_of_sets(og_pairs)
  
  g <- make_empty_graph()
  g <- g + vertices(model@react_id, color = "green")
  g <- rxn_set_edges(g, og_rxn_set, "grey")
  
  added_rxns <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$gained_pairs[comparison_num]))
  lost_rxns <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$lost_pairs[comparison_num]))
  
  g <- rxn_set_edges(g, added_rxns, "green")
  g <- rxn_set_edges(g, lost_rxns, "red")
  plot(g)
}

compare_correlation_sets(12)
