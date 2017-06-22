data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=5000
nPnts=10000
steps=20

opt <- fluxVar(model, percentage = 99)
model_fva <- opt@lp_obj
fva_min <- model_fva[1:95]
fva_max <- model_fva[96:190]

S <- model@S

sample <- sampler(model)
rescaled_sample <- rescale_sample(sample)

# sample <- sampler_lm_fitting(model, 37, binary = TRUE)
# coupling1 <- flux_coupling_specific(rescale_sample(sample,37), 37, binary = TRUE)

plotMaxDiff <- function(maxDiff, idx){
  hist(maxDiff[,idx], breaks = seq(0, max(maxDiff[,idx])+0.005, 0.0001))
}

hist(rescaled_sample_[,40], breaks = seq(-1, 1, 0.001))
