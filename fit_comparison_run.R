data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=1000
nPnts=10000
steps=20

# sample <- sampler_lm_fitting(model, 37, binary = TRUE)
# coupling1 <- flux_coupling_specific(rescale_sample(sample,37), 37, binary = TRUE)

sample_og <- sampler(model)
print(maxDiff_dist(sample_og))

# sample_suppr <- sampler(suppressed_model(model, 37))
# coupling_og <- flux_coupling(rescale_sample(sample_og), binary = FALSE)
# coupling_suppr <- flux_coupling(rescale_sample(sample_suppr), binary = FALSE)
# d_couples1 <- coupling_change(coupling_og, coupling_suppr)


#write.csv(coupling_og[1, ,], "coupling_og2.csv")
#write.csv(coupling_suppr[1, ,], "coupling_suppr2.csv")
#write.csv(coupling_suppr[1, ,] - coupling_og[1, ,], "coupling_diff2.csv")