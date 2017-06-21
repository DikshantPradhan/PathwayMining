data(Ec_core);
model=Ec_core;
#solver="cplexAPI"
solver="glpkAPI"
W=250
nPnts=1000
steps=10

# sample <- sampler_lm_fitting(model, 37, binary = TRUE)
# coupling1 <- flux_coupling_specific(rescale_sample(sample,37), 37, binary = TRUE)

# sample_og <- sampler(model)
# maxDiff1 <- maxDiff_dist(sample_og)
# 
# steps=50
# sample_og <- sampler(model)
# maxDiff2 <- maxDiff_dist(sample_og)
# 
# steps=500
# sample_og <- sampler(model)
# maxDiff3 <- maxDiff_dist(sample_og)
# 
# steps=1000
# sample_og <- sampler(model)
# maxDiff4 <- maxDiff_dist(sample_og)

W = 200
sample <- sampler(model)
diff_w1 <- maxDiff_dist(sample)
W = 500
sample <- sampler(model)
diff_w2 <- maxDiff_dist(sample)
W = 1000
sample <- sampler(model)
diff_w3 <- maxDiff_dist(sample)
W = 2000
sample <- sampler(model)
diff_w4 <- maxDiff_dist(sample)

# steps = 50
# 
# sample_og <- sampler(model)
# sample_suppr <- sampler(suppressed_model(model, 37))
# 
# coupling_og <- flux_coupling(rescale_sample(sample_og), binary = FALSE)
# coupling_suppr <- flux_coupling(rescale_sample(sample_suppr), binary = FALSE)
# 
# coupling_og_cor <- flux_coupling_cor(rescale_sample(sample_og))
# coupling_suppr_cor <- flux_coupling_cor(rescale_sample(sample_suppr))

write.csv(coupling_suppr - coupling_og, "coupling_diff.csv")
write.csv(coupling_suppr_cor - coupling_og_cor, "coupling_diff_cor.csv")


plotMaxDiff <- function(maxDiff, idx){
  hist(maxDiff[,idx], breaks = seq(0, max(maxDiff[,idx])+0.005, 0.0001))
}

idx = 55
par(mfrow=c(2,2))
plotMaxDiff(diff_w1, idx)
plotMaxDiff(diff_w2, idx)
plotMaxDiff(diff_w4, idx)
plotMaxDiff(diff_w6, idx)

par(mfrow=c(2,2))
plotMaxDiff(maxDiff1, idx)
plotMaxDiff(maxDiff2, idx)
plotMaxDiff(maxDiff3, idx)
plotMaxDiff(maxDiff4, idx)