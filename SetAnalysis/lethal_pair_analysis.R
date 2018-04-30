library(readr)
library(Matrix)
SGD_synthetic_lethal_pairs <- read_table2("~/GitHub/PathwayMining/data/SGD_synthetic_lethal_pairs.csv")

source('~/GitHub/PathwayMining/set_tools.R')
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g1_sets.RData")
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g0_sets.RData")


identify_captured_pairs <- function(pairs, sets){ # 0: not capture, -1: one gene in pair is not in set, 1: captured
  captured <- Matrix(data = 0, nrow = nrow(pairs), ncol = 1, sparse = TRUE)
  for (i in 1:nrow(pairs)){
    idx1 <- get_set_idx(pairs[i,1], sets)
    idx2 <- get_set_idx(pairs[i,2], sets)
    if (length(idx1) < 1 | length(idx2) < 1){
      captured[i] <- -1
      next
    }
    if (idx1 == idx2){captured[i] <- 1}
  }
  return(captured)
}

pairs <- SGD_synthetic_lethal_pairs
new_pairs <- matrix(nrow = nrow(pairs), ncol = ncol(pairs))
for (i in 1:nrow(pairs)){
  for (j in 1:ncol(pairs)){
    new_pairs[i,j] <- paste('Ex_a_', pairs[i,j], sep = '')
  }
}

# yeast_g1_captured_lethal_pairs <- identify_captured_pairs(new_pairs, yeast_g1_sets) # 56 captured, 906 not captured, 51526 not in sets
yeast_g0_captured_lethal_pairs  <- identify_captured_pairs(new_pairs, g0_sets) # 26 captured, 936 not captured, 31526 not in sets
