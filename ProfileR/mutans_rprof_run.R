library(gurobi)
library(sybil)
library(dplyr)
library(tictoc)

source('GitHub/PathwayMining/set_tools.R')
source('GitHub/PathwayMining/gurobi_tools.R')

Rprof()

model <- gurobi_read('GitHub/PathwayMining/data/mutans.lp')
coupled <- flux_coupling_raptor(model, directional_coupling=TRUE, partial_coupling=TRUE)$coupled

Rprof(NULL)
