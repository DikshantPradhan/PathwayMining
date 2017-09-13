# scripts for loading different models

# HOW TO HANDLE YEAST MODEL
# read txt to data frame
# use minval to convert to tsv for sybil (writeTSVmod)
# use sybil to read tsv (readTSVmod)


## NECESSARY PRESETS FOR USAGE:
# model
# S

## ECOLI MODEL

data(Ec_core);
model=Ec_core;
model <- changeBounds(model, 11, lb = 0) # this and next lead to no PGI in GND coset
# model <- changeBounds(model, 13, lb = 0, ub = 0)
model <- rmReact(model = model, react = 13)
for (i in findExchReact(model)@react_pos){
  model <- changeBounds(model, i, lb = -1000, ub = 1000)
  # if (model@lowbnd[i] == 0){
  #   model <- changeBounds(model, i, lb = -1000)
  # }
}

#solver="cplexAPI"
# solver="glpkAPI"
# W=8000
# warmup=8000
# nPnts=3000
# steps=2

# opt <- fluxVar(model, percentage = 99)
# model_fva <- opt@lp_obj
# fva_min <- model_fva[1:95]
# fva_max <- model_fva[96:190]

## YEAST MODEL

setwd("~/GitHub/PathwayMining/data/yeast_model")
yeast_model <- readTSVmod(reactList = "Y7_test_react.tsv", metList = "Y7_met.tsv")
yeast_model <- rmReact(model = yeast_model, react = 1606)
yeast_model <- rmReact(model = yeast_model, react = 1590)

setwd("~/GitHub/PathwayMining/")

# model <- yeast_model

## ALL

# S <- model@S
