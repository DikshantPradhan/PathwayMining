library(readr)
# scripts for loading different models

# HOW TO HANDLE YEAST MODEL
# read txt to data frame
# use minval to convert to tsv for sybil (writeTSVmod)
# use sybil to read tsv (readTSVmod)


## NECESSARY PRESETS FOR USAGE:
# model
# S

## ECOLI MODEL

get_ecoli_model <- function(){
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
  return(model)
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

get_yeast_model <- function(){

  setwd("~/GitHub/PathwayMining/data/yeast_model")
  # yeast_model <- readTSVmod(reactList = "Y7_test_react.tsv", metList = "Y7_met.tsv")
  # library(readr)
  # Y7_react_names <- read_delim("~/GitHub/PathwayMining/data/yeast_model/Y7_react_names.csv", ";", escape_double = FALSE, trim_ws = TRUE)
  # yeast_model@react_name <- Y7_react_names$Description
  # yeast_model <- rmReact(model = yeast_model, react = 1606)
  # yeast_model <- rmReact(model = yeast_model, react = 1590)
  #
  # yeast_exch_rxns <- grep("exchange", yeast_model@react_name)
  # for (i in yeast_exch_rxns){
  #   yeast_model <- changeBounds(yeast_model, i, lb = -1000, ub = 1000)
  # }

  yeast_model <- readTSVmod(reactList = "Y4_05_noCompart_react.tsv", metList = "Y4_05_noCompart_met.tsv")
  # yeast_4_05_compound <- read_delim("~/Documents/yeast_model/yeast_4/yeast_4_05_compound.csv", "\t", escape_double = FALSE, trim_ws = TRUE)

  yeast_4_05_noCompartments_compound <- read_delim("~/GitHub/PathwayMining/data/yeast_model/yeast_4_05_noCompartments_compound.csv", "\t", escape_double = FALSE, trim_ws = TRUE)

  # yeast_4_05_noCompartments_compound <- read.delim("~/GitHub/PathwayMining/data/yeast_model/yeast_4_05_noCompartments_compound.csv", "\t")

  for (i in 1:length(yeast_model@met_name)){
    # print(c(i, yeast_model@met_id[i], yeast_4_05_noCompartments_compound$NAME[which(yeast_4_05_noCompartments_compound$ID == yeast_model@met_id[i])]))
    yeast_model@met_name[i] <- yeast_4_05_noCompartments_compound$NAME[which(yeast_4_05_noCompartments_compound$ID == yeast_model@met_id[i])]
    # print(yeast_4_05_compound$Name[which(yeast_4_05_compound$ID == yeast_model@met_id[i])])
    # print(which(yeast_4_05_compound$Name == yeast_model@met_id[i]))
  }
  
  exch <- findExchReact(yeast_model)
  for (i in rev(exch@react_pos)){
    yeast_model@lowbnd[i] <- -1000
    yeast_model@uppbnd[i] <- 1000
  }
  
  blocked <- get_blocked(yeast_model)
  
  yeast_model@lowbnd[blocked] <- -1000
  yeast_model@uppbnd[blocked] <- 1000
  
  yeast_model@lowbnd[223] <- 0
  
  #yeast_model <- rmReact(model = yeast_model, react = 1297)
  #yeast_model <- rmReact(model = yeast_model, react = 1295)
  #yeast_model <- rmReact(model = yeast_model, react = 1293)
  
  # library(readr)
  #yeast_model <- rmReact(model = yeast_model, react = 1606)
  #yeast_model <- rmReact(model = yeast_model, react = 1590)

  # Y4_react_names <- read_delim("~/GitHub/PathwayMining/data/yeast_model/Y4_react_names.tsv", ";", escape_double = FALSE, trim_ws = TRUE)
  # for (i in 1:length(yeast_model@react_id)){
  #   num <- as.numeric(strsplit(yeast_model@react_id[i], split = "_")[[1]][2])
  #   # print(num)
  #   yeast_model@react_name[i] <- Y4_react_names$Description[num]
  # }
  #
  # yeast_model <- rmReact(model = yeast_model, react = get_rxn_idx(yeast_model@react_name, "r_1812"))
  # yeast_model <- rmReact(model = yeast_model, react = get_rxn_idx(yeast_model@react_name, "r_1815")) # not sure if this one warrants removal

  # yeast_exch_rxns <- grep("exchange", yeast_model@react_name)
  # for (i in yeast_exch_rxns){
  #   yeast_model <- changeBounds(yeast_model, i, lb = -1000, ub = 1000)
  # }
  setwd("~/GitHub/PathwayMining/")
  return(yeast_model)
}


##

# for (i in 1:nrow(yeast_4_05_noCompartments_reaction)){
#   match <- which(Y4$ID == yeast_4_05_noCompartments_reaction$X.ID[i])
#   yeast_4_05_noCompartments_reaction$GPR[i] <- Y4$GPR[match]
#   yeast_4_05_noCompartments_reaction$LOWER.BOUND[i] <- Y4$LOWER.BOUND[match]
#   yeast_4_05_noCompartments_reaction$UPPER.BOUND[i] <- Y4$UPPER.BOUND[match]
#   yeast_4_05_noCompartments_reaction$OBJECTIVE[i] <- Y4$OBJECTIVE[match]
# }

get_blocked <- function(model){
  lb <- which(model@lowbnd > -1000)
  ub <- which(model@uppbnd < 1000)
  
  potential_blocked <- intersect(lb, ub)
  blocked <- c()
  for (i in potential_blocked){
    if (model@lowbnd[i] == 0 & model@uppbnd[i] == 0){
      blocked <- c(blocked, i)
    }
  }
  
  return(blocked)
}

# model <- yeast_model

## ALL

# S <- model@S
