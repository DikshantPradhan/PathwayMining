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

get_yeast_no_compart_model <- function(){

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
    yeast_model@met_name[i] <- yeast_4_05_noCompartments_compound$NAME[which(yeast_4_05_noCompartments_compound$ID == yeast_model@met_id[i])]
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

get_yeast_compart_model <- function(){
  
  setwd("~/GitHub/PathwayMining/data/yeast_model")

  yeast_model <- readTSVmod(reactList = "Y4_05_open_react.tsv", metList = "Y4_05_met.tsv")
  yeast_4_05_compound <- read_delim("~/Documents/yeast_model/yeast_4/yeast_4_05_compound.csv", "\t") # , escape_double = FALSE, trim_ws = TRUE
  
  for (i in 1:length(yeast_model@met_name)){
    # print(c(i, yeast_model@met_id[i], yeast_4_05_compound$Name[which(yeast_4_05_compound$ID == yeast_model@met_id[i])]))
    yeast_model@met_name[i] <- yeast_4_05_compound$Name[which(yeast_4_05_compound$ID == yeast_model@met_id[i])]
    # print(yeast_4_05_compound$Name[which(yeast_4_05_compound$ID == yeast_model@met_id[i])])
    # print(which(yeast_4_05_compound$Name == yeast_model@met_id[i]))
  }

  for (i in 1:length(yeast_model@react_id)){
    yeast_model@lowbnd[i] <- -1000
    yeast_model@uppbnd[i] <- 1000
  }
  
  remove <- which(yeast_model@react_name %in% c("growth", "biomass production", "lipid production"))
  
  yeast_model <- rmReact(model = yeast_model, react = remove[3])
  yeast_model <- rmReact(model = yeast_model, react = remove[2])
  yeast_model <- rmReact(model = yeast_model, react = remove[1])
  
  setwd("~/GitHub/PathwayMining/")
  return(yeast_model)
}

get_yeast_open_model <- function(){
  
  yeast_model <- readTSVmod(reactList = "Y4open_reactions.csv", metList = "Y4open_metabolites.csv")
  
  remove <- which(yeast_model@react_name %in% c("growth", "biomass production", "lipid production"))
  
  yeast_model <- rmReact(model = yeast_model, react = remove[3])
  yeast_model <- rmReact(model = yeast_model, react = remove[2])
  yeast_model <- rmReact(model = yeast_model, react = remove[1])
  
  return(yeast_model)
}

## MUTANS MODEL

get_mutans_model <- function(){
  mutans_model <- readTSVmod(reactList = "mutans_model_test.csv", metList = "mutans_model_met.csv")
  
  # print('biomass')
  # print(mutans_model@S[which(mutans_model@S[,477] != 0), 477])
  
  mutans_model@met_name[401] <- "DAP-type peptidoglycan"
  mutans_model@met_name[422] <- "Lys-type peptidoglycan"
  
  # make sure biomass split is implemented
  exch_idxs <- which(mutans_model@S[,477] != 0) # biomass
  # biom_consumed <- which(mutans_model@S[,477] < 0)
  # biom_produced <- which(mutans_model@S[,477] > 0)
  # # print(paste(mutans_model@met_id[exch_idxs], mutans_model@met_name[exch_idxs]))
  # print('consumed:')
  # print(paste(mutans_model@met_id[biom_consumed], mutans_model@met_name[biom_consumed], mutans_model@S[biom_consumed, 477]))
  # print('produced:')
  # print(paste(mutans_model@met_id[biom_produced], mutans_model@met_name[biom_produced], mutans_model@S[biom_produced, 477]))
  
  exch <- findExchReact(mutans_model)
  add_exch <- c()
  for (i in exch_idxs){
    if ((mutans_model@met_id[i] %in% exch@met_id) | (paste(mutans_model@met_id[i], '[e]', sep = "") %in% exch@met_id)){
      #print(mutans_model@react_id[i] %in% exch@met_id)
      next
    }
    
    else {
      # print(generate_exch_rxn(mutans_model, i))
      # add_exch <- c(add_exch, i)
      mutans_model <- addReact(mutans_model, paste('new_exch', i, sep = "_"), 
                               met = mutans_model@met_id[i], Scoef = c(-1*mutans_model@S[i, 477]), reversible = FALSE,
                               lb = 0, ub = 1000)
    }
  }

  # for (i in add_exch){
  #   mutans_model <- addReact(mutans_model, paste('new_exch', i, sep = "_"), 
  #                            met = mutans_model@met_id[i], Scoef = c(-1*mutans_model@S[i, 477]), reversible = FALSE,
  #                            lb = 0, ub = 1000)
  # }
  
  # ex <- findExchReact(mutans_model)
  # for (idx in 100:144){
  #   if(ex[idx]@met_id %in% c('C00013', 'C00009', 'C00080', 'C00008', 'C00035')){
  #     mutans_model@lowbnd[ex[idx]@react_pos] <- 0
  #     mutans_model@uppbnd[ex[idx]@react_pos] <- 1000
  #   }
  #   #mutans_model@uppbnd[rxn@react_pos] <- 0
  #   else{
  #     mutans_model@uppbnd[ex[idx]@react_pos] <- 0
  #     mutans_model@lowbnd[ex[idx]@react_pos] <- -1000
  #   }
  #   
  #   mutans_model@react_rev <- FALSE
  #   mutans_model@S[which(mutans_model@met_id == ex[idx]@met_id),idx] <- abs(mutans_model@S[which(mutans_model@met_id == ex[idx]@met_id), 477])
  # }
    
  print('removing:')
  print(mutans_model@react_name[477])
  mutans_model <- rmReact(model = mutans_model, react = 477)
  print(mutans_model@react_name[616])
  mutans_model <- rmReact(model = mutans_model, react = 616)
  print(mutans_model@react_name[366])
  mutans_model <- rmReact(model = mutans_model, react = 366)
  print(mutans_model@react_name[364])
  mutans_model <- rmReact(model = mutans_model, react = 364)
  
  return(mutans_model)
}

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

generate_exch_rxn <- function(model, rxn_idx){
  rxn_id <- model@met_id[rxn_idx]
  rxn_name <- model@met_name[rxn_idx]
  paste("exc00000	", rxn_name," exchange	", " <=> ", rxn_id, "					1	0	0	0		Unknown		", sep = "")
}

# model <- yeast_model

## ALL

# S <- model@S
