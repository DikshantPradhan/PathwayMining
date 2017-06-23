data(Ec_core);
model=Ec_core;

S <- model@S

find_dwnst_rxns <- function(S, rxn_idx){
  dwnst_mets <- which(S[,rxn_idx] > 0)
  
  dwnst_rxns <- c()
  
  for (i in dwnst_mets){
    dwnst_rxns <- which(S[i,] < 0)
  }
  
  return(unique(dwnst_rxns))
}

find_upst_rxns <- function(S, rxn_idx){
  upst_mets <- which(S[,rxn_idx] < 0)
  
  upst_rxns <- c()
  
  for (i in upst_mets){
    upst_rxns <- which(S[i,] > 0)
  }
  
  return(unique(upst_rxns))
}

### NEED TO CHANGE THIS
## get list of reactions leading to production of species
get_prod_path <- function(S, rxn_idx){
  #list of reactions to return
  rxns <- c()
  #spcs <- c()
  s <- stack$new()
  
  s$push(rxn_idx)
  
  while(!s$is_empty()){
    rxn = s$pop()
    spcs = c(spcs, spc)
    
    # new reactants and species to add
    new_rxns = get_producing_rxns(spc)
    new_spcs = c()
    
    # get reactants of new reactions
    for(new_rxn in new_rxns){
      # check to see if reaction is already in list of reactions to return
      if (!(new_rxn %in% rxns)){
        new_spcs = c(new_spcs, get_rxn_rcts(new_rxn))
        rxns <- c(rxns, new_rxn)
      }
    }
    
    # add reactants to stack
    for(new_spc in new_spcs){
      # check to see if reactants have already been checked
      if (!(new_spc %in% spcs)){
        spcs <- c(spcs, new_spc)
        s$push(new_spc)
      }
    }
  }
  
  return(rxns) # reaction names
}