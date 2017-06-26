data(Ec_core);
model=Ec_core;

S <- model@S

get_rxn_idx <- function(rxn_id){
  return(which(model@react_id == rxn_id))
}

get_rxn_name_from_id <- function(rxn_id){
  return(model@react_name[which(model@react_id == rxn_id)])
}

get_rxn_name_from_idx <- function(rxn_idx){
  return(model@react_name[rxn_idx])
}

find_dwnst_rxns <- function(rxn_idx){
  dwnst_mets <- which(S[,rxn_idx] > 0)
  
  dwnst_rxns <- c()
  
  for (i in dwnst_mets){
    dwnst_rxns <- which(S[i,] < 0)
  }
  
  return(unique(dwnst_rxns))
}

find_upst_rxns <- function(rxn_idx){
  upst_mets <- which(S[,rxn_idx] < 0)
  
  upst_rxns <- c()
  
  for (i in upst_mets){
    upst_rxns <- which(S[i,] > 0)
  }
  
  return(unique(upst_rxns))
}

get_dwnst_paths <- function(rxn_idx){
  paths = get_paths(rxn_idx, downstream = TRUE)
  return(paths)
}

get_upst_paths <- function(rxn_idx){
  paths = get_paths(rxn_idx, downstream = FALSE)
  return(paths)
}

## get list of reactions leading to production of species
get_paths <- function(rxn_idx, downstream = TRUE){
  #list of reactions to return
  rxns <- c()
  #spcs <- c()
  s <- stack$new()
  
  s$push(rxn_idx)
  
  while(!s$is_empty()){
    rxn = s$pop()
    #rxns = c(rxns, rxn)
    
    new_rxns = c()
    
    # new reactants and species to add
    if (downstream == TRUE){
      new_rxns = find_dwnst_rxns(rxn)
    }
    else {
      new_rxns = find_upst_rxns(rxn)
    }
    #new_spcs = c()
    
    # get reactants of new reactions
    for(new_rxn in new_rxns){
      # check to see if reaction is already in list of reactions to return
      if (!(new_rxn %in% rxns)){
        #new_spcs = c(new_spcs, get_rxn_rcts(new_rxn))
        rxns <- c(rxns, new_rxn)
        s$push(new_rxn)
      }
    }
    
    # # add reactions to stack
    # for(new_spc in new_spcs){
    #   # check to see if reactants have already been checked
    #   if (!(new_spc %in% spcs)){
    #     spcs <- c(spcs, new_spc)
    #     s$push(new_spc)
    #   }
    # }
  }
  
  return(rxns) # reaction names
}

find_coupling_connections <- function(coupling_list){
  for (i in 1:nrow(d_coupling)){
    rxn1 <- get_rxn_idx(coupling_list[i,1])
    rxn2 <- get_rxn_idx(coupling_list[i,2])
    
    if (rxn1 != rxn2){
      overlap <- c()
      overlap <- intersect(get_upst_paths(rxn1), get_dwnst_paths(rxn2))
      if (length(overlap) == 0){
        overlap <- intersect(get_upst_paths(rxn2), get_dwnst_paths(rxn1))
      }
      if (length(overlap) == 0){
        overlap <- intersect(get_upst_paths(rxn2), get_upst_paths(rxn1))
      }
      if (length(overlap) == 0){
        overlap <- intersect(get_dwnst_paths(rxn2), get_dwnst_paths(rxn1))
      }
      
      print(c(paste(rxn1, " & ", rxn2, "(", length(overlap), ")"), overlap))
    }
  }
}

find_lost_pairs <- function(coupling_list){
  lost_pairs <- find_changed_pairs(coupling_list, gained = FALSE)
  return(lost_pairs)
}

find_gained_pairs <- function(coupling_list){
  gained_pairs <- find_changed_pairs(coupling_list, gained = TRUE)
  return(gained_pairs)
}

find_changed_pairs <- function(coupling_list, gained = TRUE){
  pairs <- c()
  for (i in 1:nrow(coupling_list)){
    if ((gained == TRUE & abs(as.numeric(coupling_list[i,4])) > 0.85) | (gained == FALSE & abs(as.numeric(coupling_list[i,4])) < 0.15)){
      pairs <- c(pairs, paste(coupling_list[i,1], coupling_list[i, 2], sep = " & "))
    }
  }
  
  if (length(pairs) == 0){
    return("__")
  }
  
  return(pairs)
}
