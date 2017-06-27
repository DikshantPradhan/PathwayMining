library(rstack)

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

get_rxn_id_from_idx <- function(rxn_idx){
  return(model@react_id[rxn_idx])
}

get_dwnst_rxns <- function(rxn_idx){
  dwnst_mets <- which(S[,rxn_idx] > 0)
  
  if (model@react_rev[rxn_idx]){
    dwnst_mets <- c(dwnst_mets, which(S[,rxn_idx] < 0))
  }

  # print(c("mets: ", dwnst_mets))
  
  dwnst_rxns <- c()
  
  for (i in dwnst_mets){
    # print(paste("i: ", i))
    dwnst_rxns <- c(dwnst_rxns, which(S[i,] < 0))
    # print(dwnst_rxns)
  }
  
  return(unique(dwnst_rxns))
}

get_upst_rxns <- function(rxn_idx){
  upst_mets <- which(S[,rxn_idx] < 0)
  
  if (model@react_rev[rxn_idx]){
    upst_mets <- c(upst_mets, which(S[,rxn_idx] > 0))
  }
  
  upst_rxns <- c()
  
  for (i in upst_mets){
    upst_rxns <- c(upst_rxns, which(S[i,] > 0)) 
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
      new_rxns = get_dwnst_rxns(rxn)
    }
    else {
      new_rxns = get_upst_rxns(rxn)
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
      
      # print(c(paste(rxn1, " & ", rxn2, "(", length(overlap), ")"), overlap))
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

find_replacement_EX_rxns <- function(gained_pairs){
  replacements <- c()
  for (i in 1:length(gained_pairs)){
    repl <- c()
    for (j in gained_pairs[[i]]){
      # print(j)
      rxns <- strsplit(j, split = " & ")[[1]]
      # print(rxns)
      if (j != "__" & rxns[1] == rxns[2] & grepl("EX_", rxns[1])){
        repl <- c(repl, rxns[1])
      }
      #print(rxns)
    }
    replacements[i] <- list(repl)
  }
  return(replacements)
}

find_new_essential_rxns <- function(gained_pairs){
  essential <- c()
  for (i in 1:length(gained_pairs)){
    ess <- c()
    for (j in gained_pairs[[i]]){
      # print(j)
      rxns <- strsplit(j, split = " & ")[[1]]
      # print(rxns)
      if (j != "__" & rxns[1] == "Biomass_Ecoli_core_w_GAM"){
        ess <- c(ess, rxns[2])
      }
      if (j != "__" & rxns[2] == "Biomass_Ecoli_core_w_GAM"){
        ess <- c(ess, rxns[1])
      }
      #print(rxns)
    }
    essential[i] <- list(ess)
  }
  return(essential)
}

convert_pair_strings_to_vector <- function(pair_string_list){
  rxn1 <- c()
  rxn2 <- c()
  
  for (i in pair_string_list[[1]]){
    #print(i)
    if (i == "__"){
      rxn1 <- c(rxn1, NA)
      rxn2 <- c(rxn2, NA)
    }
    else{
      rxns <- strsplit(i, split = " & ")[[1]]
      rxn1 <- c(rxn1, rxns[1])
      rxn2 <- c(rxn2, rxns[2])
    }
  }
  
  pair_vec <- cbind(rxn1, rxn2)
  #names(pair_vec) <- c("rxn1", "rxn2")
  return(pair_vec)
}

core_rxn_id <- function(rxn_id){ # rxn id with parenthesis
  return(strsplit(rxn_id, split = "\\(")[[1]][1])
}

get_list_of_sets <- function(pairs){ #2d columns
  rxns <- unique(union(pairs[,1], pairs[,2]))
  rxns_list <- c()
  
  for (i in rxns){
    #print(i)
    rxns_list <- c(rxns_list, list(i))
  }
  
  # rxn_list <- c(rxn_list[-c(1, 13)], list(union(rxn_list[1], rxn_list[13])))
  # grep("MALt2_2", rxn_list)
  
  for (i in 1:nrow(pairs)){
    # print(paste("round ",i))
    idx1 <- grep(core_rxn_id(pairs[i,1]), rxns_list)
    idx2 <- grep(core_rxn_id(pairs[i,2]), rxns_list)
    
    # print(idx1)
    # print(idx2)
    
    for (j in idx1){
      # print(j)
      # print(rxns_list[j])
      if (pairs[i,1] %in% rxns_list[[j]]){
        idx1 <- c()
        idx1 <- j
        # print(paste("found: ", j))
      }
    }
    for (j in idx2){
      # print(j)
      # print(rxns_list[j])
      if (pairs[i,2] %in% rxns_list[[j]]){
        idx2 <- c()
        idx2 <- j
        # print(paste("found: ", j))
      }
    }
    
    # if (length(idx1) == 0){
    #   idx1 <- which(rxns_list == pairs[i,1])
    # }
    # if (length(idx2) == 0){
    #   idx2 <- which(rxns_list == pairs[i,2])
    # }
    # print(i)
    # print(pairs[i,])
    # print(idx1)
    # print(idx2)
    # print(rxns_list)
    rxns_list <- c(rxns_list[-c(idx1, idx2)], list(union(unlist(rxns_list[idx1]), unlist(rxns_list[idx2]))))
  }
  #print(rxns_list)
  return(rxns_list)
}