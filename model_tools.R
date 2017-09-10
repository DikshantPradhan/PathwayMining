# functions related to model access and usage

library(rstack)

# data(Ec_core);
# model=Ec_core;

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

get_dwnst_rxns <- function(rxn_idx, sample = NULL){
  # fwd_dwnst_mets <- which(S[,rxn_idx] > 0)

  # if (length(model@react_rev[rxn_idx]) == 0){
  #   print(rxn_idx)
  # }
  # if (model@react_rev[rxn_idx]){ # add reactants to downstream metabolites if rxn is reversible
  #   dwnst_mets <- c(dwnst_mets, which(S[,rxn_idx] < 0))
  # }

  dwnst_rxns <- c()
  for (i in which(S[,rxn_idx] > 0)){ # considering forward flux: products

    if (!is.null(sample)){
      if (length(which(sample[, rxn_idx] > 0)) > 0){ # rxn has fwd flux

        for (j in which(S[i,] < 0)){ # reactions w downstream metabolite as reactant
          if (length(which(sample[, j] > 0)) != 0){ # rxn has fwd flux & ds_rxn has fwd flux
            # print(paste(get_rxn_id_from_idx(rxn_idx), "fwd,", get_rxn_id_from_idx(j), "fwd:", model@met_id[i]))
            dwnst_rxns <- c(dwnst_rxns, j)
          }
        }

        for (j in which(S[i,] > 0)){ # reactions w downstream metabolite as product
          if (length(which(sample[, j] < 0)) != 0){ # rxn has fwd flux & ds_rxn has rev flux
            # print(paste(get_rxn_id_from_idx(rxn_idx), "fwd,", get_rxn_id_from_idx(j), "rev:", model@met_id[i]))
            dwnst_rxns <- c(dwnst_rxns, j)
          }
        }

      }
    }
    else {
      # print("null sample")
      for (j in which(S[i,] < 0)){ # reactions w downstream metabolite as reactant
        dwnst_rxns <- c(dwnst_rxns, j)
      }

      for (j in which(S[i,] > 0)){ # reactions w downstream metabolite as product
        if (model@react_rev[j]){ # downstream rxn is reversible
          dwnst_rxns <- c(dwnst_rxns, j)
        }
      }
      # dwnst_rxns <- c(dwnst_rxns, j)
    }
  }

  for (i in which(S[,rxn_idx] < 0)){ # considering reverse flux: reactants

    if (!is.null(sample)){
      if (length(which(sample[, rxn_idx] < 0)) > 0){ # rxn has rev flux

        for (j in which(S[i,] < 0)){ # reactions w upstream metabolite as reactant
          if (length(which(sample[, j] > 0)) != 0){ # rxn has rev flux & ds_rxn has fwd flux
            # print(paste(get_rxn_id_from_idx(rxn_idx), "rev,", get_rxn_id_from_idx(j), "fwd:", model@met_id[i]))
            dwnst_rxns <- c(dwnst_rxns, j)
          }
        }

        for (j in which(S[i,] > 0)){ # reactions w upstream metabolite as product
          if (length(which(sample[, j] < 0)) != 0){ # rxn has rev flux & ds_rxn has rev flux
            # print(paste(get_rxn_id_from_idx(rxn_idx), "rev,", get_rxn_id_from_idx(j), "rev:", model@met_id[i]))
            dwnst_rxns <- c(dwnst_rxns, j)
          }
        }

      }
    } else {
      # print("null reverse")
      if (model@react_rev[rxn_idx]){ # rxn has rev flux
        # print("null sample")
        for (j in which(S[i,] < 0)){ # reactions w downstream metabolite as reactant
          dwnst_rxns <- c(dwnst_rxns, j)
        }

        for (j in which(S[i,] > 0)){ # reactions w downstream metabolite as product
          if (model@react_rev[j]){ # downstream rxn is reversible
            dwnst_rxns <- c(dwnst_rxns, j)
          }
        }
        # dwnst_rxns <- c(dwnst_rxns, j)

      }
    }
  }

  return(unique(dwnst_rxns))
}

get_dwnst_rxns_2 <- function(rxn_idx, sample = NULL){
  fwd_dwnst_mets <- which(S[,rxn_idx] > 0)

  # if (length(model@react_rev[rxn_idx]) == 0){
  #   print(rxn_idx)
  # }
  rev_dwnst_mets <- c()
  if (model@react_rev[rxn_idx]){
    rev_dwnst_mets <- which(S[,rxn_idx] < 0)
  }

  dwnst_mets <- c()

  # if sample data is given, isolate products if flux is positive and reactants if flux is negative; both if flux goes both ways
  if (!is.null(sample)){
    flux <- sample[, rxn_idx]
    if (length(which(flux < 0)) == 0){ # all positive
      # print(paste("fwd:", get_rxn_id_from_idx(rxn_idx)))
      dwnst_mets <- fwd_dwnst_mets
    }
    if (length(which(flux > 0)) == 0){ # all negative
      # print(paste("rev:", get_rxn_id_from_idx(rxn_idx)))
      dwnst_mets <- rev_dwnst_mets
    }
  }
  else {
    # print(paste("both ways:", get_rxn_id_from_idx(rxn_idx)))
    dwnst_mets <- c(fwd_dwnst_mets, rev_dwnst_mets)
  }

  dwnst_rxns <- c()
  # dwnst_search <- c()
  fwd_rxns <- c()
  rev_rxns <- c()
  for (i in dwnst_mets){
    # dwnst_rxns <- c(dwnst_rxns, which(S[i,] < 0))
    for (j in which(S[i,] < 0)){ # reactions w downstream metabolite as reactant
      if (!is.null(sample)){
        flux <- sample[, j]
        if (length(which(sample[, j] < 0)) == 0){ # all positive
          dwnst_rxns <- c(dwnst_rxns, i)
          print(paste(get_rxn_id_from_idx(rxn_idx), "(fwd):", get_rxn_id_from_idx(i)))
        }
      }
      else {
        print(paste(get_rxn_id_from_idx(rxn_idx), "(fwd_):", get_rxn_id_from_idx(i)))
        dwnst_rxns <- c(dwnst_rxns, fwd_rxns)
      }
    }

    fwd_rxns <- c(fwd_rxns, which(S[i,] < 0)) # reactions w downstream metabolite as reactant
    rev_rxns <- c(rev_rxns, which(S[i,] > 0)) # reactions w downstream metabolite as product (helpful if reversible)
    # dwnst_search <- c(dwnst_search, paste(model@met_id[i], ": ", get_rxn_id_from_idx(which(S[i,] < 0)), sep = ""))
  }

  for (i in fwd_rxns){
    # if sample data is given,
    if (!is.null(sample)){
      flux <- sample[, i]
      if (length(which(flux < 0)) == 0){ # all positive
        dwnst_rxns <- c(dwnst_rxns, i)
        print(paste(get_rxn_id_from_idx(rxn_idx), "(fwd):", get_rxn_id_from_idx(i)))
      }
    }
    else {
      print(paste(get_rxn_id_from_idx(rxn_idx), "(fwd_):", get_rxn_id_from_idx(i)))
      dwnst_rxns <- c(dwnst_rxns, fwd_rxns)
    }
  }
  for (i in rev_rxns){
    if (!is.null(sample)){
      flux <- sample[, i]
      if (length(which(flux > 0)) == 0){ # all negative
        print(paste(get_rxn_id_from_idx(rxn_idx), "(rev):", get_rxn_id_from_idx(i)))
        dwnst_rxns <- c(dwnst_rxns, i)
      }
    }
    else {
      if (model@react_rev[i]){
        print(paste(get_rxn_id_from_idx(rxn_idx), "(rev_):", get_rxn_id_from_idx(i)))
        dwnst_rxns <- c(dwnst_rxns, i)
      }
    }
  }

  # print(paste(get_rxn_id_from_idx(rxn_idx), ":", get_rxn_id_from_idx(dwnst_rxns)))
  # print(dwnst_search)
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
  s <- stack$new()

  s$push(rxn_idx)

  while(!s$is_empty()){
    rxn = s$pop()
    new_rxns = c()

    # new reactants and species to add
    if (downstream == TRUE){
      new_rxns = get_dwnst_rxns(rxn)
    }
    else {
      new_rxns = get_upst_rxns(rxn)
    }

    # get reactants of new reactions
    for(new_rxn in new_rxns){
      # check to see if reaction is already in list of reactions to return
      if (!(new_rxn %in% rxns)){
        #new_spcs = c(new_spcs, get_rxn_rcts(new_rxn))
        rxns <- c(rxns, new_rxn)
        s$push(new_rxn)
      }
    }
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

convert_pair_vector_to_string <- function(coupling_list){
  pairs <- c()
  for (i in 1:nrow(coupling_list)){
    pairs <- c(pairs, paste(coupling_list[i,1], coupling_list[i, 2], sep = " & "))
  }
  return(pairs)
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

group_sets <- function(set_list, rxn1, rxn2){
  idx1 <- get_set_idx(rxn1, set_list)
  idx2 <- get_set_idx(rxn2, set_list)

  set_list <- c(set_list[-c(idx1, idx2)], list(union(unlist(set_list[idx1]), unlist(set_list[idx2]))))

  return(set_list)
}

get_list_of_sets <- function(pairs, rxns_list = c()){ #2d columns

  if (length(rxns_list) == 0){
    print("new rxn list")
    rxns <- unique(union(pairs[,1], pairs[,2]))
    # rxns_list <- c()

    for (i in rxns){
      rxns_list <- c(rxns_list, list(i))
    }
  }

  # rxn_list <- c(rxn_list[-c(1, 13)], list(union(rxn_list[1], rxn_list[13])))
  # grep("MALt2_2", rxn_list)

  for (i in 1:nrow(pairs)){
    idx1 <- get_set_idx(pairs[i,1], rxns_list) #grep(core_rxn_id(pairs[i,1]), rxns_list)
    idx2 <- get_set_idx(pairs[i,2], rxns_list) #grep(core_rxn_id(pairs[i,2]), rxns_list)

    # print(paste(pairs[i, 1], "&", pairs[i,2], ":", idx1, idx2))

    rxns_list <- c(rxns_list[-c(idx1, idx2)], list(union(unlist(rxns_list[idx1]), unlist(rxns_list[idx2]))))
  }

  return(rxns_list)
}

correlating_sets_from_sample <- function(sample){
  pairs <- return_couples(flux_coupling_cor(sample))
  rxn_set <- get_list_of_sets(pairs)
  return(rxn_set)
}

generate_pair_lists <- function(suppression_idxs){
  pair_lists <- c()
  for (i in suppression_idxs){
    pair_lists[i] <- list(return_couples(flux_coupling_cor(sampler(suppressed_model(model, i)))))
  }
  return(pair_lists)
}

generate_set_lists <- function(suppression_idxs){
  set_lists <- c()
  for (i in suppression_idxs){
    set_lists[i] <- list(correlating_sets_from_sample(sampler(suppressed_model(model, i))))
  }
  return(set_lists)
}

generate_og_set_list <- function(){
  return(correlating_sets_from_sample(sampler(model)))
}

generate_og_pair_list <- function(){
  return(return_couples(flux_coupling_cor(sampler(suppressed_model(model, i)))))
}

get_union_set_from_degen_pairs <- function(set_list = model@react_id, pair_lists){

  for (i in 1:length(pair_lists)){
    set_list <- get_list_of_sets(pair_lists[[i]], rxns_list = set_list)
    # print(set_list)
    print(paste(i, ": ", length(set_list)))
  }

  return(set_list)
}
