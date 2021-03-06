# functions related to model access and usage

library(rstack)
library(datastructures)

# data(Ec_core);
# model=Ec_core;

get_rxn_idx <- function(vars, rxn_id){
  return(which(vars == rxn_id))
}

get_rxn_name_from_id <- function(vars, rxn_id){
  return(model@react_name[which(vars == rxn_id)])
}

get_rxn_name_from_idx <- function(model, rxn_idx){
  return(model@react_name[rxn_idx])
}

get_rxn_id_from_idx <- function(vars, rxn_idx){
  return(vars[rxn_idx])
}

get_path_mtx_between_reactions <- function(S, start, end,
    active_mets = matrix(data = TRUE, nrow = nrow(S), ncol = 1)){ # start and end are rows (metabolites) of S

  q = queue()
  q <- insert(q, as.character(start))

  tracing_mtx <- matrix(data = FALSE, nrow = nrow(S), ncol = nrow(S)) # rows are parents, columns are children

  traverse <- function(q, index, active_mets, tracing_mtx){
    # print(index)
    linked_rxns <- which(S[index,] != 0)
    linked_mets <- unique(unlist(lapply(linked_rxns, function(x) which(S[,x] != 0))))
    new_linked_mets <- intersect(which(active_mets), linked_mets)
    active_mets[new_linked_mets] <- FALSE
    tracing_mtx[index, new_linked_mets] <- TRUE
    # print(which(active_mets))
    for (met in new_linked_mets){
      # print(met)
      q <- insert(q, as.character(met))
    }

    list(q = q, active_mets = active_mets, tracing_mtx = tracing_mtx)
  }

  while (size(q) > 0 & active_mets[end]){
    next_met <- as.numeric(pop(q))
    # if (next_met == end){break}
    traverse_iter <- traverse(q, next_met, active_mets, tracing_mtx)

    q <- traverse_iter$q
    active_mets <- traverse_iter$active_mets
    tracing_mtx <- traverse_iter$tracing_mtx
  }

  return(tracing_mtx) # path matrix; TRUE indicates that row idx is parent of col idx
}

trace_path_mtx_between_reactions <- function(mtx, start, end){ # start is leaf of tree encoded in matrix
  # path <- c(end)
  path <- c()
  # print(paste(start, end))
  # print(start != end)
  while (start != end){
    start <- which(mtx[,start])
    # print(start)
    path <- c(start, path)
  }

  return(path)
}

get_dwnst_rxns <- function(rxn_idx, model, sample = NULL){
  S <- model@S
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

get_upst_rxns <- function(rxn_idx, model){
  S <- model@S

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

print_mets_and_coeffs <- function(model, rxn_idx){
  print(model@react_id[rxn_idx])
  met_idxs <- which(model@S[,rxn_idx] != 0)
  print(paste(test_model@met_id[met_idxs], test_model@S[met_idxs, rxn_idx]))
}

check_for_duplicate_reactions <- function(model){
  S <- model@S
  n_rxn <- dim(S)[2]

  for (i in 1:(n_rxn-1)){
    test <- S[,i]
    for (j in (i+1):n_rxn){
      test_2 <- S[,j]
      if (identical(test_2, test) | identical(test_2, -1*test)){
        print(paste('error', i, j))
        # print(test_2)
        # print(test)
      }
    }
  }
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

group_sets <- function(set_list, rxn1, rxn2){
  idx1 <- get_set_idx(rxn1, set_list)
  idx2 <- get_set_idx(rxn2, set_list)

  set_list <- c(set_list[-c(idx1, idx2)], list(union(unlist(set_list[idx1]), unlist(set_list[idx2]))))

  return(set_list)
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

get_union_set_from_degen_pairs <- function(vars, pair_lists){
  # set_list <- vars

  for (i in 1:length(pair_lists)){
    set_list <- get_list_of_sets(pair_lists[[i]])
    # print(set_list)
    print(paste(i, ": ", length(set_list)))
  }

  return(set_list)
}

replace_names_in_set <- function(set, og_names, new_names){
  for (i in 1:length(set)){
    for (j in 1:length(set[[i]])){
      idx <- which(og_names == set[[i]][[j]])
      set[[i]][[j]] <- new_names[idx]
    }
  }

  return(set)
}

pair_data <- function(og_sets, set_lists){
  all_rxn_pairs <- return_pairs_from_set(unlist(og_sets))
  r0_rxn_pairs <- return_pairs_from_set_list(og_sets)
  new_r1_pairs <- isolate_new_pairs_from_sets(og_sets, set_lists)
  r1_pairs <- append_pair_lists(r0_rxn_pairs, new_r1_pairs)

  all_gene_pairs <- get_gene_pairs_from_rxn_pair_list(all_rxn_pairs)
  r0_gene_pairs <- get_gene_pairs_from_rxn_pair_list(r0_rxn_pairs)
  r1_gene_pairs <- get_gene_pairs_from_rxn_pair_list(r1_rxn_pairs)
  new_r1_gene_pairs <- get_gene_pairs_from_rxn_pair_list(new_r1_rxn_pairs)

  pair_data <- c()
  pair_data$all_rxn_pairs <- all_rxn_pairs
  pair_data$r0_rxn_pairs <- r0_rxn_pairs
  pair_data$r1_rxn_pairs <- r1_rxn_pairs
  pair_data$new_r1_rxn_pairs <- new_r1_rxn_pairs

  pair_data$all_gene_pairs <- all_gene_pairs
  pair_data$r0_gene_pairs <- r0_gene_pairs
  pair_data$r1_gene_pairs <- r1_gene_pairs
  pair_data$new_r1_gene_pairs <- new_r1_gene_pairs

  return(pair_data)
}

get_react_coeffs <- function(model, idx){
  met_idx <- which(model@S[,idx] != 0)
  met_coeff <- model@S[met_idx, idx]
  met_id <- model@met_id[met_idx]
  print(paste(met_coeff, met_id))
  print(paste(model@lowbnd[idx], model@uppbnd[idx]))
}
