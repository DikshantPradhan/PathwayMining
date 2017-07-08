get_set_idx <- function(rxn, rxns_list){
  idx <- grep(core_rxn_id(rxn), rxns_list)
  
  for (j in idx){
    if (rxn %in% rxns_list[[j]]){
      # idx <- c()
      # idx <- j
      return(j)
    }
  }
  return(integer(0))
}

check_sets_for_containing <- function(rxns, set_list){
  
  for (i in 1:length(set_list)){
    if (all(rxns %in% set_list[[i]])){
      return(TRUE)
    }
  }
  
  return(FALSE)
}

check_set_list_for_containing <- function(rxns, set_lists){
  lists <- c()
  
  for (i in 1:length(set_lists)){
    if (check_sets_for_containing(rxns, set_lists[[i]])){
      lists <- c(lists, i)
    }
  }
  
  return(lists)
}

find_all_sets_for_rxn <- function(rxn_id, set_lists){
  sets <- c()

  for (i in 1:95){
    set <- set_lists[[i]][get_set_idx(rxn_id, set_lists[[i]])]

    # print(paste(i, ":"))
    # print(set)

    if (length(set) > 0){
      sets[i] <- set
    }
  }
  
  return(sets)
}

find_all_sets_for_rxns <- function(rxns, set_lists){
  sets <- c()
  
  for (i in check_set_list_for_containing(rxns, set_lists)){
    # print(get_set_idx(rxns[1], set_lists[[i]]))
    # for (j in which(check_sets_for_containing(rxns, set_lists[[i]]))){
    #   set <- set_lists[[i]][[j]]
    #   if (length(set) > 0){
    #     sets[i] <- set
    #   }
    # }
    # print(i)
    set <- list(set_lists[[i]][[get_set_idx(rxns[1], set_lists[[i]])]])
    if (length(set) > 0){
      sets[i] <- set
    }
  }
  
  return(sets)
}

total_union <- function(sets){
  u <- c()
  
  for (i in 1:(length(sets))){ # length sets - 1
    u <- union(u, sets[[i]])
  }
  
  return(u)
}

compare_r1_sets <- function(og_set_list, set_lists){
  for (i in 1:length(set_lists)){
    print(i)
    containment <- "TRUE"
    for (j in 1:length(og_set_list)){
      if (check_sets_for_containing(og_set_list[[j]], set_lists[[i]]) == FALSE){
        print(og_set_list[[j]])
        for (rxn in og_set_list[[j]]){
          if (check_sets_for_containing(rxn, set_lists[[i]]) == TRUE){
            print(paste(rxn, "exists"))
          }
        }
        containment <- "FALSE"
      }
    }
    # print(containment)
  }
}