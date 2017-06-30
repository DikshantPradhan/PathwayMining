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