get_set_idx <- function(rxn, rxns_list){
  idx <- grep(core_rxn_id(rxn), rxns_list)
  # print(idx)
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

find_recurring_rxns <- function(rxns, set_lists){
  all_sets <- find_all_sets_for_rxns(rxns, set_lists)
  sets <- c()
  for (i in 1:length(all_sets)){
    if (length(all_sets[[i]]) > 0){
      sets <- c(sets, list(all_sets[[i]]))
    }
  }
  
  redundant <- c()
  
  for (i in 1:length(sets)){
    diff <- setdiff(sets[[i]], rxns)
    if (length(diff) > 0){
      redundant <- c(redundant, list(setdiff(sets[[i]], rxns)))
    }
  }
  redundant <- redundant[which(duplicated(redundant) == FALSE)]
  return(redundant)
}

total_union <- function(sets){
  u <- c()
  
  for (i in 1:(length(sets))){ # length sets - 1
    u <- union(u, sets[[i]])
  }
  
  return(u)
}

compare_r1_sets <- function(og_set_list, set_lists){ # see which og_sets don't appear in each set_list
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

return_composition_sets <- function(og_set_list, set_lists){
  composition <- c()
  for (i in 1:length(set_lists)){
    print(paste(i, ":"))
    sets <- c()
    for (j in 1:length(set_lists[[i]])){
      composing <- find_composing_sets(set_lists[[i]][[j]], og_set_list)
      if (length(composing) > 1){
        print(paste(composing))
        sets <- c(sets, list(composing))
      }
    }
    if (length(sets) > 0){
      composition[i] <- list(sets)
    }
  }
  
  return(composition)
}

find_composing_sets <- function(rxns, sets){
  composition <- c()
  for (i in 1:length(sets)){
    if (sets[[i]] %in% rxns){
      composition <- c(composition, i)
    }
  }
  
  return(composition)
}

find_redundancies <- function(composition_set){ # composition set is list of joined sets at each rxn deletion
  # redundancies <- matrix(nrow = 54, ncol = 54)
  # redundancies2 <- c()
  red1 <- c()
  red2 <- c()
  
  
  for (i in 1:length(composition_set)){ # deleted reaction
    a <- get_set_idx(get_rxn_id_from_idx(i), og_set_list) # containing set in og_set_list
    # print(paste(a, "/", i, ":"))
    if (length(a) > 0 & length(composition_set[[i]]) > 0){ # check to make sure there are any changes due to this deletion
      print(paste(a, "/", i, ":"))
      for (j in 1:length(composition_set[[i]])){ # newly created sets
        for (k in 1:length(composition_set[[i]][[j]])){ # sets composing new sets
          set <- composition_set[[i]][[j]][k] # number
          rxns <- get_rxn_idx(og_set_list[[set]][1]) # first reaction in each set
          # print(rxns)
          # new_compositions <- unlist(composition_set[[rxns]]) # sets creating by deletion of deleted rxn/set
          # print(rxns)
          # print(new_compositions[1])
          print(paste(set, "--"))
          for (new_set in composition_set[[rxns]]){
            if (a %in% new_set){
              # print(paste(set, "-"))
              print(new_set)
              # redundancies[a, set] <- as.character(list(new_set))
              # print(paste(set, new_compositions, sep = "- "))
              # print(paste("lengths:", length(new_set), length(composition_set[[i]][[j]][k]), length(set), length(a)))
              # print(new_set[-which(new_set == a)])
              # print(composition_set[[i]][[j]][-which(composition_set[[i]][[j]] == set)])
              if (all(new_set[-which(new_set == a)] == composition_set[[i]][[j]][-which(composition_set[[i]][[j]] == set)])){
                print(paste("~redundant~", a, "&", set))
                # redundancies2 <- c(redundancies2, paste(a, "&", set))
                red1 <- c(red1, a)
                red2 <- c(red2, set)
              }
            }
          }
          
          # if (a %in% new_compositions){
          #   print(paste(set, "-"))
          #   print(new_compositions)
          #   # print(paste(set, new_compositions, sep = "- "))
          # }
          
          
        }
        print("-- end set")
      }
    }
  }
  
  # redundancies <- cbind(red1, red2)
  
  # remove duplicates
  delete <- c()
  print(paste("removing deuplicates; ", length(red1), length(red2)))
  for (i in 1:(length(red1)-1)){
    for (j in (i+1):length(red1)){
      # print(i)
      # print(j)
      # print(paste(red1[i], red2[i], red1[j], red2[j]))
      # print(((red1[j] == red1[i]) & (red2[j] == red2[i])) | ((red1[j] == red2[i]) & (red2[j] == red1[i])))
      if (((red1[j] == red1[i]) & (red2[j] == red2[i])) | ((red1[j] == red2[i]) & (red2[j] == red1[i]))){
        delete <- c(delete, j)
        # print(paste("~",j))
        # red1 <- red1[-j]
        # red2 <- red2[-j]
      }
      # if ((red1[j] == red2[i]) & (red2[j] == red1[i])){
      #   # delete <- j
      #   print(j)
      #   red1 <- red1[-j]
      #   red2 <- red2[-j]
      # }
      
    }
  }
  
  red1 <- red1[-delete]
  red2 <- red2[-delete]
  redundancies <- cbind(red1, red2)
  # print(redundancies2)
  return(redundancies)
}