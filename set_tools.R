# functions related to identification and analysis of reaction sets

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

check_sets_for_deletion <- function(rxn_list, sets){
  deleted <- c()
  for (i in rxn_list){
    if (!check_sets_for_containing(i, sets)){
      deleted <- c(deleted, i)
    }
  }

  return(deleted)
}

check_set_list_for_deletion <- function(rxn_list, set_list){
  deletions <- c()
  for (i in 1:length(set_list)){
    deletions[[i]] <- check_sets_for_deletion(rxn_list, set_list[[i]])
  }

  return(deletions)
}

find_all_sets_for_rxn <- function(rxn_id, set_lists){
  sets <- c()

  for (i in 1:length(set_lists)){
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

compare_sets <- function(set_1, set_2){
  if (length(set_1) != length(set_2)){
    return(FALSE)
  }
  
  for (i in 1:length(set_1)){
    if (!all.equal(set_1[[i]], set_2[[i]])){
      return(FALSE)
    }
  }
  
  return(TRUE)
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

return_composition_sets <- function(og_set_list, set_lists, model){
  composition <- c()

  error <- c() # 94 is not a set number

  for (i in 1:length(model$get_names()$VarName)){
    error[i] <- paste(model$get_names()$VarName[i], ":", sep = "")
  }

  for (i in 1:length(set_lists)){
    print(paste(i, ":"))
    sets <- c()
    for (j in 1:length(set_lists[[i]])){
      composing <- find_composing_sets(set_lists[[i]][[j]], og_set_list)
      if (length(composing) > 1){
        print(paste(composing))
        sets <- c(sets, list(composing))
      }
      if (length(composing) == 0){
        print(c("error", i, j, set_lists[[i]][[j]]))
        error[get_rxn_idx(model$get_names()$VarName, set_lists[[i]][[j]])] <- paste(error[get_rxn_idx(model$get_names()$VarName, set_lists[[i]][[j]])], get_rxn_id_from_idx(model, i), sep = " ")
      }
    }
    if (length(sets) > 0){
      composition[i] <- list(sets)
    }
  }

  #rownames(error) <- ecoli$get_names()$VarName
  #print(error)

  composition_sets <- c()

  composition_sets$composition <- composition
  composition_sets$error <- error

  return(composition_sets)
}

find_composing_sets <- function(rxns, sets){
  composition <- c()
  for (i in 1:length(sets)){
    if (all(sets[[i]] %in% rxns)){
      composition <- c(composition, i)
    }
  }

  return(composition)
}

find_redundancies <- function(og_set_list, composition_set){ # composition set is list of joined sets at each rxn deletion
  # redundancies <- matrix(nrow = 54, ncol = 54)
  # redundancies2 <- c()
  red1 <- c()
  red2 <- c()

  for (i in 1:length(composition_set)){ # deleted reaction
    a <- get_set_idx(get_rxn_id_from_idx(i), og_set_list) # containing set in og_set_list
    print(paste(a, "/", i, ":"))
    if (length(a) > 0 & length(composition_set[[i]]) > 0){ # check to make sure there are any changes due to this deletion
      for (j in 1:length(composition_set[[i]])){ # newly created sets
        for (k in 1:length(composition_set[[i]][[j]])){ # sets composing new sets
          set <- composition_set[[i]][[j]][k] # number
          rxns <- get_rxn_idx(og_set_list[[set]][1]) # first reaction in each set
          print(paste(set, "--"))
          for (new_set in composition_set[[rxns]]){
            if (a %in% new_set){
              print(new_set)
              if (all(new_set[-which(new_set == a)] == composition_set[[i]][[j]][-which(composition_set[[i]][[j]] == set)])){
                print(paste("~redundant~", a, "&", set))
                red1 <- c(red1, a)
                red2 <- c(red2, set)
              }
            }
          }
        }
        print("-- end set")
      }
    }
  }

  # remove duplicates
  delete <- c()
  print(paste("removing deuplicates; ", length(red1), length(red2)))
  for (i in 1:(length(red1)-1)){
    for (j in (i+1):length(red1)){
      if (((red1[j] == red1[i]) & (red2[j] == red2[i])) | ((red1[j] == red2[i]) & (red2[j] == red1[i]))){
        delete <- c(delete, j)
      }

    }
  }

  red1 <- red1[-delete]
  red2 <- red2[-delete]
  redundancies <- cbind(red1, red2)
  # print(redundancies2)
  return(redundancies)
}

optimize_suppression_idxs <- function(model, og_set_list){
  idxs <- c()
  for (i in og_set_list){
    idxs <- c(idxs, GRB_get_rxn_idx(model, i[[1]]))
  }

  return(idxs)
}

check_for_pairs <- function(pair, pair_list){
  #print(pair[1])
  #print(length(pair_list[,1]))
  if (length(pair_list[,1]) < 1){
    return(FALSE)
  }
  for (i in 1:length(pair_list[,1])){
    #print("~")
    #print(pair)
    #print(pair_list[i,])
    if ((pair[1] == pair_list[i, 1] & pair[2] == pair_list[i, 2]) | (pair[1] == pair_list[i, 2] & pair[2] == pair_list[i, 1])){
      return(TRUE)
    }
  }
  return(FALSE)
}

isolate_new_pairs <- function(og_pairs, pair_lists){
  new_rxn1 <- c()
  new_rxn2 <- c()

  for (i in 1:length(pair_lists)){
    for (j in 1:(length(pair_lists[[i]])/2)){
      #print(pair_lists[[i]][j,])
      if ((pair_lists[[i]][j,1] != pair_lists[[i]][j,2]) & !check_for_pairs(pair_lists[[i]][j,], og_pairs) & !check_for_pairs(pair_lists[[i]][j,], cbind(new_rxn1, new_rxn2))){
        new_rxn1 <- c(new_rxn1, pair_lists[[i]][j,1])
        new_rxn2 <- c(new_rxn2, pair_lists[[i]][j,2])
      }
    }
  }

  return(cbind(new_rxn1, new_rxn2))
}

isolate_new_pairs_from_sets <- function(og_set_list, full_set_lists){
  new_rxn1 <- c()
  new_rxn2 <- c()

  og_pairs <- return_pairs_from_set_list(og_set_list)

  for (i in 1:length(full_set_lists)){
    pairs <- return_pairs_from_set_list(full_set_lists[[i]])

    for (j in 1:length(pairs[,1])){
      if(!check_for_pairs(pairs[j,], og_pairs) & !check_for_pairs(pairs[j,], cbind(new_rxn1, new_rxn2))){
        new_rxn1 <- c(new_rxn1, pairs[j,1])
        new_rxn2 <- c(new_rxn2, pairs[j,2])
      }
    }
  }

  return(cbind(new_rxn1, new_rxn2))
}

return_pairs_from_set <- function(set){
  #print(length(set))
  rxn1 <- c()
  rxn2 <- c()
  if (length(set) == 1){
    return(cbind(c(set), c(set)))
  }
  for (i in 1:(length(set)-1)){
    for (j in (i+1):length(set)){
      #print(c(i,j))
      #print(c(set[i], set[j]))
      if (!is.na(set[i]) & !is.na(set[j])){
        #print(c(set[i], set[j]))
        rxn1 <- c(rxn1, set[i])
        rxn2 <- c(rxn2, set[j])
      }
    }
  }

  return(cbind(rxn1, rxn2))
}

return_pairs_from_set_list <- function(set_list){

  rxn1 <- c()
  rxn2 <- c()

  for (i in 1:length(set_list)){
    pairs <- return_pairs_from_set(set_list[[i]])
    rxn1 <- c(rxn1, pairs[,1])
    rxn2 <- c(rxn2, pairs[,2])
  }

  return(cbind(rxn1, rxn2))
}

return_pair_lists_from_set_lists <- function(set_lists){
  pair_lists <- c()

  for (i in 1:length(set_lists)){
    pair_lists[[i]] <- return_pairs_from_set(set_lists[[i]])
  }

  return(pair_lists)
}
