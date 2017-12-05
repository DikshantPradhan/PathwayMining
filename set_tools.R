# functions related to identification and analysis of reaction sets

get_set_idx <- function(rxn, rxns_list){
  idx <- grep(core_rxn_id(rxn), rxns_list)
  #print(idx)
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

check_all_sets_for_containing <- function(comparison_sets, target_sets){
  for (set in comparison_sets){
    if (!check_sets_for_containing(set, target_sets)){
      print(set)
      return(FALSE)
    }
  }
  return(TRUE)
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

# return T/F based on whether or not two set-lists are equivalent
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

# find the sets combined in each set list formed from a rxn deletion
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
        error[get_rxn_idx(model$get_names()$VarName, set_lists[[i]][[j]])] <- paste(error[get_rxn_idx(model$get_names()$VarName, set_lists[[i]][[j]])], get_rxn_id_from_idx(model$get_names()$VarName, i), sep = " ")
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

# find the sets contained by 'rxns'
find_composing_sets <- function(rxns, sets){
  composition <- c()
  for (i in 1:length(sets)){
    if (length(sets[[i]]) == 0){next}
    if (all(sets[[i]] %in% rxns)){
      composition <- c(composition, i)
    }
  }

  return(composition)
}

# find the r0 sets combined to form the r1 sets
find_set_list_composition <- function(new_set_list, og_set_list){
  composition <- c()

  for (i in 1:length(new_set_list)){
    composing <- list(find_composing_sets(new_set_list[[i]], og_set_list))
    #print(composing)
    composition[i] <- composing
  }

  return(composition)
}

find_redundancies <- function(vars, og_set_list, composition_set){ # composition set is list of joined sets at each rxn deletion
  # redundancies <- matrix(nrow = 54, ncol = 54)
  # redundancies2 <- c()
  red1 <- c()
  red2 <- c()

  for (i in 1:length(composition_set)){ # deleted reaction
    # print(i)
    a <- get_set_idx(get_rxn_id_from_idx(vars, i), og_set_list) # containing set in og_set_list
    print(paste(a, "/", i, ":"))
    if (length(a) > 0 & length(composition_set[[i]]) > 0){ # check to make sure there are any changes due to this deletion
      for (j in 1:length(composition_set[[i]])){ # newly created sets
        for (k in 1:length(composition_set[[i]][[j]])){ # sets composing new sets
          set <- composition_set[[i]][[j]][k] # number
          rxns <- get_rxn_idx(vars, og_set_list[[set]][1]) # first reaction in each set
          if (rxns > length(composition_set)){next}
          # print(rxns)
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

remove_duplicate_pairs <- function(pairs){
  red1 <- pairs[,1]
  red2 <- pairs[,2]
  delete <- c()

  for (i in 1:(length(red1)-1)){
    for (j in (i+1):length(red1)){
      if (((red1[j] == red1[i]) & (red2[j] == red2[i])) | ((red1[j] == red2[i]) & (red2[j] == red1[i]))){
        delete <- c(delete, j)
      }

    }
  }

  if (length(delete) > 0){
	red1 <- red1[-delete]
  	red2 <- red2[-delete]
  }
  #red1 <- red1[-delete]
  #red2 <- red2[-delete]

  new_pairs <- cbind(red1, red2)
  # print(redundancies2)
  return(new_pairs)
}

optimize_suppression_idxs <- function(model, og_set_list){
  idxs <- c()
  for (i in og_set_list){
    idxs <- c(idxs, GRB_get_rxn_idx(model, i[[1]]))
  }

  return(idxs)
}

# return a T/F based on whether or not a pair exists in a pair-list
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

recurring_pairs <- function(target_pairs, comparison_pairs){
  elem_1 <- c()
  elem_2 <- c()

  if (length(target_pairs[,1]) < 1){
    return(cbind(elem_1, elem_2))
  }

  for (i in 1:length(target_pairs[,1])){
    if (check_for_pairs(target_pairs[i,], comparison_pairs)){
      elem_1 <- c(elem_1, target_pairs[i,1])
      elem_2 <- c(elem_2, target_pairs[i,2])
    }
  }

  return(cbind(elem_1, elem_2))
}

isolate_new_pairs_from_pair_lists <- function(og_pairs, pair_lists){
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

# extract the pairs in 'pairs' that are not in 'og pairs'
isolate_new_pairs <- function(og_pairs, pairs){
  new_rxn1 <- c()
  new_rxn2 <- c()

  for (i in 1:(length(pairs)/2)){
    if (!check_for_pairs(pairs[i,], og_pairs) & !check_for_pairs(pairs[i,], cbind(new_rxn1, new_rxn2))){
      new_rxn1 <- c(new_rxn1, pairs[i,1])
      new_rxn2 <- c(new_rxn2, pairs[i,2])
    }
  }

  return(cbind(new_rxn1, new_rxn2))
}

# extract the new pairs from the set lists not in the original set list
isolate_new_pairs_from_sets <- function(og_set_list, full_set_lists){
  new_rxn1 <- c()
  new_rxn2 <- c()

  og_pairs <- return_pairs_from_set_list(og_set_list)

  for (i in 1:length(full_set_lists)){
    print(i)
    if (length(full_set_lists[[i]]) == 0){next}
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

# extract only the new pairs formed from combining sets as specified in the composition
new_pairs_from_composition <- function(og_set_list, composition){
  new_rxn1 <- c()
  new_rxn2 <- c()

  for (i in 1:length(composition)){ # each deletion
    if (length(composition[[i]]) == 0){next}

    for (j in 1:length(composition[[i]])){ # each group of combined sets
      combined_sets <- composition[[i]][[j]]

      # isolate each pair in each composition
      for (k in 1:(length(combined_sets)-1)){
        set1_idx <- combined_sets[k]
        set1 <- og_set_list[[set1_idx]]
        for (l in (k+1):length(combined_sets)){
          set2_idx <- combined_sets[l]
          set2 <- og_set_list[[set2_idx]]

          # each reaction in each set
          for (m in set1){
            for (n in set2){
              if (!check_for_pairs(c(m,n), cbind(new_rxn1, new_rxn2))){
                new_rxn1 <- c(new_rxn1, m)
                new_rxn2 <- c(new_rxn2, n)
              }

            }
          }

        }
      }

    }
  }

  return(cbind(new_rxn1, new_rxn2))
}

return_pairs_from_set <- function(set){
  # print(set)
  #print(length(set))
  rxn1 <- c()
  rxn2 <- c()

  if (length(set) == 0){
    return(cbind(c(), c()))
  }

  if (length(set) == 1){
    return(cbind(c(set), c(set)))
  }

  # if (length(set) == 1){
  #   return(NULL)
  # }

  for (i in 1:(length(set))){
    for (j in (i):length(set)){
      #print(c(i,j))
      # print(c(set[i], set[j]))
      if (!is.na(set[i]) & !is.na(set[j]) & !is.null(set[i]) & !is.null(set[j])){ # & (set[i] != set[j])
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
    # print(set_list[[i]])
    # print(set_list[i])
    pairs <- return_pairs_from_set(set_list[[i]])
    # print(pairs)
    if (length(pairs) == 0){next}
    rxn1 <- c(rxn1, pairs[,1])
    rxn2 <- c(rxn2, pairs[,2])
    # if (pairs[1] != pairs[2]){
    #   rxn1 <- c(rxn1, pairs[,1])
    #   rxn2 <- c(rxn2, pairs[,2])
    # }

  }

  return(cbind(rxn1, rxn2))
}

# get pairs from multiple set lists
return_pair_lists_from_set_lists <- function(set_lists){
  pair_lists <- c()

  for (i in 1:length(set_lists)){
    pair_lists[[i]] <- return_pairs_from_set_list(set_lists[[i]])
  }

  return(pair_lists)
}

# combine two pair lists
append_pair_lists <- function(pairs1, pairs2){
  p1 <- c(pairs1[,1], pairs2[,1])
  p2 <- c(pairs1[,2], pairs2[,2])

  return(cbind(p1,p2))
}

# extract pairs based on a list of which sets have been combined
get_rxn_pairs_from_set_pairs <- function(og_sets, set_pairs){
  rxn1 <- c()
  rxn2 <- c()

  for (i in 1:nrow(set_pairs)){

    new_set <- union(og_sets[[set_pairs[i,1]]], og_sets[[set_pairs[i,2]]])
    # print(i)
    # print(new_set)
    new_pairs <- return_pairs_from_set(new_set)

    # print(new_pairs)

    rxn1 <- c(rxn1, new_pairs[,1])
    rxn2 <- c(rxn2, new_pairs[,2])
  }

  pairs <- cbind(rxn1, rxn2)
  pairs <- remove_duplicate_pairs(pairs)

  return(pairs)
}

# find out which of the target pairs is included in the pairs
isolate_pairs <- function(targets, pairs){
  item1 <- c()
  item2 <- c()

  for (i in 1:nrow(pairs)){
    for (j in 1:length(targets)){
      if (pairs[i,1] == targets[j] | pairs[i,2] == targets[j]){
        item1 <- c(item1, pairs[i,1])
        item2 <- c(item2, pairs[i,2])
      }
    }
  }

  return(cbind(item1, item2))
}

# list of length(set_list) where each entry is length of set at that idx
get_size_list <- function(set_list){
  size_list <- matrix(data = c(0), nrow = length(set_list), ncol = 1)

  for (i in 1:length(set_list)){
    size_list[i] <- length(set_list[[i]])
  }

  return(size_list)
}

# number of sets of each length
get_size_distribution <- function(set_list){

  size_list <- get_size_list(set_list)

  size_hist <- matrix(data = c(0), nrow = max(size_list), ncol = 1)

  # entry at index i is number of sets of length = i
  for (i in 1:length(size_list)){
    idx <- size_list[i]
    size_hist[idx] <- size_hist[idx] + 1
  }

  return(size_hist)
}

get_composition_size_distribution <- function(set_list, composition_set){
  composition_size <- c()

  for (i in 1:length(composition_set)){
    #print(i)
    if (length(composition_set[[i]]) == 0){next}

    comp_temp <- c()
    for (j in composition_set[[i]]){
      comp_temp <- c(comp_temp, length(set_list[[j]]))
      #idx <- composition_set[[i]][[j]]
      #for (k in idx){
      #  comp_temp <- c(comp_temp, length(set_list[[k]]))
      #}
      #print(idx)
    }
    #print(comp_temp)
    composition_size[i] <- list(comp_temp)
  }

  return(composition_size)
}

sample_sets_to_distribution <- function(elements, distribution, replacement = FALSE){

  sample <- sample(elements, size = length(elements), replace = replacement)

  set_list <- c()
  set_idx <- 1
  for (i in 1:length(distribution)){
    if (distribution[i] == 0){next}
    for (j in 1:distribution[i]){
      # print(paste(i, j, set_idx))
      # print(sample[1:i])
      # print(length(sample))
      set_list[set_idx] <- list(sample[1:i])
      sample <- sample[-c(1:i)]
      set_idx <- set_idx + 1

    }
  }

  return(set_list)
}

classify_sets_by_size <- function(set_list){
  size_dist <- get_size_distribution(set_list)
  size_list <- get_size_list(set_list)

  dist <- c()

  # index of list is size of all sets contained at that index
  for (i in 1:length(size_dist)){
    dist[i] <- list(which(size_list == i)) #c(dist[[l]], list(i))
  }

  return(dist)
}

sample_sets_to_composition <- function(size_class, size_composition){

  # randomize order of sets classified by size
  for (i in 1:length(size_class)){
    if (length(size_class[[i]]) == 1){ # passing in a single entry to sample returns i from 1:x instead of i == x
      size_class[[i]] <- size_class[[i]]
    }
    else {
      size_class[[i]] <- sample(size_class[[i]], size = length(size_class[[i]]), replace = FALSE)
    }
  }

  # new sampled set lit
  sampled_set_list <- c()

  for (i in 1:length(size_composition)){
    comp <- size_composition[[i]]
    #print(comp)

    sampled_set <- c()
    for (j in comp){
      #print(j)
      # take 1st from size_class and delete
      sampled_set <- c(sampled_set, size_class[[j]][1])
      size_class[[j]] <- size_class[[j]][-c(1)]
    }

    sampled_set_list[i] <- list(sampled_set)
  }

  # size_class should be empty
  for (i in size_class){if (length(i) > 0){print('len error')}}

  return(sampled_set_list)
}

sample_multiple_sets_to_distribution <- function(n_samples, elements, distribution, replacement = FALSE){
  set_lists <- c()

  for (i in 1:n_samples){
    set_lists[[i]] <- sample_sets_to_distribution(elements, distribution, replacement)
  }

  return(set_lists)
}

sample_multiple_sets_to_composition <- function(n_samples, size_class, size_composition){
  set_lists <- c()

  for (i in 1:n_samples){
    set_lists[[i]] <- sample_sets_to_composition(size_class, size_composition)
  }

  return(set_lists)
}
