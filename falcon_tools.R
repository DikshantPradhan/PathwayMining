library(pryr)
library(data.tree)
source('~/GitHub/PathwayMining/logic_tools.R')

generate_falcon_model <- function(model, gene_sets = c(), rxn_sets = c()){
  # gprRules <- model@gprRules
  # genes <- model@genes
  og_dim <- dim(model@S)
  og_react_id <- model@react_id # model$get_names()$VarName
  og_met_id <- model@met_id
  og_genes <- model@genes # may need to pass this in
  og_allGenes <- model@allGenes # may need to pass this in
  rxnGeneMat <- model@rxnGeneMat # may not need
  colnames(rxnGeneMat) <- model@allGenes # may not need
  rownames(rxnGeneMat) <- model@react_id # may not need

  gene_set_bool <- FALSE
  if ((length(gene_sets) > 0) & (length(rxn_sets) > 0)){
    gene_set_bool <- TRUE
  }

  marked_genes <- matrix(FALSE, nrow = length(og_allGenes), ncol = 1)
  marked_rxns <- matrix(FALSE, nrow = length(og_react_id), ncol = 1)
  rownames(marked_genes) <- model@allGenes
  rownames(marked_rxns) <- model@react_id

  # clean gene sets, if passed in
  if (gene_set_bool){
    for (i in 1:length(gene_sets)){
      set <- gene_sets[[i]]
      removal <- which(set == "")
      if (length(removal) > 0){
        gene_sets[[i]] <- set[-c(which(set == ""))]
      }
    }
  }

  # add all exchange reactions
  for (gene in og_allGenes){
    new_met <- paste('a', gene, sep = '_')
    model <- addExchReact(model, met = new_met, lb = -1000, ub = 1000) # may need to alter this for grb
  }

  # rxn_removal_ids <- c()

  genes_from_path <- function(path, rxn_idx){
    gene_idxs <- gprRule_to_idx(path)
    genes <- og_genes[rxn_idx]
    return(genes[[1]][gene_idxs])
  }

  # helper function for simple_add and or_add
  normal_add <- function(model, new_met_list, rxn_id, simple = FALSE, addExch = FALSE, identifier = NULL){

    rxn_idx <- which(model@react_id == rxn_id)
    exch <- findExchReact(model)

    # metabolites of existing reaction
    old_met_idxs <- which(model@S[(1:og_dim[1]), rxn_idx] != 0)
    old_met_list <- og_met_id[old_met_idxs]
    old_met_coeff <- model@S[old_met_idxs, rxn_idx]

    # add exchange reactions, if needed
    # SHOULDN'T HAVE TO ADD ANY NEW EXCH REACTIONS
    # if (addExch){
    #   for (met in new_met_list){
    #     if (met %in% exch@met_id){next}
    #     model <- addExchReact(model, met, -1000, 1000)
    #   }
    # }

    met_list <- c(unlist(old_met_list), unlist(new_met_list))
    met_list <- met_list[!is.na(met_list)]

    if (!simple){ # add reverse reaction if needed # change from 0 -> 1000 & -1000 -> 0
      lowbnd <- model@lowbnd[rxn_idx]
      uppbnd <- model@uppbnd[rxn_idx]
      model <- addReact(model, paste(rxn_id, identifier, 'fwd', sep = '_'), met = met_list,
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = 0, ub = uppbnd, reversible = FALSE)
      model <- addReact(model, paste(rxn_id, identifier, 'rev', sep = '_'), met = met_list,
                        Scoef = c(old_met_coeff, rep(1, length(new_met_list))), lb = lowbnd, ub = 0, reversible = FALSE)
      # rxn_removal_ids <- c(rxn_removal_ids, rxn_id)

      ## NEED CONSTRAINTS TO PREVENT MODEL FROM PUSHING FLUX THROUGH BOTH DIRECTIONS AT ONCE
    }
    else { # changed from [-1000 to 1000] to [lowbnd to uppbnd]
      lowbnd <- model@lowbnd[rxn_idx]
      uppbnd <- model@uppbnd[rxn_idx]
      model <- addReact(model, rxn_id, met = met_list,
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = lowbnd, ub = uppbnd, reversible = TRUE)
    }

    # rxn_removal_ids <- c(rxn_removal_ids, rxn_id)
    return(model)
  }

  # used for direct addition of gene/enzyme to reaction; not used when reaction_activity is needed (or case)
  simple_add <- function(model, new_met_list, rxn_id, simple = TRUE){
    # print(c(rxn_id, ":", new_met_list))
    model <- normal_add(model, new_met_list, rxn_id, simple = TRUE)
    return(model)
  }

  # used for multiple gene combinations; reaction_activity needed
  or_add <- function(model, path_list, rxn_id, simple = TRUE, split = TRUE){
    rxn_idx <- which(og_react_id == rxn_id)
    rxn_activity <- paste('a', rxn_id, sep = "_")

    # add conversion for path to react_activity
    identifier <- 1 # need to differentiate breakdown of rxn in name
    for (mets in path_list){ # CHECK THIS FUNCTION IF EVERYTHING BREAKS
      genes <- genes_from_path(mets, rxn_idx)
      new_mets <- paste('a_', genes, sep = '')
      met_list <- c(unlist(new_mets), rxn_activity) # genes first, then reaction_activity
      met_list <- met_list[!is.na(met_list)] # need to get rid of this
      coeff_list <- c(unlist(rep(-1, length(new_mets))), 1)
      if (split){ # genes are always consumed, reaction_activity is consumed or produced; reactions proceed forwards
        model <- addReact(model, paste(rxn_activity, 'fwd conversion', identifier, sep = ' '), met = met_list,
                          Scoef = c(unlist(rep(-1, length(new_mets))), 1), lb = 0, ub = 1000, reversible = FALSE)
        model <- addReact(model, paste(rxn_activity, 'rev conversion', identifier, sep = ' '), met = met_list,
                          Scoef = c(unlist(rep(-1, length(new_mets))), -1), lb = 0, ub = 1000, reversible = FALSE)
      }
      else {
        model <- addReact(model, paste(rxn_activity, 'conversion', identifier, sep = ' '), met = met_list,
                          Scoef = coeff_list, lb = -1000, ub = 1000, reversible = TRUE)
      }
      identifier <- identifier + 1
    }

    # add activity specific to react
    model <- normal_add(model, new_met_list = c(rxn_activity), rxn_id, simple = simple, addExch = FALSE)

    return(model)
  }

  # list of all rxns w which a gene participates
  gene_rxn_recurrence <- c()
  for (i in 1:length(model@allGenes)){
    gene_rxn_recurrence[i] <- list(which(rxnGeneMat[,i] == TRUE))
  }
  names(gene_rxn_recurrence) <- model@allGenes
  gene_rxn_promiscuity <- sapply(gene_rxn_recurrence, function(x) length(x))

  # loyal genes and reactions

  #print('RXN EXCLUSIVITY')

  rxn_exclusive_genes <- which(gene_rxn_promiscuity == 1)
  loyal_rxns <- c()
  for (i in gene_rxn_recurrence[rxn_exclusive_genes]){
    if (all(sapply(model@genes[i], function(x) x %in% names(rxn_exclusive_genes)))){
      loyal_rxns <- c(loyal_rxns, i)
    }
  }
  loyal_rxn_idxs <- unique(loyal_rxns)

  #print(paste('loyal rxns:', length(loyal_rxn_idxs)))
  #print(paste('rxn exclusive genes:', length(rxn_exclusive_genes)))

  # all reactions which are comprised of non-promiscuous genes
  for (i in loyal_rxn_idxs){
    rxn_id <- og_react_id[i]
    gpr_rule <- model@gprRules[i]
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (length(gpr_paths) == 1){
      genes <- genes_from_path(gpr_paths[[1]], i)
      model <- simple_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      marked_rxns[i] <- TRUE
    }
    else {
      genes <- c()
      for (path in gpr_paths){
        genes <- c(genes, genes_from_path(path, i))
      }
      genes <- unlist(genes)
      if (genes[1] == ''){next}
      model <- or_add(model, gpr_paths, og_react_id[i], split = TRUE)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      marked_rxns[i] <- TRUE
    }
  }

  remaining_genes <- model@allGenes[-c(which(marked_genes == TRUE))]
  remaining_rxns <- og_react_id[-c(which(marked_rxns == TRUE))]

  #print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  #print(paste('remaining genes:', length(which(marked_genes == FALSE))))

  # IF NO GENE SET INFORMATION, THEN SKIP THIS
  if (gene_set_bool){
  # list of gene sets in which gene participates (idxs same as rxn sets)
  gene_set_recurrence <- c()
  for (i in 1:length(model@allGenes)){
    recur <- c()
    for (j in 1:length(gene_sets)){
      if (model@allGenes[i] %in% gene_sets[[j]]){
        recur <- c(recur, j)
      }
    }
    gene_set_recurrence[i] <- list(recur)
  }
  names(gene_set_recurrence) <- model@allGenes

  gene_set_promiscuity <- sapply(gene_set_recurrence, function(x) length(x))
  set_exclusive_genes <- which(gene_set_promiscuity == 1)
  loyal_sets <- c()
  for (i in gene_set_recurrence[set_exclusive_genes]){
    # if (marked_genes[i]){next}
    if (all(sapply(gene_sets[[i]], function(x) x %in% names(set_exclusive_genes)))){
      loyal_sets <- c(loyal_sets, i)
    }
  }
  loyal_set_idxs <- unique(loyal_sets)

  # genes exclusive to a single set, but not a single gene

  #print('SET EXCLUSIVITY')

  #print(paste('loyal sets:', length(loyal_set_idxs)))
  #print(paste('set exclusive genes:', length(set_exclusive_genes)))

  for (i in loyal_set_idxs){
    for (j in rxn_sets[[i]]){
      rxn_idx <- which(model@react_id == j)
      gpr_rule <- model@gprRules[rxn_idx]
      gpr_paths <- find_gpr_paths(gpr_rule)
      genes <- c()
      for (path in gpr_paths){
        genes <- c(genes, genes_from_path(path, rxn_idx))
      }
      genes <- unlist(genes)
      if (which(og_react_id == j) %in% which(marked_rxns)){
        next
      }
      if (nchar(gpr_paths[1]) == 0){
        next
      }
      model <- or_add(model, gpr_paths, j)
      marked_rxns[rxn_idx] <- TRUE
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
    }
  }

  #print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  #print(paste('remaining genes:', length(which(marked_genes == FALSE))))

  }

  # NORMAL ADD FOR ALL REMAINING REACTIONS

  for (rxn_idx in which(marked_rxns == FALSE)){
    # print(rxn_idx)
    rxn_id <- og_react_id[rxn_idx]
    # print(rxn_id)
    gpr_rule <- model@gprRules[rxn_idx]
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (nchar(gpr_paths[1]) == 0){next}

    #print('testing mets')
    #print(which(model@S[(1:og_dim[1]),rxn_idx] != 0))
    model <- or_add(model, gpr_paths, rxn_id, split = TRUE)

    # id <- 1
    # for (path in gpr_paths){
    #   genes <- genes_from_path(path, rxn_idx)
    #   # print(genes)
    #   # print(rxn_id)
    #   model <- normal_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id, identifier = id)
    #   rxn_removal_ids <- c(rxn_removal_ids, rxn_id)
    #   marked_genes[which(model@allGenes %in% genes)] <- TRUE
    #   id <- id + 1
    # }
    marked_rxns[rxn_idx] <- TRUE
  }

  # if gene_activity is only consumed, then exchange reaction should only flow inwards??

  #print('FINAL COUNT')
  #print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  #print(paste('remaining genes:', length(which(marked_genes == FALSE))))

  # print(marked_genes)

  # print('RXN REMOVAL')
  # # print(rxn_removal_ids)
  # print(paste('met num:', length(model@met_id)))
  # for (rxn_id in unique(rxn_removal_ids)){
  #   # print(rxn_id)
  #   model <- rmReact(model, react = rxn_id, rm_met = FALSE)
  # }
  # print(paste('met num:', length(model@met_id)))

  return(model)
}

gprRule_to_idx <- function(list){
  for (i in 1:length(list)){
    list[i] <- gsub("[^0-9]", "", list[i])
  }

  return(as.numeric(list))
}

clean_rxn_names_in_set <- function(set_list){

  clean_ex_a <- function(name){
    new_name <- strsplit(name, 'Ex_a_')[[1]]
    return(new_name[2])
  }
  
  new_sets <- c()

  for (i in 1:length(set_list)){
    k <- 1
    new_set <- c()
    for (j in 1:length(set_list[[i]])){
      new_name <- clean_ex_a(set_list[[i]][j])
      if (!is.na(new_name)){
        # set_list[[i]][j] <- clean_ex_a(set_list[[i]][j])
        new_set[k] <- clean_ex_a(set_list[[i]][j])
        k <- k+1
      }
      else {
        # set_list[[i]][j] <- ""
        next
      }
    }
    if (length(new_set) > 0){
      new_sets[i] <- list(new_set)
    }
  }

  return(new_sets)
}

isolate_gene_matrix <- function(coupling_matrix){
  row_genes <- which(grepl('Ex_a', rownames(coupling_matrix)))
  col_genes <- which(grepl('Ex_a', colnames(coupling_matrix)))

  gene_matrix <- coupling_matrix[row_genes, col_genes]

  return(gene_matrix)
}

isolate_rxn_matrix <- function(coupling_matrix){
  row_rxns <- which(!grepl('a_', rownames(coupling_matrix)))
  col_rxns <- which(!grepl('a_', colnames(coupling_matrix)))

  rxn_matrix <- coupling_matrix[row_rxns, col_rxns]

  return(rxn_matrix)
}

# for (i in 1:length(mutans_model@react_id)){
#   if (mutans_model@lowbnd[i] != mutans_falcon@lowbnd[i]){print(c(i, "wrong lowbnd"))}
#   if (mutans_model@uppbnd[i] != mutans_falcon@uppbnd[i]){print(c(i, "wrong uppbnd"))}
# }

yeast_falcon_compar_react <- function(react_id){
  idx <- which(yeast_model@react_id == react_id)
  f_idx <- which(yeast_falcon@react_id == react_id)

  met_idx <- which(yeast_model@S[,idx] != 0)
  f_met_idx <- which(yeast_falcon@S[,f_idx] != 0)

  met_coeff <- yeast_model@S[met_idx, idx]
  f_met_coeff <- yeast_falcon@S[f_met_idx, f_idx]

  met_id <- yeast_model@met_id[met_idx]
  f_met_id <- yeast_falcon@met_id[f_met_idx]

  print(react_id)
  print('Yeast Model')
  print(paste(met_coeff, met_id))
  print(paste(yeast_model@lowbnd[idx], yeast_model@uppbnd[idx]))
  print('Yeast Falcon')
  print(paste(f_met_coeff, f_met_id))
  print(paste(yeast_falcon@lowbnd[f_idx], yeast_falcon@uppbnd[f_idx]))
}

yeast_falcon_compar_met <- function(met_id){
  met_idx <- which(yeast_model@met_id == met_id)
  f_met_idx <- which(yeast_falcon@met_id == met_id)

  r_idxs <- which(yeast_model@S[met_idx,] != 0)
  f_r_idxs <- which(yeast_falcon@S[f_met_idx,] != 0)

  print(met_id)
  print('Yeast Model')
  print(yeast_model@react_id[r_idxs])
  print('Yeast Falcon')
  print(yeast_falcon@react_id[f_r_idxs])
}

# yeast_falcon_compar_met("s_0097")
# yeast_falcon_compar_met("s_0102")
# yeast_falcon_compar_react("r_0036")
# yeast_falcon_compar_react("r_0084")
# yeast_falcon_compar_react("r_0360")
# yeast_falcon_compar_react("r_1414")
#yeast_falcon_compar_react("EX_s_0199")

#model@react_id[3484]
#biomass_mets <- which(model@S[,3484] != 0)
#for (met in biomass_mets){
#  print(paste(model@met_id[met], length(which(model@S[met,] != 0))))
#}
