library(pryr)
library(data.tree)

generate_falcon_model <- function(model, gene_sets, rxn_sets){
  # gprRules <- model@gprRules
  # genes <- model@genes
  og_dim <- dim(model@S)
  og_react_id <- model@react_id
  og_genes <- model@genes
  rxnGeneMat <- model@rxnGeneMat
  colnames(rxnGeneMat) <- model@allGenes
  rownames(rxnGeneMat) <- model@react_id
  
  marked_genes <- matrix(FALSE, nrow = length(model@allGenes), ncol = 1)
  marked_rxns <- matrix(FALSE, nrow = length(model@react_id), ncol = 1)
  rownames(marked_genes) <- model@allGenes
  rownames(marked_rxns) <- model@react_id
  
  # clean gene sets
  for (i in 1:length(gene_sets)){
    set <- gene_sets[[i]]
    removal <- which(set == "")
    if (length(removal) > 0){
      gene_sets[[i]] <- set[-c(which(set == ""))]
    }
  }
  
  # add all exchange reactions
  for (gene in model@allGenes){
    new_met <- paste('a', gene, sep = '_')
    model <- addExchReact(model, met = new_met, lb = -1000, ub = 1000)
  }
  
  rxn_removal_idxs <- c()
  
  genes_from_path <- function(path, rxn_idx){
    gene_idxs <- gprRule_to_idx(path)
    genes <- og_genes[rxn_idx]
    return(genes[[1]][gene_idxs])
  } 
  
  normal_add <- function(model, new_met_list, rxn_id, simple = FALSE, addExch = FALSE){
    rxn_idx <- which(og_react_id == rxn_id)
    exch <- findExchReact(model)
    
    # metabolites of existing reaction
    old_met_idxs <- which(model@S[, rxn_idx] != 0)
    old_met_list <- model@met_id[old_met_idxs]
    old_met_coeff <- model@S[old_met_idxs, rxn_idx]
    
    # add exchange reactions, if needed
    if (addExch){
      for (met in new_met_list){
        if (met %in% exch@met_id){next}
        model <- addExchReact(model, met, -1000, 1000)
      }
    }
    
    # check to see if rxn already exists
    new_react <- matrix(0, nrow = dim(model@S)[1], ncol = 1)
    rownames(new_react) <- model@met_id
    new_react <- model@S[, rxn_idx]
    new_react[new_met_list] <- -1 # indexing error here
      
    for (i in 1:dim(model@S)[2]){
      if (identical(new_react, model@S[, i])){
        return(model)
      }
    }
    
    
    if (!simple){ # add reverse reaction if needed
      model <- addReact(model, rxn_id, met = c(old_met_list, new_met_list), 
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = 0, ub = 1000)
      model <- addReact(model, paste(rxn_id, 'rev', sep = '_'), met = c(old_met_list, new_met_list), 
                        Scoef = c(old_met_coeff, rep(1, length(new_met_list))), lb = -1000, ub = 0)
      ## need constraints to prevent model from pushing flux through both directions at once
    }
    else {
      model <- addReact(model, rxn_id, met = c(old_met_list, new_met_list), 
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = -1000, ub = 1000)
    }
    
    rxn_removal_idxs <- c(rxn_removal_idxs, rxn_idx)
    return(model)
  }
  
  simple_add <- function(model, new_met_list, rxn_id, simple = FALSE){
    model <- normal_add(model, new_met_list, rxn_id, simple = TRUE)
    return(model)
  }
  
  or_add <- function(model, path_list, rxn_id, simple = TRUE){
    rxn_idx <- which(og_react_id == rxn_id)
    rxn_activity <- paste('a', rxn_id, sep = "_")
  
    # add conversion for path to react_activity
    for (mets in path_list){
      genes <- genes_from_path(mets, rxn_idx)
      new_mets <- paste('a_', genes, sep = '')
      met_list <- c(unlist(new_mets), rxn_activity)
      coeff_list <- c(unlist(rep(-1, length(new_mets))), 1)
      model <- addReact(model, paste(rxn_activity, 'conversion', sep = ' '), met = met_list, 
                        Scoef = coeff_list, lb = -1000, ub = 1000)
    }
    
    # add activity specific to react
    model <- normal_add(model, new_met_list = c(rxn_activity), rxn_id, simple = simple, addExch = FALSE)
    
    ## NEED TO ADD FLUX CONSTRAINT TO THIS
    
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
  
  print('RXN EXCLUSIVITY')
  
  rxn_exclusive_genes <- which(gene_rxn_promiscuity == 1)
  loyal_rxns <- c()
  for (i in gene_rxn_recurrence[rxn_exclusive_genes]){
    if (all(sapply(model@genes[i], function(x) x %in% names(rxn_exclusive_genes)))){
      loyal_rxns <- c(loyal_rxns, i)
    }
  }
  loyal_rxn_idxs <- unique(loyal_rxns)
  
  print(paste('loyal rxns:', length(loyal_rxn_idxs)))
  print(paste('rxn exclusive genes:', length(rxn_exclusive_genes)))
  
  # all reactions which are comprised of non-promiscuous genes
  for (i in loyal_rxn_idxs){
    gpr_rule <- model@gprRules[i]
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (length(gpr_paths) == 1){
      genes <- genes_from_path(gpr_paths[[1]], i)
      model <- simple_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      marked_rxns[i] <- TRUE
    }
    else {
      model <- or_add(model, gpr_paths, og_react_id[i])
      # for (path in gpr_paths){
      #   genes <- genes_from_path(path, i)
      #   # print(genes)
      #   model <- simple_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      marked_rxns[i] <- TRUE
      # }
    }
  }
  
  remaining_genes <- model@allGenes[-c(which(marked_genes == TRUE))]
  remaining_rxns <- og_react_id[-c(which(marked_rxns == TRUE))]
  
  print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  print(paste('remaining genes:', length(which(marked_genes == FALSE))))
  
  # IF NO GENE SET INFORMATION, THEN SKIP THIS
  
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
  
  print('SET EXCLUSIVITY')
  
  for (i in loyal_set_idxs){
    # print(paste('new set:', length(rxn_sets[[i]])))
    # print(gene_sets[[i]])
    for (j in rxn_sets[[i]]){
      rxn_idx <- which(og_react_id == j)
      
      if (which(og_react_id == j) %in% which(marked_rxns)){
        # print(paste('loyal rxn;', model@gpr[rxn_idx]))
        next
      }
      print('set excl')
      gpr_rule <- model@gprRules[rxn_idx]
      gpr_paths <- find_gpr_paths(gpr_rule)
      # genes <- genes_from_path(gpr_paths[[1]], rxn_idx)
      
      # print(model@gpr[rxn_idx])
      model <- or_add(model, gpr_paths, og_react_id[i])
    }
    
    marked_genes[gene_sets[[i]]] <- TRUE
    marked_rxns[rxn_sets[[i]]] <- TRUE
  }
  
  print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  print(paste('remaining genes:', length(which(marked_genes == FALSE))))
  
  # # genes exclusive to a single set, but not a single gene
  # for (i in which(gene_set_promiscuity == 1)){
  #   if (marked_genes[i]){next}
  #   gene <- names(gene_set_promiscuity)[i]
  #   k <- gene_set_recurrence[[i]]
  #   set <- rxn_sets[[k]]
  #   # identify if rxn is loyal to genes exclusive to this set: if yes then simple_add and mark gene
  #   # if gene is already marked, remove from paths (watch for newly empty paths)
  #   for (j in set){
  #     if (marked_rxns[which(og_react_id == j)]){next}
  #     # print(model@genes[which(og_react_id == j)])
  #   }
  # }
  
  # NORMAL ADD FOR ALL REMAINING REACTIONS
  
  for (rxn_idx in which(marked_rxns == FALSE)){
    rxn_id <- og_react_id[rxn_idx]
    gpr_rule <- model@gprRules[rxn_idx]
    gpr_paths <- find_gpr_paths(gpr_rule)
    for (path in gpr_paths){
      genes <- genes_from_path(path, rxn_idx)
      # print(genes)
      model <- normal_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id)
    }
    marked_rxns[rxn_idx] <- TRUE
  }
  
  print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  print(paste('remaining genes:', length(which(marked_genes == FALSE))))
  
  return(model)
}

gprRule_to_idx <- function(list){
  for (i in 1:length(list)){
    list[i] <- gsub("[^0-9]", "", list[i])
  }
  
  return(as.numeric(list))
}

# print(gprRule_to_idx(c('x[2]', 'x[5]', 'x[7]')))
model <- get_ecoli_model()
gene_sets <- ecoli_r0_gene_set
rxn_sets <- ecoli_og_set_list
test_model <- generate_falcon_model(model, ecoli_r0_gene_set, ecoli_og_set_list)
# model_2 <- generate_falcon_model(model, ecoli_r0_gene_set, ecoli_og_set_list)