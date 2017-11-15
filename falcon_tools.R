library(pryr)
library(data.tree)
source('~/GitHub/PathwayMining/logic_tools.R')

generate_falcon_model <- function(model, gene_sets = c(), rxn_sets = c()){
  # gprRules <- model@gprRules
  # genes <- model@genes
  og_dim <- dim(model@S) # not necessary
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

  rxn_removal_ids <- c()

  genes_from_path <- function(path, rxn_idx){
    gene_idxs <- gprRule_to_idx(path)
    genes <- og_genes[rxn_idx]
    return(genes[[1]][gene_idxs])
  }

  normal_add <- function(model, new_met_list, rxn_id, simple = FALSE, addExch = FALSE, identifier = NULL){
    # print(rxn_id)
    # print(new_met_list)
    rxn_idx <- which(model@react_id == rxn_id)
    # print(rxn_idx)
    exch <- findExchReact(model)

    # metabolites of existing reaction
    old_met_idxs <- which(model@S[, rxn_idx] != 0)
    old_met_idxs <- old_met_idxs[which(old_met_idxs <= og_dim[1])]
    old_met_list <- og_met_id[old_met_idxs]
    old_met_coeff <- model@S[old_met_idxs, rxn_idx]
    
    # add exchange reactions, if needed
    if (addExch){
      for (met in new_met_list){
        if (met %in% exch@met_id){next}
        model <- addExchReact(model, met, -1000, 1000)
      }
    }

    # check to see if rxn already exists
    # new_react <- matrix(0, nrow = dim(model@S)[1], ncol = 1)
    # rownames(new_react) <- model@met_id
    # new_react <- model@S[, rxn_idx]
    # new_react[new_met_list] <- -1 # indexing error here

    # for (i in 1:dim(model@S)[2]){
    #   if (identical(new_react, model@S[, i])){
    #     return(model)
    #   }
    # }
    met_list <- c(unlist(old_met_list), unlist(new_met_list))
    met_list <- met_list[!is.na(met_list)]
    # print(met_list)
    # print(new_met_list)
    # print(old_met_coeff)
    
    if (!simple){ # add reverse reaction if needed
      # print(rxn_id)
      # print(c(old_met_list, new_met_list))
      model <- addReact(model, paste(rxn_id, identifier, 'fwd', sep = '_'), met = met_list,
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = 0, ub = 1000, reversible = FALSE)
      model <- addReact(model, paste(rxn_id, identifier, 'rev', sep = '_'), met = met_list,
                        Scoef = c(old_met_coeff, rep(1, length(new_met_list))), lb = -1000, ub = 0, reversible = FALSE)
      rxn_removal_ids <- c(rxn_removal_ids, rxn_id)
      # if (rxn_id %in% model@react_id){
      #   model <- rmReact(model = model, rxn_id, rm_met = FALSE)
      # }
      ## NEED CONSTRAINTS TO PREVENT MODEL FROM PUSHING FLUX THROUGH BOTH DIRECTIONS AT ONCE
    }
    else {
      model <- addReact(model, rxn_id, met = met_list,
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = -1000, ub = 1000, reversible = TRUE)
    }

    # rxn_removal_ids <- c(rxn_removal_ids, rxn_id)
    return(model)
  }

  simple_add <- function(model, new_met_list, rxn_id, simple = FALSE){
    # print(c(rxn_id, ":", new_met_list))
    model <- normal_add(model, new_met_list, rxn_id, simple = TRUE)
    return(model)
  }

  or_add <- function(model, path_list, rxn_id, simple = TRUE){
    # print(paste('or add:', rxn_id))
    rxn_idx <- which(og_react_id == rxn_id)
    rxn_activity <- paste('a', rxn_id, sep = "_")
    # if (rxn_activity == 'a_'){print(paste('blank a test:', rxn_id))}
    # print(path_list)
    # add conversion for path to react_activity
    identifier <- 1
    for (mets in path_list){
      # print(mets)
      genes <- genes_from_path(mets, rxn_idx)
      # print(genes)
      new_mets <- paste('a_', genes, sep = '')
      # print(new_mets)
      # if ('a_' %in% new_mets){print(rxn_id)}
      met_list <- c(unlist(new_mets), rxn_activity)
      met_list <- met_list[!is.na(met_list)]
      coeff_list <- c(unlist(rep(-1, length(new_mets))), 1)
      model <- addReact(model, paste(rxn_activity, 'conversion', identifier, sep = ' '), met = met_list,
                        Scoef = coeff_list, lb = -1000, ub = 1000, reversible = TRUE)
      identifier <- identifier + 1
    }

    # add activity specific to react
    model <- normal_add(model, new_met_list = c(rxn_activity), rxn_id, simple = simple, addExch = FALSE)

    ## NEED TO ADD FLUX CONSTRAINT TO THIS?

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
    rxn_id <- og_react_id[i]
    # print(model@react_id[i])
    gpr_rule <- model@gprRules[i]
    # print(gpr_rule)
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (length(gpr_paths) == 1){
      genes <- genes_from_path(gpr_paths[[1]], i)
      # print(genes)
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
      # print(genes)
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

  print('SET EXCLUSIVITY')

  print(paste('loyal sets:', length(loyal_set_idxs)))
  print(paste('set exclusive genes:', length(set_exclusive_genes)))

  for (i in loyal_set_idxs){
    # print(paste('new set:', length(rxn_sets[[i]])))
    # print(gene_sets[[i]])
    for (j in rxn_sets[[i]]){
      # print(j)
      rxn_idx <- which(model@react_id == j)
      # print('set excl')
      gpr_rule <- model@gprRules[rxn_idx]
      gpr_paths <- find_gpr_paths(gpr_rule)
      # print(gpr_paths)
      genes <- c()
      for (path in gpr_paths){
        genes <- c(genes, genes_from_path(path, rxn_idx))
      }
      genes <- unlist(genes)
      # genes <- genes_from_path(gpr_paths[[1]], rxn_idx)
      if (which(og_react_id == j) %in% which(marked_rxns)){
        # print(paste('loyal rxn;', model@gpr[rxn_idx]))
        next
      }
      if (nchar(gpr_paths[1]) == 0){
        # print(gpr_paths)
        next
      }
      # print(genes)
      # print(model@gpr[rxn_idx])
      model <- or_add(model, gpr_paths, j)
      marked_rxns[rxn_idx] <- TRUE
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
    }

    # marked_genes[which(model@allGenes %in% gene_sets[[i]])] <- TRUE
    # marked_rxns[which(og_react_id %in% rxn_sets[[i]])] <- TRUE
  }

  print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  print(paste('remaining genes:', length(which(marked_genes == FALSE))))
  
  }
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
    # print(rxn_idx)
    rxn_id <- og_react_id[rxn_idx]
    # print(rxn_id)
    gpr_rule <- model@gprRules[rxn_idx]
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (nchar(gpr_paths[1]) == 0){next}
    id <- 1
    for (path in gpr_paths){
      genes <- genes_from_path(path, rxn_idx)
      # print(genes)
      # print(rxn_id)
      model <- normal_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id, identifier = id)
      rxn_removal_ids <- c(rxn_removal_ids, rxn_id)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      id <- id + 1
    }
    marked_rxns[rxn_idx] <- TRUE
  }

  print('FINAL COUNT')
  print(paste('remaining reactions:', length(which(marked_rxns == FALSE))))
  print(paste('remaining genes:', length(which(marked_genes == FALSE))))

  # print(marked_genes)
  
  print('RXN REMOVAL')
  # print(rxn_removal_ids)
  print(paste('met num:', length(model@met_id)))
  for (rxn_id in unique(rxn_removal_ids)){
    # print(rxn_id)
    model <- rmReact(model, react = rxn_id, rm_met = FALSE)
  }
  print(paste('met num:', length(model@met_id)))

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
test_model2 <- generate_falcon_model(model)
a <- apply(test_model@S, 1, function(x) sum(abs(x)))
print(test_model@met_id[which(a == 0)])

# test_falcon <- function(model, idx1, idx2 = 0){
#   if (idx2 == 0){idx2 <- idx1 + 1}
#   
#   rxn1_all <- which(test_model@S[,idx1] != 0)
#   rxn2_all <- which(test_model@S[,idx2] != 0)
#   
#   rxn1_pdt <- which(test_model@S[,idx1] > 0)
#   rxn2_pdt <- which(test_model@S[,idx2] > 0)
#   
#   rxn1_rct <- which(test_model@S[,idx1] < 0)
#   rxn2_rct <- which(test_model@S[,idx2] < 0)
#   
#   rct_diff <- unique(union(setdiff(rxn1_rct, rxn2_rct), setdiff(rxn2_rct, rxn1_rct)))
#   pdt_diff <- unique(union(setdiff(rxn2_pdt, rxn1_pdt), setdiff(rxn1_pdt, rxn2_pdt)))
#   
#   if (!identical(rxn1_all, rxn2_all)){return(FALSE)}
#   
#   if (!identical(rct_diff, pdt_diff)){return(FALSE)}
#   
#   return(TRUE)
# }
# 
# test_seq <- seq(267, 321, 2)
# for (idx in test_seq){
#   print(test_falcon(test_model, idx))
# }
# 
# print('TESTING CONVERSION REACTIONS')
# 
# conv_rxns <- c("ACKr","ATPS4r","CYTBD","D_LACt2","FBA","FBP","FUM","GLNS","LDH_D","MALS","PFK","PFL","PGM","PIt2r","PTAr",   
#                "PYK","RPE","RPI","TALA","ACONTa","ACONTb","TKT1","TKT2")
# 
# # REACTANTS AND PRODUCTS SHOULD SHARE SAME DIFFERING METABOLITE
# for (rxn in conv_rxns){
#   print(paste('new rxn:', rxn))
#   og_idx <- which(model@react_id == rxn)
#   print(model@gpr[og_idx])
#   gpr_rule <- model@gprRules[og_idx]
#   gpr_paths <- find_gpr_paths(gpr_rule)
#   
#   new_rxns <- grep(paste(rxn, 'conversion'), test_model@react_id)
#   if (length(new_rxns) != length(gpr_paths)){print(paste('error', rxn))}
#   
#   for (new_rxn in new_rxns){
#     met_idxs <- which(test_model@S[,new_rxn] != 0)
#     print(paste(test_model@met_id[met_idxs], test_model@S[met_idxs, new_rxn]))
#   }
#   
#   new_idx <- which(test_model@react_id == rxn)
#   met_idxs <- which(test_model@S[,new_idx] != 0)
#   print(paste(test_model@met_id[met_idxs], test_model@S[met_idxs, new_idx]))
# }
# 
# print('TEST ALL SIMPLE RXNS')
# # PRINT ALL METABOLTES AND COEFFICIENTS
# for (rxn in model@react_id){
#   print(paste('new rxn:', rxn))
#   old_idx <- which(model@react_id == rxn)
#   print(model@gpr[old_idx])
#   
#   new_idx <- which(test_model@react_id == rxn)
#   if (length(new_idx) > 0){
#     met_idxs <- which(test_model@S[,new_idx] != 0)
#     print(paste(test_model@met_id[met_idxs], test_model@S[met_idxs, new_idx]))
#   }
#   else {
#     print(paste('deleted rxn:', (length(grep(rxn, test_model@react_id[267:324])) > 1)))
#   }
# }
# 
# print('METABOLITE TEST')
# # ORIGINAL REACTIONS SHOULD SHARE SAME METABOLITE AND COEFFICIENTS
# for (rxn in model@react_id){
#   print(paste('new rxn:', rxn))
#   old_idx <- which(model@react_id == rxn)
#   old_met_idxs <- which(model@S[,old_idx] != 0)
#   old_met_coeffs <- model@S[old_met_idxs, old_idx]
#   
#   new_idxs <- grep(rxn, test_model@react_id)
#   
#   for (new_idx in new_idxs){
#     new_rxn_id <- test_model@react_id[new_idx]
#     # print(new_rxn_id)
#     new_met_idxs <- which(test_model@S[,new_idx] != 0)
#     new_met_coeffs <- test_model@S[old_met_idxs, new_idx]
#     if (!all(old_met_idxs %in% new_met_idxs)){
#       print(paste(new_rxn_id, 'error'))
#       print(old_met_idxs)
#       print(new_met_idxs)
#     }
#     if (!identical(old_met_coeffs, new_met_coeffs)){
#       print(paste(new_rxn_id, 'error2'))
#       print(old_met_coeffs)
#       print(new_met_coeffs)
#     }
#   }
# }

# model_2 <- generate_falcon_model(model, ecoli_r0_gene_set, ecoli_og_set_list)
