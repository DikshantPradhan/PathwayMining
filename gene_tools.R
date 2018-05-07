library(dplyr)

#for (i in 1:length(Y_OPEN_GPR$react)){
#  genes <- strsplit(Y_OPEN_GPR$GPR[i], split = c("\\(|\\)|and|or| "))
#  genes <- genes[[1]][which(nchar(genes[[1]]) > 1)]
#  print(genes)
#  Y_OPEN_GPR$GENE[i] <- list(genes)
#}

# GPR <- mutans_gpr # Y_OPEN_GPR #Y4_GPR

generate_gpr <- function(model){
  GPR <- c()
  GPR$react <- model@react_id
  GPR$GENE <- model@genes
  GPR$GPR <- model@gpr

  return(GPR)
}

get_genes_from_react_id <- function(GPR, rxn_id){
  idx <- which(GPR$react == rxn_id) # react was previously Abbreviation
  #print(paste(idx, GPR$GENE[idx]))
  return(GPR$GENE[idx][[1]])
}

get_gene_pairs_from_rxn_pair <- function(rxn1, rxn2, GPR){
  g1 <- c()
  g2 <- c()

  gene1 <- get_genes_from_react_id(GPR,rxn1)
  gene2 <- get_genes_from_react_id(GPR,rxn2)

  if ((length(gene1) > 0) & (length(gene2) > 0)){
    for (i in gene1){
      for (j in gene2){
        if (!check_for_pairs(c(i, j), cbind(g1, g2))){
          g1 <- c(g1, i)
          g2 <- c(g2, j)
        }
      }
    }
  }

  return(cbind(g1, g2))
}

get_gene_pairs_from_rxn_pair_list <- function(pairs, GPR){
  g1 <- c()
  g2 <- c()

  for (i in 1:nrow(pairs)){
    #print(i)
    new_pairs <- get_gene_pairs_from_rxn_pair(pairs[i,1], pairs[i,2], GPR)
    if (length(new_pairs) > 0){
      for (j in 1:nrow(new_pairs)){
        #g1 <- c(g1, new_pairs[j,1])
        #g2 <- c(g2, new_pairs[j,2])
        if (!check_for_pairs(new_pairs[j,], cbind(g1,g2))){
          g1 <- c(g1, new_pairs[j,1])
          g2 <- c(g2, new_pairs[j,2])
        }
      }
    }
  }

  #genes <- remove_duplicate_pairs(cbind(g1, g2))
  genes <- cbind(g1, g2)

  return(genes)
}

check_for_enrichment <- function(gene_pairs, gi_e_matrix, threshold = 0.05){
  e_1 <- c()
  e_2 <- c()
  e <- c()

  for (i in 1:nrow(gene_pairs)){
    #print(i)

    gene_1 <- gene_pairs[i,1]
    gene_2 <- gene_pairs[i,2]

    gene_1_query <- grep(gene_1, rownames(gi_e_matrix))
    gene_2_query <- grep(gene_2, rownames(gi_e_matrix))
    gene_1_array <- grep(gene_1, colnames(gi_e_matrix))
    gene_2_array <- grep(gene_2, colnames(gi_e_matrix))

    e_ <- 0

    if (length(gene_1_query) > 0 & length(gene_2_array) > 0){
      for (i in gene_1_query){
        for (j in gene_2_array){
          e_ <- max(e_, abs(gi_e_matrix[gene_1_query, gene_2_array]))
        }
      }
    }

    if (length(gene_2_query) > 0 & length(gene_1_array) > 0){
      for (i in gene_2_query){
        for (j in gene_1_array){
          e_ <- max(e_, abs(gi_e_matrix[gene_2_query, gene_1_array]))
        }
      }
    }

    #e <- max(abs(gi_e_matrix[gene_1, gene_2]), abs(gi_e_matrix[gene_2, gene_1]))
    if (e_ >= threshold){
      e_1 <- c(e_1, gene_1)
      e_2 <- c(e_2, gene_2)
      e <- c(e, e_)
    }
  }

  percent <- length(e)/nrow(gene_pairs)
  mean <- mean(as.numeric(e))

  print(paste('percent: ', percent, ', mean: ', mean))

  return(cbind(e_1, e_2, e))

}

gene_set_from_rxn_set <- function(GPR, rxn_set_list){
  gene_set_list <- c()

  for (i in 1:length(rxn_set_list)){
    #print(i)
    gene_set <- c()
    for (rxn in rxn_set_list[[i]]){
      gene_set <- union(gene_set, get_genes_from_react_id(GPR,rxn))
    }
    if (length(gene_set) > 0){
      gene_set_list[i] <- list(gene_set)
    }
    else{
      print('empty gene set')
    }
  }

  return(gene_set_list)
}

gene_set_from_composition <- function(og_gene_set_list, composition){

  gene_set_list <- c()

  for (i in 1:length(composition)){
    gene_set <- c()

    for (j in composition[[i]]){
      if (j > length(og_gene_set_list)){next}
      gene_set <- c(gene_set, og_gene_set_list[[j]])
    }

    if (length(gene_set) > 0){
      gene_set_list[i] <- list(unique(gene_set))
    }
  }

  return(gene_set_list)

}

find_recurring_genes_in_set_list <- function(gene_set_list){
  recurring <- c()
  for (i in unique(unlist(gene_set_list))){
    recurring[i] <- 0
  }

  for (i in 1:length(gene_set_list)){
    for (j in unique(gene_set_list[[i]])){
      recurring[j] <- recurring[j] + 1
    }
  }

  return(recurring)
}

isolate_pairs_from_recurring_genes <- function(gene_pairs, gene_set_list, tol = 0){
  recurring_genes <- find_recurring_genes_in_set_list(gene_set_list)
  pairs <- isolate_pairs(names(recurring_genes)[which(recurring_genes >= tol)], gene_pairs)
  return(pairs)
}

# used one time to load data into usable form for R from file
build_gi_matrix <- function(gi_ExE, gi_NxN, gi_ExN_NxE){
  e_query_genes <- unique(gi_ExE$"Query Strain ID")
  n_query_genes <- unique(gi_NxN$"Query Strain ID")
  en_query_genes <- unique(gi_ExN_NxE$"Query Strain ID")

  e_array_genes <- unique(gi_ExE$"Array Strain ID")
  n_array_genes <- unique(gi_NxN$"Array Strain ID")
  en_array_genes <- unique(gi_ExN_NxE$"Array Strain ID")

  all_query_genes <- unique(c(e_query_genes, en_query_genes, n_query_genes))
  all_array_genes <- unique(c(e_array_genes, en_array_genes, n_array_genes))

  #total_n <- length(all_genes)

  gi_e_matrix <- matrix(0, nrow = length(all_query_genes), ncol = length(all_array_genes))
  gi_p_matrix <- matrix(0, nrow = length(all_query_genes), ncol = length(all_array_genes))

  rownames(gi_e_matrix) <- all_query_genes
  rownames(gi_p_matrix) <- all_query_genes
  colnames(gi_e_matrix) <- all_array_genes
  colnames(gi_p_matrix) <- all_array_genes

  print("ExE")
  for (i in 1:nrow(gi_ExE)){
    # print(i)
    query <- gi_ExE$"Query Strain ID"[i]
    array <- gi_ExE$"Array Strain ID"[i]
    e <- gi_ExE$"Genetic interaction score (ε)"[i]
    p <- gi_ExE$"P-value"[i]
    if (abs(e) > abs(gi_e_matrix[query, array])){
      gi_e_matrix[query, array] <- e
      gi_p_matrix[query, array] <- p
    }
  }

  print("NxN")
  for (i in 1:nrow(gi_NxN)){
    # print(i)
    query <- gi_NxN$"Query Strain ID"[i]
    array <- gi_NxN$"Array Strain ID"[i]
    e <- gi_NxN$"Genetic interaction score (ε)"[i]
    p <- gi_NxN$"P-value"[i]
    if (abs(e) > abs(gi_e_matrix[query, array])){
      gi_e_matrix[query, array] <- e
      gi_p_matrix[query, array] <- p
    }
  }

  print("ExN_NxE")
  for (i in 1:nrow(gi_ExN_NxE)){
    # print(i)
    query <- gi_ExN_NxE$"Query Strain ID"[i]
    array <- gi_ExN_NxE$"Array Strain ID"[i]
    e <- gi_ExN_NxE$"Genetic interaction score (ε)"[i]
    p <- gi_ExN_NxE$"P-value"[i]
    if (abs(e) > abs(gi_e_matrix[query, array])){
      gi_e_matrix[query, array] <- e
      gi_p_matrix[query, array] <- p
    }
  }

  gi_matrix <- c()
  gi_matrix$e <- gi_e_matrix
  gi_matrix$p <- gi_p_matrix

  return(gi_matrix)
}

# synthetic lethal pairs
build_essential_gi_matrix <- function(genes, lethal_pairs){
  genes <- c(genes, lethal_pairs[,1], lethal_pairs[,2])
  genes <- unique(genes)
  lethal_mtx <- matrix(data = 0, nrow = length(genes), ncol = length(genes))
  rownames(lethal_mtx) <- genes
  colnames(lethal_mtx) <- genes

  for (i in 1:(dim(lethal_pairs)[1])){
    #print(i)
    lethal_mtx[lethal_pairs[i,1], lethal_pairs[i,2]] <- 1
    lethal_mtx[lethal_pairs[i,2], lethal_pairs[i,1]] <- 1
  }

  return(lethal_mtx)
}

# build matrix of relevant genes for quicker computation
build_gi_matrix_from_pairs <- function(genes, gene_pairs){
  #gene_pairs <- remove_duplicate_pairs(gene_pairs)
  pairs <- gene_pairs[,1:2]
  e <- as.numeric(gene_pairs[,3])

  #genes <- unique(unlist(pairs))

  gi_mtx <- matrix(data = c(0), nrow = length(genes), ncol = length(genes))
  rownames(gi_mtx) <- genes
  colnames(gi_mtx) <- genes

  for (i in 1:length(pairs[,1])){
    pair <- pairs[i,]

    gi_mtx[pair[1], pair[2]] <- e[i]
    gi_mtx[pair[2], pair[1]] <- e[i]
  }

  return(gi_mtx)
}

# matrix of binary values (1, 0) based on whether or not gi pair is above or below threshold
binary_e_matrix <- function(full_matrix, threshold){
  e_matrix <- matrix(as.numeric(full_matrix >= threshold), nrow = nrow(full_matrix), ncol = ncol(full_matrix))
  rownames(e_matrix) <- rownames(full_matrix)
  colnames(e_matrix) <- colnames(full_matrix)
  return(e_matrix)
}

# generate and load relevant data into sigle structure
generate_gene_data <- function(og_set_list, set_lists){
  data <- c()

  all_pairs <- return_pairs_from_set(unique(unlist(og_set_list)))
  r0_pairs <- return_pairs_from_set_list(og_set_list)
  new_r1_pairs <- isolate_new_pairs_from_sets(og_set_list, set_lists)
  r1_pairs <- append_pair_lists(r0_pairs, new_r1_pairs)

  r1_set_list <- get_list_of_sets(r1_pairs)

  r0_gene_set_list <- gene_set_from_rxn_set(og_set_list)
  r1_gene_set_list <- gene_set_from_rxn_set(r1_set_list)

  all_gene_pairs <- return_pairs_from_set(unique(unlist(r0_gene_set_list)))
  r0_gene_pairs <- return_pairs_from_set_list(r0_gene_set_list)
  r1_gene_pairs <- return_pairs_from_set_list(r1_gene_set_list)
  new_r1_pairs <- isolate_new_pairs(r0_gene_pairs, r1_gene_pairs)

  #r0_gene_pairs <- get_gene_pairs_from_rxn_pair_list(r0_pairs)
  #r1_gene_pairs <- get_gene_pairs_from_rxn_pair_list(r1_pairs)
  #new_r1_gene_pairs <- get_gene_pairs_from_rxn_pair_list(new_r1_pairs)

  # load data
  data$r0_rxn_set_list <- og_set_list
  data$r1_rxn_set_list <- r1_set_list

  data$r0_gene_set_list <- r0_gene_set_list
  data$r1_gene_set_list <- r1_gene_set_list

  data$all_rxn_pairs <-
  data$r0_rxn_pairs <- r0_pairs
  data$r1_rxn_pairs <- r1_pairs
  data$new_r1_rxn_pairs <- new_r1_pairs

  data$all_gene_pairs <- all_gene_pairs
  data$r0_gene_pairs <- r0_gene_pairs
  data$r1_gene_pairs <- r1_gene_pairs
  data$new_r1_gene_pairs <- new_r1_gene_pairs

  return(data)
}

# find percent and mean enrichment values
multiple_enrichment_analysis <- function(all, r0, r1, new_r1, gi_e_matrix, e){
  all_enriched <- check_for_enrichment(all, gi_e_matrix, e)
  r0_enriched <- check_for_enrichment(r0, gi_e_matrix, e)
  r1_enriched <- check_for_enrichment(r1, gi_e_matrix, e)
  new_r1_enriched <- check_for_enrichment(new_r1, gi_e_matrix, e)

  all_percent <- nrow(all_enriched)/nrow(all)
  all_mean <- mean(as.numeric(all_enriched[,3]))

  r0_percent <- nrow(r0_enriched)/nrow(r0)
  r0_mean <- mean(as.numeric(r0_enriched[,3]))

  r1_percent <- nrow(r1_enriched)/nrow(r1)
  r1_mean <- mean(as.numeric(r1_enriched[,3]))

  new_r1_percent <- nrow(new_r1_enriched)/nrow(new_r1)
  new_r1_mean <- mean(as.numeric(new_r1_enriched[,3]))

  print(c('all gene pairs~', paste('total:', nrow(all_enriched), '/', nrow(all)),
    paste('percent:', all_percent),
    paste('mean:', all_mean)
  ))

  print(c('r0 gene pairs~', paste('total:', nrow(r0_enriched), '/', nrow(r0)),
    paste('percent:', r0_percent),
    paste('mean:', r0_mean)
  ))

  print(c('r1 gene pairs~', paste('total:', nrow(r1_enriched), '/', nrow(r1)),
    paste('percent:', r1_percent),
    paste('mean:', r1_mean)
  ))

  print(c('new_r1 gene pairs~', paste('total:', nrow(new_r1_enriched), '/', nrow(new_r1)),
    paste('percent:', new_r1_percent),
    paste('mean:', new_r1_mean)
  ))

  percent <- c(all_percent, r0_percent, r1_percent, new_r1_percent)
  mean <- c(all_mean, r0_mean, r1_mean, new_r1_mean)

  data <- cbind(percent, mean)

  colnames(data) <- c(paste(e, 'percent'), paste(e, 'mean'))
  rownames(data) <- c('all', 'r0', 'r1', 'new_r1')

  return(data)
}

enrichment_test_seq <- function(gene_data, gi_e_matrix, e_vals){
  all <- gene_data$all_pairs
  r0 <- gene_data$r0_pairs
  r1 <- gene_data$r1_pairs
  new_r1 <- gene_data$new_r1_pairs

  enrichment_data <- c()

  for (i in 1:length(e_vals)){
    e <- e_vals[i]
    data <- multiple_enrichment_analysis(all, r0, r1, new_r1, gi_e_matrix, e)
    enrichment_data[[i]] <- data
  }

  return(enrichment_data)
}

get_num_interactions_in_set <- function(set, e_matrix, clear = FALSE){
  row <- which(set %in% rownames(e_matrix))
  col <- which(set %in% colnames(e_matrix))
  #print(row_remove)
  #print(col_remove)
  row_set <- set[row]
  col_set <- set[col]
  if (length(row_set) < 2 & length(col_set) < 2){return(list(matrix = e_matrix, total = 0))}
  # Mp <- (abs(e_matrix) > threshold)
  # print(row_set)
  # print(col_set)
  sub_Mp <- e_matrix[row_set, col_set]
  # print(sub_Mp)
  # total <- sum(colSums(sub_Mp))# / 2
  total <- sum(sub_Mp)
  if (clear){
    e_matrix[row_set, col_set] <- 0
  }
  # Mp[r,r] <- 0
  # return(total)
  list(matrix = e_matrix, total = total)
}

# quicker running function for finding number of gene interactions
# e matrix is mtx of binary values, needs to differ for each threshold
get_num_interactions <- function(set_list, e_matrix, threshold = 0.0, clear = FALSE){

  # Mp <- e_matrix
  Mp <- (e_matrix >= threshold)
  
  total_int <- 0

  for (r in set_list) {
    # if (length(r) < 2){next}
    output <- get_num_interactions_in_set(r, Mp, clear = clear)
    total_int <- total_int + output$total
    Mp <- output$matrix
    # rm <- which(!(r %in% rownames(Mp)))
    # r <- r[-c(rm)]
    #print(r)
    # sub_Mp <- Mp[r,r]
    #print(sub_Mp)
    # total_int <- total_int + (sum(colSums(sub_Mp)) / 2)
    # Mp[r,r] <- 0
  }

  return(total_int)

}

# find genes with high expression change but low fitness change
identify_fitn_expr_relation <- function(fit_expr_mtx, fitn_thresh = 1, expr_thresh = 1){
  # col 1: genes, col 2: fitness fold change, col 3: expression fold change
  # assume log2 scaling

  interesting_genes <- matrix(data = FALSE, nrow = nrow(fit_expr_mtx), ncol = 1)

  for (i in 1:nrow(fit_expr_mtx)){
    fitn <- as.numeric(fit_expr_mtx[i,2])
    expr <- as.numeric(fit_expr_mtx[i,3])

    if (is.na(fitn) | is.na(expr)){next}

    if (near(fitn, 0, tol = fitn_thresh) & !near(expr, 0, tol = expr_thresh)){
      interesting_genes[i] <- TRUE
    }
  }

  return(fit_expr_mtx[which(interesting_genes),1])
}
