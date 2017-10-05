
# for (i in 1:nrow(Y4_GPR)){ 
#   genes <- strsplit(Y4_GPR$GPR[i], split = c("\\(|\\)|and|or"))
#   genes <- genes[[1]][which(nchar(genes[[1]]) > 1)]
#   # print(genes)
#   Y4_GPR$GENE[i] <- list(genes)
# }

GPR <- Y4_GPR

get_genes_from_react_id <- function(rxn_id){
  idx <- which(GPR$Abbreviation == rxn_id)
  return(GPR$GENE[idx][[1]])
}

get_gene_pairs_from_rxn_pair <- function(rxn1, rxn2){
  g1 <- c()
  g2 <- c()
  
  gene1 <- get_genes_from_react_id(rxn1)
  gene2 <- get_genes_from_react_id(rxn2)
  
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

get_gene_pairs_from_rxn_pair_list <- function(pairs){
  g1 <- c()
  g2 <- c()
  
  for (i in 1:nrow(pairs)){
    print(i)
    new_pairs <- get_gene_pairs_from_rxn_pair(pairs[i,1], pairs[i,2])
    if (length(new_pairs) > 0){
      for (j in 1:nrow(new_pairs)){
        g1 <- c(g1, new_pairs[j,1])
        g2 <- c(g2, new_pairs[j,2])
        # if (!check_for_pairs(new_pairs[j,], cbind(g1,g2))){
        #   g1 <- c(g1, new_pairs[j,1])
        #   g2 <- c(g2, new_pairs[j,2])
        # }
      }
    }
  }
  
  genes <- remove_duplicate_pairs(cbind(g1, g2))
  
  return(genes)
}

check_for_enrichment <- function(gene_pairs, lethal_pairs){
  e_1 <- c()
  e_2 <- c()
  
  for (i in 1:nrow(gene_pairs)){
    print(i)
    if (check_for_pairs(gene_pairs[i,], lethal_pairs)){
      e_1 <- c(e_1, gene_pairs[i,1])
      e_2 <- c(e_2, gene_pairs[i,2])
    }
  }
  
  return(cbind(e_1, e_2))
}

gene_set_from_rxn_set <- function(rxn_set_list){
  gene_set_list <- c()
  
  for (i in 1:length(rxn_set_list)){
    gene_set <- c()
    for (rxn in rxn_set_list[[i]]){
      gene_set <- union(gene_set, get_genes_from_react_id(rxn))
    }
    if (length(gene_set) > 0){
      gene_set_list[i] <- list(gene_set)
    }
  }
  
  return(gene_set_list)
}

find_recurring_genes_in_set_list <- function(gene_set_list){
  recurring <- c()
  for (i in unique(unlist(gene_set_list))){
    recurring[i] <- 0
  }
  
  for (i in 1:(length(gene_set_list)-1)){
    for (j in (i+1):length(gene_set_list)){
      repeating_genes <- intersect(gene_set_list[[i]], gene_set_list[[j]])
      
      if (length(repeating_genes) > 0){
        # print(repeating_genes)
        for (k in repeating_genes){
          recurring[k] <- recurring[k] + 1
        }
      }
      
    }
  }
  
  return(recurring)
}