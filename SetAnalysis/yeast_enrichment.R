# load("~/GitHub/PathwayMining/gi_matrix.RData")
load('new_gi_matrix.RData')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g1_sets.RData")
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g0_sets.RData")
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_gr0_sets.RData")


## GENE DELETION FROM MODEL

#pairwiseGeneDel <- doubleGeneDel(yeast_model, allComb = TRUE)

## YEAST MODEL

# load("~/GitHub/PathwayMining/data/yeast_model/Maranas_model/maranas_model_lipid_exch.RData")
# model_genes <- yeast_model@allGenes

## ENRICHMENT ANALYSIS

# CLEAN E MATRIX
e_matrix <- gi_matrix$e

positive_e_matrix <- e_matrix
# negative_e_matrix <- -1*e_matrix
# abs_e_value <- abs(e_matrix)

clean_row_names <- function(names){
  for (i in 1:length(names)){
    names[i] <- strsplit(names[i], split = "_")[[1]][1]
  }
  return(names)
}

rownames(e_matrix) <- clean_row_names(rownames(e_matrix))
colnames(e_matrix) <- clean_row_names(colnames(e_matrix))

clean_duplicate_entries <- function(e_matrix, keep_unique = TRUE){
  
  intersect_genes <- intersect(unique(colnames(e_matrix)), unique(rownames(e_matrix)))
  print(length(intersect_genes))
  # reduced_row_e_matrix <- matrix(data = 0, nrow = length(intersect_genes), ncol = length(colnames(e_matrix)))
  
  if (!keep_unique){ # remove genes that are not in both rows and columns
    keep_rows <- which(rownames(e_matrix) %in% intersect_genes)
    keep_cols <- which(colnames(e_matrix) %in% intersect_genes)
    e_matrix <- e_matrix[keep_rows, keep_cols]
  }
  
  for (gene in intersect_genes){
    rows <- which(rownames(e_matrix) == gene)
    # global_max <- pmax(global_max, sol$X)
    if (length(rows) < 2){next}
    main <- rows[1]
    for (row in rows[2:length(rows)]){
      e_matrix[main,] <- pmax(e_matrix[main,], e_matrix[row])
    }
    
    e_matrix <- e_matrix[-rows[2:length(rows)],]
  }
  
  for (gene in intersect_genes){
    cols <- which(colnames(e_matrix) == gene)
    # global_max <- pmax(global_max, sol$X)
    if (length(cols) < 2){next}
    main <- cols[1]
    for (col in cols[2:length(cols)]){
      e_matrix[main,] <- pmax(e_matrix[main,], e_matrix[col])
    }
    e_matrix <- e_matrix[,-cols[2:length(cols)]]
  }
  
  # isolate recurring
  intersect_e <- e_matrix[intersect_genes, intersect_genes]
  intersect_e <- pmax(intersect_e, t(intersect_e))
  intersect_e[lower.tri(intersect_e,diag=TRUE)] <- 0
  
  e_matrix[intersect_genes, intersect_genes] <- intersect_e
  
  return(e_matrix)
}

isolate_relevant_interactions <- function(gi_matrix, sets){
  genes <- unique(unlist(sets))
  
  row_genes <- intersect(rownames(gi_matrix), genes)
  col_genes <- intersect(colnames(gi_matrix), genes)
  
  return(gi_matrix[row_genes, col_genes])
}

count_interactions <- function(sets, e_matrices, thresholds = seq(from = 0.1, to = 1, by = 0.1)){
  
  capture <- function(mtx, sets, thresholds){
    capture_rate <- matrix(data = NA, nrow = length(thresholds), ncol = 3)
    rownames(capture_rate) <- as.character(thresholds)
    colnames(capture_rate) <- c('total', 'captured', 'percent')
    
    for (i in 1:length(thresholds)){
      capture_rate[i,1] <- length(which(mtx > thresholds[i]))
      capture_rate[i,2] <- get_num_interactions(sets, mtx, thresholds[i], clear = TRUE)
      capture_rate[i,3] <- capture_rate[i,2]/capture_rate[i,1]
    }
    
    return(capture_rate)
  }
  
  pos_e_mtx <- isolate_relevant_interactions(e_matrices$positive, sets)
  neg_e_mtx <- isolate_relevant_interactions(e_matrices$negative, sets)
  abs_e_mtx <- isolate_relevant_interactions(e_matrices$absolute, sets)
  
  pos_capture <- capture(pos_e_mtx, sets, thresholds)
  neg_capture <- capture(neg_e_mtx, sets, thresholds)
  abs_capture <- capture(abs_e_mtx, sets, thresholds)
  
  list(positive = pos_capture, negatve = neg_capture, absolute = abs_capture)
}

positive_e_matrix <- clean_duplicate_entries(e_matrix)
negative_e_matrix <- clean_duplicate_entries(-1*e_matrix)
absolute_e_matrix <- clean_duplicate_entries(abs(e_matrix))

yeast_e_matrices <- list(positive = positive_e_matrix, negative = negative_e_matrix, absolute = absolute_e_matrix)
e_matrices <- yeast_e_matrices
# model and set intersect
set_genes <- unique(unlist(g1))
mtx_genes <- unique(c(rownames(positive_e_matrix), colnames(positive_e_matrix)))

common <- intersect(mtx_genes, set_genes)
# print(set_genes[which(!(set_genes %in% common))])
print(length(grep("_", set_genes[which(!(set_genes %in% common))])))

g0 <- clean_rxn_names_in_set(g0_sets)
gr0 <- yeast_gr0_sets
g1 <- clean_rxn_names_in_set(yeast_g1_sets)

g0_capture <- count_interactions(g0, yeast_e_matrices)
g1_capture <- count_interactions(g1, yeast_e_matrices)
gr0_capture <- count_interactions(gr0, yeast_e_matrices)

yeast_interaction_capture <- list(g0 = g0_capture, g1 = g1_capture, gr0 = gr0_capture)

# 
g1_interact <- get_num_interactions(g1, positive_e_matrix, 0.1)
print(g1_interact)
g1_interact <- get_num_interactions(g1, negative_e_matrix, 0.1)
print(g1_interact)
g1_interact <- get_num_interactions(g1, absolute_e_matrix, 0.1)
print(g1_interact)

test <- get_num_interactions_in_set(g1[[986]], (positive_e_matrix >= 0))

print(get_num_interactions_in_set(g1[[986]], (positive_e_matrix >= 0)))

# 
# pairs <- return_pairs_from_set_list(g1)
# enrichment <- check_for_enrichment(pairs, negative_e_matrix, threshold = 0.1)
# 
