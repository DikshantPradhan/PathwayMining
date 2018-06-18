# PACT Analysis

# functions

g1_architecture_measurement <- function(coupling_mtx, df_frac, de_frac){
  edge_list <- matrix(nrow = length(which(coupling_mtx)), ncol = 3)
  set_idxs <- as.numeric(rownames(coupling_mtx))
  
  edge_idx <- 1
  for (i in 1:nrow(coupling_mtx)){
    for (j in which(coupling_mtx[i,])){
      set_idx_1 <- set_idxs[i]
      set_idx_2 <- set_idxs[j]
      df_diff <- abs(df_frac[set_idx_1] - df_frac[set_idx_2])
      de_diff <- abs(de_frac[set_idx_1] - de_frac[set_idx_2])
      euclidian_diff <- sqrt((df_frac[set_idx_1] - df_frac[set_idx_2])^2 + (de_frac[set_idx_1] - de_frac[set_idx_2])^2)
      
      # if (is.nan(df_diff)){
      #   df_diff <- 0
      #   de_diff <- 0
      #   euclidian_diff <- 0
      # }
      
      edge_list[edge_idx,1] <- df_diff
      edge_list[edge_idx,2] <- de_diff
      edge_list[edge_idx,3] <- euclidian_diff
      
      edge_idx <- edge_idx + 1
    }
  }
  
  return(edge_list)
}

g1_architecture_bootstrap <- function(coupling_mtx, df_frac, de_frac, n = 1000){
  sets <- rownames(coupling_mtx)
  
  boot_sets <- sample(sets, size = length(sets))
  rownames(coupling_mtx) <- boot_sets
  colnames(coupling_mtx) <- boot_sets
  
  edge_list <- g1_architecture_measurement(coupling_mtx, df_frac, de_frac)
  
  for (i in 2:n){
    boot_sets <- sample(sets, size = length(sets))
    rownames(coupling_mtx) <- boot_sets
    colnames(coupling_mtx) <- boot_sets
    
    edge_list <- cbind(edge_list, g1_architecture_measurement(coupling_mtx, df_frac, de_frac))
  }
  
  return(edge_list)
}

# initialization
source('~/GitHub/PathwayMining/data_tools.R')
# pao_g1_coupling <- read_coupling_csv('~/GitHub/PathwayMining/scripts/pao/pao_g1_coupling.csv')
# pao_g1_coupling_list <- pao_g1_coupling$coupling_vector
# 
# # pao_falcon <- generate_falcon_model(pao_model)
# pao_falcon <- GRB_pao_falcon_model()
# vars <- pao_falcon$get_names()$VarName #pao_falcon@react_id
# 
# full_pao_coupling_mtx <- coupling_matrix_from_coupling_vector_list(pao_g1_coupling_list, n_react = length(vars), vars = vars)
# fullish_pao_coupling_mtx <- full_ish_coupling_matrix_from_coupling_vector_list(pao_g1_coupling_list, n_react = length(vars), vars = vars, init_sets = g0_sets)
# 
# uncoupled_mtx <- identify_intermediate_uncoupled(full_pao_coupling_mtx, fullish_pao_coupling_mtx, n_react = length(vars))
# 
# uncoupled_pairs <- c()
# for (i in 1:nrow(uncoupled_mtx)){
#   for (j in which(uncoupled_mtx[i,])){
#     rxn_1 <- rownames(uncoupled_mtx)[i]
#     rxn_2 <- colnames(uncoupled_mtx)[j]
#     
#     # uncoupled_pairs <- c(uncoupled_pairs, list(c(rownames(uncoupled_mtx)[i], colnames(uncoupled_mtx)[j])))
#     if (grepl('PA', rxn_1) & grepl('PA', rxn_2)){
#       uncoupled_pairs <- c(uncoupled_pairs, list(c(rownames(uncoupled_mtx)[i], colnames(uncoupled_mtx)[j])))
#     }
#   }
#   # for (j in 1:ncol(uncoupled_mtx)){
#   #   if (uncoupled_mtx[i,j]){
#   #     uncoupled_pairs <- c(uncoupled_pairs, list(c(rownames(uncoupled_mtx)[i], colnames(uncoupled_mtx)[j])))
#   #   }
#   # }
# }
# 
# uncoupled_pairs <- clean_rxn_names_in_set(uncoupled_pairs)
# uncoupled_genes <- unique(unlist(uncoupled_pairs))

pao_data <- pao1#pa14[which(pa14$condition == unique(pa14$condition)[1]),]

# if using pa14 dataset, swap pa14 genes for pao1 homologs
# load("~/GitHub/PathwayMining/EnzymeCoupling/data/pao1_key.Rdata")
# pao_data$gene <- pao1_key[pao_data$gene]

# g0_set_df_og <- gene_set_dataframe(g0_sets)
# g1_set_df_og <- gene_set_dataframe(g1_sets)

genes_in_model <- unique(unlist(g0_df$clean.sets))

df_genes <- pao_data$gene[which(pao_data$sig_df)]
df_genes <- df_genes[which(df_genes %in% genes_in_model)]
de_genes <- pao_data$gene[which(pao_data$sig_de)]
de_genes <- de_genes[which(de_genes %in% genes_in_model)]
dfde_genes <- pao_data$gene[which(pao_data$sig_de & pao_data$sig_df)]
dfde_genes <- dfde_genes[which(dfde_genes %in% genes_in_model)]
num_df <- length(df_genes)
num_de <- length(de_genes)
num_dfde <- length(dfde_genes)

# coupled_genes <- setdiff(genes_in_model, uncoupled_genes)

composing_sets_list <- lapply(g1_df$sets, function(x){find_composing_sets(x, g0_df$sets)})
composing_sets_length_list <- lapply(composing_sets_list, function(x){length(x)})

g1_df$composing_g0_sets <- composing_sets_list
g1_df$num_composing_sets <- composing_sets_length_list

sets <- unlist(g1_df$composing_g0_sets[which(g1_df$num_composing_sets > 1)])

rxns <- unlist(g0_df$sets[sets])

# rxns <- c()
# set_names <- c()
# for (i in sets){
#   set <- unlist(g0_df$sets[i])
#   rxns <- c(rxns, set)
#   name <- paste(set, collapse = ', ')
#   set_names <- c(set_names, name)
# }
# 
# rxns <- unique(unlist(rxns))
g1_rxns <- unlist(g1_df$sets)
rxns <- intersect(rxns, g1_rxns)
network_coupling_mtx <- fullish_pao_coupling_mtx[rxns, rxns]

g0_set_coupling_mtx <- matrix(FALSE, nrow = length(sets), ncol = length(sets))
rownames(g0_set_coupling_mtx) <- sets
colnames(g0_set_coupling_mtx) <- sets

# for (i in 1:nrow(network_coupling_mtx)){
#   # print(i)
#   # print(length(which(network_coupling_mtx[i,])))
#   couplings <- unique(union(which(network_coupling_mtx[i,]), which(network_coupling_mtx[,i])))
#   if (length(couplings) < 2){print('nothing'); print(rxns[i])}
#   # print(length(which(network_coupling_mtx[i,])))
#   for (j in couplings){
#     v1 <- rxns[i]
#     v2 <- rxns[j]
#     if (grepl('\\(', v1)){
#       v1 <- strsplit(v1, split = '\\(')[[1]][1]
#     }
#     if (grepl('\\(', v2)){
#       v2 <- strsplit(v2, split = '\\(')[[1]][1]
#     }
#     if (grepl('\\[', v1)){
#       v1 <- strsplit(v1, split = '\\[')[[1]][1]
#     }
#     if (grepl('\\[', v2)){
#       v2 <- strsplit(v2, split = '\\[')[[1]][1]
#     }
# 
#     idx1 <- grep(grep(v1, g0_df$sets), sets)
#     idx2 <- grep(grep(v2, g0_df$sets), sets)
# 
#     g0_set_coupling_mtx[idx1, idx2] <- TRUE
#     g0_set_coupling_mtx[idx2, idx1] <- TRUE
#   }
# }

for (i in 1:length(sets)){
  for (j in i:length(sets)){
    set1 <- intersect(unlist(g0_df$sets[sets[i]]), g1_rxns)
    set2 <- intersect(unlist(g0_df$sets[sets[j]]), g1_rxns)
    
    # print('...')
    # print(set1)
    # print(set2)
    # 
    if (any(fullish_pao_coupling_mtx[set1,set2]) | any(fullish_pao_coupling_mtx[set2, set1])){
      # print('true')
      g0_set_coupling_mtx[i,j] <- TRUE
    }
  }
}

# compositions <- g1_df$composing_g0_sets[which(g1_df$num_composing_sets > 1)]
# sets <- unlist(compositions)
# g0_set_coupling_mtx <- matrix(FALSE, nrow = length(sets), ncol = length(sets))
# rownames(g0_set_coupling_mtx) <- sets
# colnames(g0_set_coupling_mtx) <- sets
# for (i in compositions){
#   idxs <- which(sets %in% i)
#   print(length(idxs))
#   g0_set_coupling_mtx[idxs, idxs] <- TRUE
# }
g0_set_coupling_mtx[lower.tri(g0_set_coupling_mtx, diag = TRUE)] <- FALSE

for (i in 1:length(sets)){
  couplings <- union(which(g0_set_coupling_mtx[i,]), which(g0_set_coupling_mtx[,i]))
  if (length(couplings) < 1){print(paste(i, 'nothing'))}
}

og_edge_list <- g1_architecture_measurement(g0_set_coupling_mtx, g0_df$df_frac, g0_df$de_frac)

edge_list <- g1_architecture_bootstrap(g0_set_coupling_mtx, g0_df$df_frac, g0_df$de_frac, 1000)

df_edges <- edge_list[, seq(from = 1, to = 3000, by = 3)]
de_edges <- edge_list[, seq(from = 2, to = 3000, by = 3)]
euc_edges <- edge_list[, seq(from = 3, to = 3000, by = 3)]

plot_density(df_edges, ylimits=c(0, 150))
freqs <- table(as.numeric(og_edge_list[,1]))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="df edges for g1 sets larger than 0 g0 sets", xlab="differential across edges", ylab="frequency")

plot_density(de_edges, ylimits=c(0, 150))
freqs <- table(as.numeric(og_edge_list[,2]))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="de edges for g1 sets larger than 0 g0 sets", xlab="differential across edges", ylab="frequency")

plot_density(euc_edges, xlimits = c(0,sqrt(2)), ylimits=c(0, 150))
freqs <- table(as.numeric(og_edge_list[,3]))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="df&de edges for g1 sets larger than 0 g0 sets", xlab="differential across edges", ylab="frequency")
