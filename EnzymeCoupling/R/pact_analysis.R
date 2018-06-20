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
      
      # g1_set <- get_set_idx(set_idxs[i], g1_composition)
      # print(g1_set)
      if (is.nan(df_diff)){
        df_diff <- -0.1
        de_diff <- -0.1
        euclidian_diff <- -0.1
      }
      
      edge_list[edge_idx,1] <- df_diff
      edge_list[edge_idx,2] <- de_diff
      edge_list[edge_idx,3] <- euclidian_diff
      # edge_list[edge_idx,4] <- g1_set
      
      edge_idx <- edge_idx + 1
    }
  }
  
  return(edge_list)
}

g1_architecture_measurement_binary <- function(coupling_mtx, df_frac, de_frac){
  edge_list <- matrix(nrow = length(which(coupling_mtx)), ncol = 3)
  set_idxs <- as.numeric(rownames(coupling_mtx))
  
  edge_idx <- 1
  for (i in 1:nrow(coupling_mtx)){
    for (j in which(coupling_mtx[i,])){
      set_idx_1 <- set_idxs[i]
      set_idx_2 <- set_idxs[j]
      df_diff <- df_frac[set_idx_1]&df_frac[set_idx_2]
      de_diff <- de_frac[set_idx_1]&de_frac[set_idx_2]
      euclidian_diff <- (df_frac[set_idx_1]|de_frac[set_idx_1])&(df_frac[set_idx_2]&de_frac[set_idx_2])
      
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
  
  edge_list <- g1_architecture_measurement_binary(coupling_mtx, df_frac, de_frac)
  
  for (i in 2:n){
    boot_sets <- sample(sets, size = length(sets))
    rownames(coupling_mtx) <- boot_sets
    colnames(coupling_mtx) <- boot_sets
    
    edge_list <- cbind(edge_list, g1_architecture_measurement(coupling_mtx, df_frac, de_frac))
  }
  
  return(edge_list)
}

binary_edge_histogram <- function(edge_lists, og_edges, xlimits = c(0,50)){
  
  bootstrapped <- c()
  for (i in 1:ncol(edge_lists)){
    bootstrapped <- c(bootstrapped, length(which(edge_lists[,i]==0)))
  }
  
  observed <- length(which(og_edges==0))
  
  qplot(bootstrapped, geom = 'histogram', bins = 50, main = 'title', xlim = xlimits) + 
    geom_vline(xintercept = observed, colour = 'red') + xlab('') + ylab('count')
}

# initialization

pao_data <- pao1#pa14[which(pa14$condition == unique(pa14$condition)[1]),]

sig_df_frac <- g0_df$df_frac >= 0.5
sig_de_frac <- g0_df$de_frac >= 0.5

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

composing_sets_list <- lapply(g1_df$sets, function(x){find_composing_sets(x, g0_df$sets)})
composing_sets_length_list <- lapply(composing_sets_list, function(x){length(x)})

g1_df$composing_g0_sets <- composing_sets_list
g1_df$num_composing_sets <- composing_sets_length_list
# g1_df <- g1_df[order(g1_df$num_composing_sets),]
sets <- unlist(g1_df$composing_g0_sets[which(g1_df$num_composing_sets > 1)])
sets <- intersect(sets, which(g0_df$X..genes >= 1))

rxns <- unlist(g0_df$sets[sets])

g1_rxns <- unlist(g1_df$sets)
rxns <- intersect(rxns, g1_rxns)
network_coupling_mtx <- fullish_pao_coupling_mtx[rxns, rxns]

g0_set_coupling_mtx <- matrix(FALSE, nrow = length(sets), ncol = length(sets))
rownames(g0_set_coupling_mtx) <- sets
colnames(g0_set_coupling_mtx) <- sets

for (i in 1:length(sets)){
  for (j in i:length(sets)){
    set1 <- intersect(unlist(g0_df$sets[sets[i]]), g1_rxns)
    set2 <- intersect(unlist(g0_df$sets[sets[j]]), g1_rxns)
    
    if (any(fullish_pao_coupling_mtx[set1,set2]) | any(fullish_pao_coupling_mtx[set2, set1])){
      g0_set_coupling_mtx[i,j] <- TRUE
    }
  }
}

g0_set_coupling_mtx[lower.tri(g0_set_coupling_mtx, diag = TRUE)] <- FALSE

# for (i in 1:length(sets)){
#   couplings <- union(which(g0_set_coupling_mtx[i,]), which(g0_set_coupling_mtx[,i]))
#   if (length(couplings) < 1){print(paste(i, 'nothing'))}
# }

# og_edge_list <- g1_architecture_measurement(g0_set_coupling_mtx, g0_df$df_frac, g0_df$de_frac)#, g1_df$composing_g0_sets)

og_edge_list <- g1_architecture_measurement(g0_set_coupling_mtx, sig_df_frac, sig_de_frac)

# edge_list <- g1_architecture_bootstrap(g0_set_coupling_mtx, g0_df$df_frac, g0_df$de_frac, 1000)
edge_list <- g1_architecture_bootstrap(g0_set_coupling_mtx, sig_df_frac, sig_de_frac, 1000)

df_edges <- edge_list[, seq(from = 1, to = 3000, by = 3)]
de_edges <- edge_list[, seq(from = 2, to = 3000, by = 3)]
euc_edges <- edge_list[, seq(from = 3, to = 3000, by = 3)]

binary_edge_histogram(df_edges, og_edge_list[,1], xlimits = c(0,100))

# plot_density(df_edges, xlimits = c(-0.15,1), ylimits=c(0, 75))
# freqs <- table(as.numeric(og_edge_list[,1]))
# freqs_x <- as.numeric(names(freqs))
# lines(freqs_x, freqs, col = 'red')
# title(main="df edges for g1 sets larger than 1 g0 sets", xlab="differential across edges", ylab="frequency")
# 
# plot_density(de_edges, xlimits = c(-0.15,1), ylimits=c(0, 75))
# freqs <- table(as.numeric(og_edge_list[,2]))
# freqs_x <- as.numeric(names(freqs))
# lines(freqs_x, freqs, col = 'red')
# title(main="de edges for g1 sets larger than 1 g0 sets", xlab="differential across edges", ylab="frequency")
# 
# plot_density(euc_edges, xlimits = c(-0.15,sqrt(2)), ylimits=c(0, 75))
# freqs <- table(as.numeric(og_edge_list[,3]))
# freqs_x <- as.numeric(names(freqs))
# lines(freqs_x, freqs, col = 'red')
# title(main="df&de edges for g1 sets larger than 1 g0 sets", xlab="differential across edges", ylab="frequency")

# og_edge_list <- og_edge_list[order(og_edge_list[,4]),]
# og_edge_list[,4] <- vapply(og_edge_list[,4], function(x){which(unique(og_edge_list[,4]) == x)}, c(1))
# og_edge_list <- og_edge_list[!is.nan(og_edge_list[,1]),]
# i <- 1
# plot(og_edge_list[,4], og_edge_list[,i])
# lines(unique(og_edge_list[,4]), vapply(unique(og_edge_list[,4]), function(x){mean(og_edge_list[which(og_edge_list[,4]==x),1])}, c(1)))
# 
# # plot.new()
# col <- 1
# points <- matrix(data = NA, nrow = length(sets), ncol = 2)
# i <- 1
# for (idx in which(g1_df$num_composing_sets > 1)){
#   for (g0_set in g1_df$composing_g0_sets[[idx]]){
#     # print(paste(idx, g0_set))
#     # print(idx)
#     points[i,1] <- idx
#     points[i,2] <- g0_df$df_frac[g0_set]
#     # points(idx, g0_df$df_frac[g0_set])
#     i <- i+1
#   }
# }
# 
# points <- points[!is.nan(points[,2]),]
# # points <- points[!(points[,2] == 0),]
# points[,1] <- vapply(points[,1], function(x){which(unique(points[,1]) == x)}, c(1))
# plot(x = points[,1], y = jitter(points[,2], amount = 0.02), main="set df_fracs and average edge diff", xlab="set", ylab="fracs and average edge diff")
# lines(unique(og_edge_list[,4]), vapply(unique(og_edge_list[,4]), function(x){mean(og_edge_list[which(og_edge_list[,4]==x),col])}, c(1)))
# 
