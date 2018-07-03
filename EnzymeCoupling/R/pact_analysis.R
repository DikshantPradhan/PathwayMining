# PACT Analysis
library(ggplot2)
source('~/GitHub/PathwayMining/set_tools.R')
load("~/GitHub/PathwayMining/EnzymeCoupling/PA_FitExp.Rdata")
load("~/GitHub/PathwayMining/EnzymeCoupling/data/pao_g1_coupling_mtxs.RData")

# initialization

pao_data <- pao1#pa14[which(pa14$condition == unique(pa14$condition)[1]),]

sig_threshold <- 0.5

sig_df_frac <- g0_df$df_frac >= sig_threshold
sig_de_frac <- g0_df$de_frac >= sig_threshold

# sig_df_frac[is.na(sig_df_frac)] <- FALSE
# sig_de_frac[is.na(sig_de_frac)] <- FALSE

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
# sets <- intersect(sets, which(g0_df$X..genes >= 1))

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

og_edge_list <- g1_architecture_measurement_binary(g0_set_coupling_mtx, sig_df_frac, sig_de_frac)

# edge_list <- g1_architecture_bootstrap(g0_set_coupling_mtx, g0_df$df_frac, g0_df$de_frac, 1000)
edge_list <- g1_architecture_bootstrap(g0_set_coupling_mtx, sig_df_frac, sig_de_frac, 1000)

df_edges <- edge_list[, seq(from = 1, to = 3000, by = 3)]
de_edges <- edge_list[, seq(from = 2, to = 3000, by = 3)]
euc_edges <- edge_list[, seq(from = 3, to = 3000, by = 3)]

binary_edge_histogram(df_edges, og_edge_list[,1], title = '', color = 'lightskyblue') # df-df_edges (threshold fraction = 0.50)
binary_edge_histogram(de_edges, og_edge_list[,2], title = '', color = 'chartreuse3') # de-de_edges (threshold fraction = 0.50)
binary_edge_histogram(euc_edges, og_edge_list[,3], title = '', color = 'mediumorchid') # df-de_edges (threshold fraction = 0.50)

df_degree <- node_degree(g0_set_coupling_mtx, g0_df$df_frac)
de_degree <- node_degree(g0_set_coupling_mtx, g0_df$de_frac)

plot(jitter(df_degree[2,]), jitter(df_degree[1,]), main = 'df node degree', xlab = 'df frac', ylab = 'degree')
plot(jitter(de_degree[2,]), jitter(de_degree[1,]), main = 'de node degree', xlab = 'de frac', ylab = 'degree')

plot_density(df_edges, xlimits = c(-0.15,1), ylimits=c(0, 75))
freqs <- table(as.numeric(og_edge_list[,1]))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="df edges for g1 sets larger than 1 g0 sets", xlab="differential across edges", ylab="frequency")

plot_density(de_edges, xlimits = c(-0.15,1), ylimits=c(0, 75))
freqs <- table(as.numeric(og_edge_list[,2]))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="de edges for g1 sets larger than 1 g0 sets", xlab="differential across edges", ylab="frequency")

plot_density(euc_edges, xlimits = c(-0.15,sqrt(2)), ylimits=c(0, 75))
freqs <- table(as.numeric(og_edge_list[,3]))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="df&de edges for g1 sets larger than 1 g0 sets", xlab="differential across edges", ylab="frequency")

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
