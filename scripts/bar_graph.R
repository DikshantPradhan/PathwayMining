
bar_graph_dist <- function(r0_dist, r1_dist){
  len <- max(nrow(r0_dist), nrow(r1_dist))
  dist <- matrix(0, nrow = 2, ncol = len)
  rownames(dist) <- c('r0', 'r1')
  colnames(dist) <- c(1:len)
  
  for (i in 1:nrow(r0_dist)){
    dist[1,i] <- r0_dist[i]
  }
  for (i in 1:nrow(r1_dist)){
    dist[2,i] <- r1_dist[i]
  }
  
  return(dist)
}

rxn_dist <- bar_graph_dist(r0_rxn_dist, r1_rxn_dist)
gene_dist <- bar_graph_dist(r0_gene_dist, r1_gene_dist)

bar_graph_dist_recurr <- function(r0_dist, r1_dist){
  len <- max(as.numeric(names(r0_dist)), as.numeric(names(r1_dist)))
  dist <- matrix(0, nrow = 2, ncol = len)
  rownames(dist) <- c('r0', 'r1')
  colnames(dist) <- c(1:len)
  
  for (i in names(r0_dist)){
    dist[1, as.numeric(i)] <- r0_dist[i]
  }
  for (i in names(r1_dist)){
    dist[2, as.numeric(i)] <- r1_dist[i]
  }
  
  return(dist)
}

recurrence_dist <- bar_graph_dist_recurr(table(r0_recurring_genes), table(r1_recurring_genes))

# barplot(gene_dist[,6:28], beside = TRUE, main = 'Gene Set Comparison', xlab = 'Size of Set', ylab = 'Frequency of Size', legend = rownames(gene_dist))
# barplot(rxn_dist[,5:52], beside = TRUE, main = 'Reaction Set Comparison', xlab = 'Size of Set', ylab = 'Frequency of Size', legend = rownames(rxn_dist))
barplot(recurrence_dist[,4:22], beside = TRUE, main = 'Gene Promiscuity Comparison', xlab = 'Degree of Promiscuity', ylab = 'Frequency of Promiscuity', legend = rownames(recurrence_dist))
