library(genbankr)
library(sybil)
library(Matrix)

# get_gene_dataframe <- function(gb, model){
#   genes <- gb@genes$gene_id
#   transcripts
# }

get_gene_order_list <- function(gb){
  start_positions <- gb@transcripts@ranges@start
  gene_positions <- sort(start_positions)
  # return(gene_positions)
  
  gene_list <- c()
  for (i in gene_positions){
    idx <- which(start_positions == i)
    new_gene <- gb@transcripts$locus_tag[idx]
    gene_list <- c(gene_list, new_gene)
  }
  
  return(gene_list)
}

get_genome_distance <- function(gene_order_list){
  get_genome_distance <- function(a, b){
    idx_a <- which(gene_order_list == a)
    idx_b <- which(gene_order_list == b)
    dist <- abs(idx_a - idx_b)
    return(dist)
  }
  
  n_genes <- length(gene_order_list)
  distance <- array(data = NA, dim = c(n_genes, n_genes))
  rownames(distance) <- gene_order_list
  colnames(distance) <- gene_order_list
  
  for (i in 1:(length(gene_order_list)-1)){
    gene_i <- gene_order_list[i]
    distance[i,i] <- 0
    for (j in (i+1):length(gene_order_list)){
      gene_j <- gene_order_list[j]
      dist <- get_genome_distance(gene_i, gene_j)
      distance[i,j] <- dist
      distance[j,i] <- dist
    }
  }
  distance[n_genes, n_genes] <- 0
  return(distance)
}

get_gene_distance_matrix <- function(gb){
  get_range <- function(x){
    start <- gb@transcripts@ranges@start[x]
    width <- gb@transcripts@ranges@width[x]
    end <- start + width
    return(c(start, end))
  }
  # data <- gb@transcripts
  genes <- gb@transcripts$gene_id
  n_genes <- length(genes)
  
  distance <- array(data = NA, dim = c(n_genes, n_genes))
  rownames(distance) <- genes
  colnames(distance) <- genes
  
  for (i in 1:(length(genes)-1)){
    gene_i <- genes[i]
    gene_i_loc <- get_range(i)
    distance[i,i] <- 0
    for (j in (i+1):length(genes)){
      gene_j <- genes[j]
      gene_j_loc <- get_range(j)
      dist <- calculate_gene_distance(gene_i_loc, gene_j_loc)
      distance[i,j] <- dist
      distance[j,i] <- dist
    }
  }
  distance[n_genes, n_genes] <- 0
  return(distance)
}

find_all_neighbors <- function(distance_matrix, idx, cluster){
  neighbors <- union(which(distance_matrix[idx,]), which(distance_matrix[,idx]))
  new_neighbors <- setdiff(neighbors, cluster)
  if (length(new_neighbors) == 0){
    return(cluster)
  }
  cluster <- union(cluster, neighbors)
  for (i in new_neighbors){
    cluster <- find_all_neighbors(distance_matrix, i, cluster)
  }
  return(cluster)
}

cluster_gene_distances <- function(distance_matrix, threshold = 50){
  distance_matrix = (distance_matrix > threshold)
  genes <- rownames(distance_matrix)
  clusters <- c()
  active <- array(data = TRUE, dim = c(length(genes)))
  rownames(active) <- genes
  for (i in 1:nrow(distance_matrix)){
    if (!active[i]){next}
    cluster <- find_all_neighbors(distance_matrix, i, c(genes[i]))
    active[cluster] <- FALSE
    clusters <- c(clusters, list(cluster))
  }
  
  return(clusters)
}

# smu <- readGenBank('~/Documents/jensn lab/mutans_data/genome/S_mutans_UA159.gb')
# pa <- readGenBank('~/Documents/jensn lab/pseudomonas aeruginosa/GCF_000006765.1_ASM676v1_genomic.gbff')
# load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')

smu_gene_order <- get_gene_order_list(smu)
smu_genome_dist_mtx <- get_genome_distance(smu_gene_order)
smu_gene_dist_mtx <- get_gene_distance_matrix(smu)
# pa_genbank_genes <- pa@genes$gene_id
# unannotated_pa_genes <- setdiff(pa_genbank_genes, pao1_active_genes)
# unannotated_functions <- sapply(unannotated_pa_genes, function(x){pa@transcripts$product[grep(x, pa@transcripts$gene_id)]})
# pao_gene_dist_mtx <- get_gene_distance_matrix(pa)
# gene_clusters <- cluster_gene_distances(distance_matrix = pao_gene_dist_mtx)
hclustfunc <- function(x, method = "complete", dmeth = "euclidean") {    
  hclust(dist(x, method = dmeth), method = method)
}

load("~/Documents/jensn lab/Enzyme Coupling/data/mutans_g0_sets.RData")
source('~/GitHub/PathwayMining/falcon_tools.R')

g0_sets <- clean_rxn_names_in_set(g0_sets)
g1_sets <- clean_rxn_names_in_set(mutans_g1_sets)
