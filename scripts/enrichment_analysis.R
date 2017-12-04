# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
# source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
# source('~/GitHub/PathwayMining/load_mod.R')

## MAIN ENRICHMENT (making edits for gene set list from falcon model)

r0_gene_set_list <- clean_rxn_names_in_set(yeast_falcon_r0_set_list)
r1_gene_set_list <- clean_rxn_names_in_set(yeast_falcon_r1_set_list)

all_gene_pairs <- return_pairs_from_set(unlist(r0_gene_set_list))

#load gi matrix
load('gi_matrix.RData')

# r0_rxn_dist <- get_size_distribution(r0_set_list)
# r1_rxn_dist <- get_size_distribution(r1_set_list)
r0_gene_dist <- get_size_distribution(r0_gene_set_list)
r1_gene_dist <- get_size_distribution(r1_gene_set_list)

# r0_recurring_genes <- find_recurring_genes_in_set_list(r0_gene_set_list)
# r1_recurring_genes <- find_recurring_genes_in_set_list(r1_gene_set_list)

e_thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)

gene_pairs_00 <- check_for_enrichment(all_gene_pairs, gi_matrix$e, 0.0)

y_open_gi_matrix <- build_gi_matrix_from_pairs(unlist(r0_gene_set_list), gene_pairs_00)

save(gene_pairs_00, y_open_gi_matrix, file = 'remote_1.RData')

e_matrix_01 <- matrix(as.numeric(abs(y_open_gi_matrix) >= 0.1), nrow = nrow(y_open_gi_matrix), ncol = ncol(y_open_gi_matrix))
rownames(e_matrix_01) <- rownames(y_open_gi_matrix)
colnames(e_matrix_01) <- colnames(y_open_gi_matrix)
e_matrix_02 <- matrix(as.numeric(abs(y_open_gi_matrix) >= 0.2), nrow = nrow(y_open_gi_matrix), ncol = ncol(y_open_gi_matrix))
rownames(e_matrix_02) <- rownames(y_open_gi_matrix)
colnames(e_matrix_02) <- colnames(y_open_gi_matrix)
e_matrix_03 <- matrix(as.numeric(abs(y_open_gi_matrix) >= 0.3), nrow = nrow(y_open_gi_matrix), ncol = ncol(y_open_gi_matrix))
rownames(e_matrix_03) <- rownames(y_open_gi_matrix)
colnames(e_matrix_03) <- colnames(y_open_gi_matrix)
e_matrix_04 <- matrix(as.numeric(abs(y_open_gi_matrix) >= 0.4), nrow = nrow(y_open_gi_matrix), ncol = ncol(y_open_gi_matrix))
rownames(e_matrix_04) <- rownames(y_open_gi_matrix)
colnames(e_matrix_04) <- colnames(y_open_gi_matrix)
e_matrix_05 <- matrix(as.numeric(abs(y_open_gi_matrix) >= 0.5), nrow = nrow(y_open_gi_matrix), ncol = ncol(y_open_gi_matrix))
rownames(e_matrix_05) <- rownames(y_open_gi_matrix)
colnames(e_matrix_05) <- colnames(y_open_gi_matrix)

e_matrix_07 <- matrix(as.numeric(abs(y_open_gi_matrix) >= 0.7), nrow = nrow(y_open_gi_matrix), ncol = ncol(y_open_gi_matrix))
rownames(e_matrix_07) <- rownames(y_open_gi_matrix)
colnames(e_matrix_07) <- colnames(y_open_gi_matrix)

save(e_matrix_01, e_matrix_02, e_matrix_03, e_matrix_04, e_matrix_05, file = 'remote_2.RData')

r0_gi_01 <- get_num_interactions(r0_gene_set_list, e_matrix_01)
r0_gi_02 <- get_num_interactions(r0_gene_set_list, e_matrix_02)
r0_gi_03 <- get_num_interactions(r0_gene_set_list, e_matrix_03)
r0_gi_04 <- get_num_interactions(r0_gene_set_list, e_matrix_04)
r0_gi_05 <- get_num_interactions(r0_gene_set_list, e_matrix_05)

r1_gi_01 <- get_num_interactions(r1_gene_set_list, e_matrix_01)
r1_gi_02 <- get_num_interactions(r1_gene_set_list, e_matrix_02)
r1_gi_03 <- get_num_interactions(r1_gene_set_list, e_matrix_03)
r1_gi_04 <- get_num_interactions(r1_gene_set_list, e_matrix_04)
r1_gi_05 <- get_num_interactions(r1_gene_set_list, e_matrix_05)

#gene_pairs_01 <- check_for_enrichment(gene_pairs_00[,1:2], y_open_gi_matrix, 0.1)
#gene_pairs_02 <- check_for_enrichment(gene_pairs_01[,1:2], y_open_gi_matrix, 0.2)
#gene_pairs_03 <- check_for_enrichment(gene_pairs_02[,1:2], y_open_gi_matrix, 0.3)
#gene_pairs_04 <- check_for_enrichment(gene_pairs_03[,1:2], y_open_gi_matrix, 0.4)
#gene_pairs_05 <- check_for_enrichment(gene_pairs_04[,1:2], y_open_gi_matrix, 0.5)

#r0_gene_pairs_01 <- recurring_pairs(r0_gene_pairs, gene_pairs_01)
#r0_gene_pairs_02 <- recurring_pairs(r0_gene_pairs_01, gene_pairs_02)
#r0_gene_pairs_03 <- recurring_pairs(r0_gene_pairs_02, gene_pairs_03)
#r0_gene_pairs_04 <- recurring_pairs(r0_gene_pairs_03, gene_pairs_04)
#r0_gene_pairs_05 <- recurring_pairs(r0_gene_pairs_04, gene_pairs_05)

#r1_gene_pairs_01 <- recurring_pairs(r1_gene_pairs, gene_pairs_01)
#r1_gene_pairs_02 <- recurring_pairs(r1_gene_pairs_01, gene_pairs_02)
#r1_gene_pairs_03 <- recurring_pairs(r1_gene_pairs_02, gene_pairs_03)
#r1_gene_pairs_04 <- recurring_pairs(r1_gene_pairs_03, gene_pairs_04)
#r1_gene_pairs_05 <- recurring_pairs(r1_gene_pairs_04, gene_pairs_05)

#print(paste('0.1: ', 'r0-', length(r0_gene_pairs_01[,1])/length(gene_pairs_01[,1]),
#  'r1-', length(r1_gene_pairs_01[,1])/length(gene_pairs_01[,1])))
#print(paste('0.2: ', 'r0-', length(r0_gene_pairs_02[,1])/length(gene_pairs_02[,1]),
#  'r1-', length(r1_gene_pairs_02[,1])/length(gene_pairs_02[,1])))
#print(paste('0.3: ', 'r0-', length(r0_gene_pairs_03[,1])/length(gene_pairs_03[,1]),
#  'r1-', length(r1_gene_pairs_03[,1])/length(gene_pairs_03[,1])))
#print(paste('0.4: ', 'r0-', length(r0_gene_pairs_04[,1])/length(gene_pairs_04[,1]),
#  'r1-', length(r1_gene_pairs_04[,1])/length(gene_pairs_04[,1])))
#print(paste('0.5: ', 'r0-', length(r0_gene_pairs_05[,1])/length(gene_pairs_05[,1]),
#  'r1-', length(r1_gene_pairs_05[,1])/length(gene_pairs_05[,1])))


## SAMPLING FOR COMPARISON (FULLY RANDOM BASED ON SIZE)

r0_gene_elements <- unlist(r0_gene_set_list)
r1_gene_elements <- unlist(r1_gene_set_list)

r0_samples <- sample_multiple_sets_to_distribution(1000, r0_gene_elements, r0_gene_dist, replacement = FALSE)
r1_samples <- sample_multiple_sets_to_distribution(1000, r1_gene_elements, r1_gene_dist, replacement = FALSE)

r0_sampled_pair_ct <- matrix(0, nrow = length(r0_samples), ncol = 5)
r1_sampled_pair_ct <- matrix(0, nrow = length(r1_samples), ncol = 5)

for (i in 1:length(r0_samples)){
  #r0_sampled_gene_pairs <- remove_duplicate_pairs(return_pairs_from_set_list(r0_samples[[i]]))

  #r0_sampled_gene_pairs_01 <- recurring_pairs(r0_sampled_gene_pairs, gene_pairs_01)
  #r0_sampled_gene_pairs_02 <- recurring_pairs(r0_sampled_gene_pairs_01, gene_pairs_02)
  #r0_sampled_gene_pairs_03 <- recurring_pairs(r0_sampled_gene_pairs_02, gene_pairs_03)
  #r0_sampled_gene_pairs_04 <- recurring_pairs(r0_sampled_gene_pairs_03, gene_pairs_04)
  #r0_sampled_gene_pairs_05 <- recurring_pairs(r0_sampled_gene_pairs_04, gene_pairs_05)

  r0_sampled_pair_ct[i,1] <- get_num_interactions(r0_samples[[i]], e_matrix_01) #length(r0_sampled_gene_pairs_01[,1])
  r0_sampled_pair_ct[i,2] <- get_num_interactions(r0_samples[[i]], e_matrix_02) #length(r0_sampled_gene_pairs_02[,1])
  r0_sampled_pair_ct[i,3] <- get_num_interactions(r0_samples[[i]], e_matrix_03) #length(r0_sampled_gene_pairs_03[,1])
  r0_sampled_pair_ct[i,4] <- get_num_interactions(r0_samples[[i]], e_matrix_04) #length(r0_sampled_gene_pairs_04[,1])
  r0_sampled_pair_ct[i,5] <- get_num_interactions(r0_samples[[i]], e_matrix_05) #length(r0_sampled_gene_pairs_05[,1])
}

for (i in 1:length(r1_samples)){
  #r1_sampled_gene_pairs <- remove_duplicate_pairs(return_pairs_from_set_list(r1_samples[[i]]))

  #r1_sampled_gene_pairs_01 <- recurring_pairs(r1_sampled_gene_pairs, gene_pairs_01)
  #r1_sampled_gene_pairs_02 <- recurring_pairs(r1_sampled_gene_pairs_01, gene_pairs_02)
  #r1_sampled_gene_pairs_03 <- recurring_pairs(r1_sampled_gene_pairs_02, gene_pairs_03)
  #r1_sampled_gene_pairs_04 <- recurring_pairs(r1_sampled_gene_pairs_03, gene_pairs_04)
  #r1_sampled_gene_pairs_05 <- recurring_pairs(r1_sampled_gene_pairs_04, gene_pairs_05)

  r1_sampled_pair_ct[i,1] <- get_num_interactions(r1_samples[[i]], e_matrix_01) #length(r1_sampled_gene_pairs_01[,1])
  r1_sampled_pair_ct[i,2] <- get_num_interactions(r1_samples[[i]], e_matrix_02) #length(r1_sampled_gene_pairs_02[,1])
  r1_sampled_pair_ct[i,3] <- get_num_interactions(r1_samples[[i]], e_matrix_03) #length(r1_sampled_gene_pairs_03[,1])
  r1_sampled_pair_ct[i,4] <- get_num_interactions(r1_samples[[i]], e_matrix_04) #length(r1_sampled_gene_pairs_04[,1])
  r1_sampled_pair_ct[i,5] <- get_num_interactions(r1_samples[[i]], e_matrix_05) #length(r1_sampled_gene_pairs_05[,1])
}

print('R0')
print(paste(mean(r0_sampled_pair_ct[,1]), ', ', mean(r0_sampled_pair_ct[,2]), ', ', mean(r0_sampled_pair_ct[,3]), ', ', mean(r0_sampled_pair_ct[,4]), ', ', mean(r0_sampled_pair_ct[,5])))

print('R1')
print(paste(mean(r1_sampled_pair_ct[,1]), ', ', mean(r1_sampled_pair_ct[,2]), ', ', mean(r1_sampled_pair_ct[,3]), ', ', mean(r1_sampled_pair_ct[,4]), ', ', mean(r1_sampled_pair_ct[,5])))

print('R0')
print(paste(sd(r0_sampled_pair_ct[,1]), ', ', sd(r0_sampled_pair_ct[,2]), ', ', sd(r0_sampled_pair_ct[,3]), ', ', sd(r0_sampled_pair_ct[,4]), ', ', sd(r0_sampled_pair_ct[,5])))

print('R1')
print(paste(sd(r1_sampled_pair_ct[,1]), ', ', sd(r1_sampled_pair_ct[,2]), ', ', sd(r1_sampled_pair_ct[,3]), ', ', sd(r1_sampled_pair_ct[,4]), ', ', sd(r1_sampled_pair_ct[,5])))

## ALTERNATE SAMPLING FOR COMPARISON (RANDOM BASED ON COMPOSITION AND SIZES OF SETS BEING COMPOSED)

r1_composition <- find_set_list_composition(r1_gene_set_list, r0_gene_set_list) # _set_list -> _gene_set_list
r1_size_composition <- get_composition_size_distribution(r0_gene_set_list, r1_composition) # _set_list -> _gene_set_list
r0_size_class <- classify_sets_by_size(r0_gene_set_list)

r1_rxn_samples <- sample_multiple_sets_to_composition(1000, r0_size_class, r1_size_composition)

r1_sampled_pair_ct <- matrix(0, nrow = length(r1_samples), ncol = 5)

for (i in 1:length(r1_rxn_samples)){
  print(i)

  r1_gene_sample <- gene_set_from_composition(r0_gene_set_list, r1_rxn_samples[[i]])

  r1_sampled_pair_ct[i,1] <- get_num_interactions(r1_gene_sample, e_matrix_01) #length(r1_sampled_gene_pairs_01[,1])
  r1_sampled_pair_ct[i,2] <- get_num_interactions(r1_gene_sample, e_matrix_02) #length(r1_sampled_gene_pairs_02[,1])
  r1_sampled_pair_ct[i,3] <- get_num_interactions(r1_gene_sample, e_matrix_03) #length(r1_sampled_gene_pairs_03[,1])
  r1_sampled_pair_ct[i,4] <- get_num_interactions(r1_gene_sample, e_matrix_04) #length(r1_sampled_gene_pairs_04[,1])
  r1_sampled_pair_ct[i,5] <- get_num_interactions(r1_gene_sample, e_matrix_05) #length(r1_sampled_gene_pairs_05[,1])
}

print('R1')
print(paste(mean(r1_sampled_pair_ct[,1]), ', ', mean(r1_sampled_pair_ct[,2]), ', ', mean(r1_sampled_pair_ct[,3]), ', ', mean(r1_sampled_pair_ct[,4]), ', ', mean(r1_sampled_pair_ct[,5])))

print('R1')
print(paste(sd(r1_sampled_pair_ct[,1]), ', ', sd(r1_sampled_pair_ct[,2]), ', ', sd(r1_sampled_pair_ct[,3]), ', ', sd(r1_sampled_pair_ct[,4]), ', ', sd(r1_sampled_pair_ct[,5])))

a <- paste(mean(r1_sampled_pair_ct[,1]), ', ', mean(r1_sampled_pair_ct[,2]), ', ', mean(r1_sampled_pair_ct[,3]), ', ', mean(r1_sampled_pair_ct[,4]), ', ', mean(r1_sampled_pair_ct[,5]))
b <- paste(sd(r1_sampled_pair_ct[,1]), ', ', sd(r1_sampled_pair_ct[,2]), ', ', sd(r1_sampled_pair_ct[,3]), ', ', sd(r1_sampled_pair_ct[,4]), ', ', sd(r1_sampled_pair_ct[,5]))

save(a, b, file = 'remote_3.RData')

print('FIN')
