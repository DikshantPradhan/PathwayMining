# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
# source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')
# source('~/GitHub/PathwayMining/load_mod.R')


r0_rxn_dist <- get_size_distribution(yeast_og_set_list)
r1_rxn_dist <- get_size_distribution(r1_set_list)
r0_gene_dist <- get_size_distribution(r0_gene_set_list)
r1_gene_dist <- get_size_distribution(r1_gene_set_list)

r0_recurring_genes <- find_recurring_genes_in_set_list(r0_gene_set_list)
r1_recurring_genes <- find_recurring_genes_in_set_list(r1_gene_set_list)

e_thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)

gene_pairs_00 <- check_for_enrichment(all_gene_pairs, gi_matrix$e, 0.0)

y_open_gi_matrix <- build_gi_matrix_from_pairs(gene_pairs_00)

gene_pairs_01 <- check_for_enrichment(gene_pairs_00[,1:2], gi_matrix$e, 0.1)
r0_gene_pairs_01 <- recurring_pairs(r0_gene_pairs, gene_pairs_01)
r1_gene_pairs_01 <- recurring_pairs(r1_gene_pairs, gene_pairs_01)

gene_pairs_02 <- check_for_enrichment(gene_pairs_01[,1:2], gi_matrix$e, 0.2)
r0_gene_pairs_02 <- recurring_pairs(r0_gene_pairs, gene_pairs_02)
r1_gene_pairs_02 <- recurring_pairs(r1_gene_pairs, gene_pairs_02)

gene_pairs_03 <- check_for_enrichment(gene_pairs_02[,1:2], gi_matrix$e, 0.3)
r0_gene_pairs_03 <- recurring_pairs(r0_gene_pairs, gene_pairs_03)
r1_gene_pairs_03 <- recurring_pairs(r1_gene_pairs, gene_pairs_03)

gene_pairs_04 <- check_for_enrichment(gene_pairs_03[,1:2], gi_matrix$e, 0.4)
r0_gene_pairs_04 <- recurring_pairs(r0_gene_pairs, gene_pairs_04)
r1_gene_pairs_04 <- recurring_pairs(r1_gene_pairs, gene_pairs_04)

gene_pairs_05 <- check_for_enrichment(gene_pairs_04[,1:2], gi_matrix$e, 0.5)
r0_gene_pairs_05 <- recurring_pairs(r0_gene_pairs, gene_pairs_05)
r1_gene_pairs_05 <- recurring_pairs(r1_gene_pairs, gene_pairs_05)


print(paste('0.1: ', 'r0-', length(r0_gene_pairs_01[,1])/length(gene_pairs_01[,1]),
  'r1-', length(r1_gene_pairs_01[,1])/length(gene_pairs_01[,1])))
print(paste('0.2: ', 'r0-', length(r0_gene_pairs_02[,1])/length(gene_pairs_02[,1]),
  'r1-', length(r1_gene_pairs_02[,1])/length(gene_pairs_02[,1])))
print(paste('0.3: ', 'r0-', length(r0_gene_pairs_03[,1])/length(gene_pairs_03[,1]),
  'r1-', length(r1_gene_pairs_03[,1])/length(gene_pairs_03[,1])))
print(paste('0.4: ', 'r0-', length(r0_gene_pairs_04[,1])/length(gene_pairs_04[,1]),
  'r1-', length(r1_gene_pairs_04[,1])/length(gene_pairs_04[,1])))
print(paste('0.5: ', 'r0-', length(r0_gene_pairs_05[,1])/length(gene_pairs_05[,1]),
  'r1-', length(r1_gene_pairs_05[,1])/length(gene_pairs_05[,1])))
