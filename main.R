# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/data_tools.R')

## get fold change of fitness from Burne paper

BHI_1_total <- sum(as.numeric(Tn_seq_Counts$`BHI 1`))
BHI_2_total <- sum(as.numeric(Tn_seq_Counts$`BHI 2`))
B_A_1_total <- sum(as.numeric(Tn_seq_Counts$`Blood agar 1`))
B_A_2_total <- sum(as.numeric(Tn_seq_Counts$`Blood agar 2`))

BHI_1_freq <- as.numeric(Tn_seq_Counts$`BHI 1`)/BHI_1_total
BHI_2_freq <- as.numeric(Tn_seq_Counts$`BHI 2`)/BHI_2_total
B_A_1_freq <- as.numeric(Tn_seq_Counts$`Blood agar 1`)/B_A_1_total
B_A_2_freq <- as.numeric(Tn_seq_Counts$`Blood agar 2`)/B_A_2_total

BHI_avg <- rowMeans(cbind(BHI_1_freq, BHI_2_freq))
B_A_avg <- rowMeans(cbind(B_A_1_freq, B_A_2_freq))

fitn_fold_change <- BHI_avg/B_A_avg

fitn_fold_change <- cbind(Tn_seq_Counts$Gene_ID, fitn_fold_change)

## combine with Baker paper
expr_fold_change <- cbind(Baker_UA159ss_pH7_Glucose_Shock$`GB acc`, Baker_UA159ss_pH7_Glucose_Shock$`Fold Change`)
fitn_expr_fold_change <- combine_tables(fitn_fold_change, expr_fold_change)
colnames(fitn_expr_fold_change) <- c('Gene_ID', 'Fitness_Fold_Change', 'Expression_Fold_Change')
fitn_expr_fold_change[,2] <- log2(as.numeric(as.character(fitn_expr_fold_change[,2])))
fitn_expr_fold_change[,3] <- log2(as.numeric(as.character(fitn_expr_fold_change[,3])))

genes_of_interest <- identify_fitn_expr_relation(fitn_expr_fold_change) #, fitn_thresh = 0.01, expr_thresh = 0.25)

# gene check

available_genes <- matrix(data = FALSE, nrow = 1, ncol = length(fitn_expr_fold_change[,1]))
for(i in 1:length(fitn_expr_fold_change[,1])){
  # avail <- fitn_expr_fold_change[i,1] %in% mutans@allGenes
  avail <- any(grepl(fitn_expr_fold_change[i,1], mutans@allGenes))
  # if (avail){print('y')}
  available_genes[i] <- avail
}

available_genes2 <- matrix(data = FALSE, nrow = 1, ncol = length(mutans@allGenes))
for(i in 1:length(mutans@allGenes)){
  # avail <- fitn_expr_fold_change[i,1] %in% mutans@allGenes
  avail <- any(grepl(mutans@allGenes[i], fitn_expr_fold_change[,1]))
  # if (avail){print('y')}
  available_genes2[i] <- avail
}

available_genes3 <- matrix(data = FALSE, nrow = 1, ncol = length(genes_of_interest))
for(i in 1:length(genes_of_interest)){
  # avail <- genes_of_interest[i] %in% mutans@allGenes
  avail <- any(grepl(genes_of_interest[i], mutans@allGenes))
  # if (avail){print('y')}
  available_genes3[i] <- avail
}

gene_map <- map_elements_to_set(genes_of_interest, clean_mutans_g1_set)
