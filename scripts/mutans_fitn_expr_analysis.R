# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/data_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')

load("~/Documents/jensn lab/mutans_data/Fitness & Expression/Tn_seq_Counts.RData")
load("~/Documents/jensn lab/mutans_data/Fitness & Expression/Baker_UA159ss_pH7_Glucose_Shock.RData")
load("~/GitHub/PathwayMining/data/mutans_model/mutans_falcon_g1_matrix.RData")
load("~/GitHub/PathwayMining/data/mutans_model/mutans_model.RData")
mutans_falcon <- generate_falcon_model(mutans)
## IDENTIFY GENES OF INTEREST WITH HIGH EXPRESSION CHANGE AND LOW FITNESS CHANGE

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

inf_idx <- which(fitn_fold_change == 'Inf')
nan_idx <- which(fitn_fold_change == 'NaN')

fitn_fold_change[inf_idx] <- '100'
fitn_fold_change <- fitn_fold_change[,-c(nan_idx)]

## combine with Baker paper
expr_fold_change <- cbind(Baker_UA159ss_pH7_Glucose_Shock$`GB acc`, Baker_UA159ss_pH7_Glucose_Shock$`Fold Change`)

fitn_expr_fold_change <- combine_tables(fitn_fold_change, expr_fold_change)
colnames(fitn_expr_fold_change) <- c('Gene_ID', 'Fitness_Fold_Change', 'Expression_Fold_Change')
fitn_expr_fold_change <- fitn_expr_fold_change[-c(which(is.na(fitn_expr_fold_change[,3]))),]
fitn_expr_fold_change <- fitn_expr_fold_change[-c(which(fitn_expr_fold_change[,2] == 'NaN')),]
log2_fitn_expr_fold_change <- matrix(nrow = nrow(fitn_expr_fold_change), ncol = ncol(fitn_expr_fold_change))
log2_fitn_expr_fold_change[,1] <- fitn_expr_fold_change[,1]
log2_fitn_expr_fold_change[,2] <- log2(as.numeric(as.character(fitn_expr_fold_change[,2])))
log2_fitn_expr_fold_change[,3] <- log2(as.numeric(as.character(fitn_expr_fold_change[,3])))

genes_of_interest <- identify_fitn_expr_relation(log2_fitn_expr_fold_change) #, fitn_thresh = 0.01, expr_thresh = 0.25)
gene_idx <- which(log2_fitn_expr_fold_change %in% genes_of_interest)
# gene check

available_genes <- matrix(data = FALSE, nrow = 1, ncol = length(log2_fitn_expr_fold_change[,1]))
for(i in 1:length(log2_fitn_expr_fold_change[,1])){
  # avail <- log2_fitn_expr_fold_change[i,1] %in% mutans@allGenes
  avail <- any(grepl(log2_fitn_expr_fold_change[i,1], mutans@allGenes))
  # if (avail){print('y')}
  available_genes[i] <- avail
}

available_genes2 <- matrix(data = FALSE, nrow = 1, ncol = length(mutans@allGenes))
for(i in 1:length(mutans@allGenes)){
  # avail <- log2_fitn_expr_fold_change[i,1] %in% mutans@allGenes
  avail <- any(grepl(mutans@allGenes[i], log2_fitn_expr_fold_change[,1]))
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

# amount of overlap
length(which(available_genes))
length(which(available_genes2))
length(which(available_genes3))

# g0
# load g0 mtx
load("~/GitHub/PathwayMining/data/mutans_model/mutans_r0_coupling_mtx.RData")

#mutans_g0_matrix <- isolate_gene_matrix(mutans_r0_coupling_mtx)
#clean_mutans_g0_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_g0_matrix)))[[1]])

g0_gene_map <- map_elements_to_set(genes_of_interest, clean_mutans_g0_set)

# gene positions
g0_genes <- g0_gene_map[which(!is.na(g0_gene_map))]
g0_dupl <- duplicated(g0_genes)
print('g0 sets containing genes of interest')
print(g0_genes)
print('sets containing multiple genes of interest')
print(genes[g0_dupl])

# g1
#mutans_g1_matrix <- isolate_gene_matrix(mutans_falcon_g1_matrix)
#clean_mutans_g1_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_g1_matrix)))[[1]])

gene_map <- map_elements_to_set(genes_of_interest, clean_mutans_g1_set)

# gene positions
genes <- gene_map[which(!is.na(gene_map))]
dupl <- duplicated(genes)
print('g1 sets containing genes of interest')
print(genes)
print('sets containing multiple genes of interest')
print(genes[dupl])
# duplicated sets: [1] 237 122 120  71 170 236 228 228 165 165 165 122  43

# plot(fitn_expr_fold_change[,2], fitn_expr_fold_change[,3], main = 'Fitness and Expression Change', xlab = 'Fitness Fold Change', ylab = 'Expression Fold Change')
plot(log2_fitn_expr_fold_change[,2], log2_fitn_expr_fold_change[,3], main = 'Fitness and Expression Change', xlab = 'Log2 Fitness Fold Change', ylab = 'Log2 Expression Fold Change')
# plot(1:nrow(log2_fitn_expr_fold_change), log2_fitn_expr_fold_change[,2], main = 'Fitness Fold Change Under Mild Nutrient Stress', xlab = 'Gene Index', ylab = 'Log2 Fold Change') # fitness
plot(1:nrow(log2_fitn_expr_fold_change), log2_fitn_expr_fold_change[,3], main = 'Expression Fold Change Under Acid Stress', xlab = 'Gene Index', ylab = 'Log2 Fold Change') # expr
# plot(1:nrow(fitn_expr_fold_change), fitn_expr_fold_change[,3]) # expression

# non_zero <- which(fitn_expr_fold_change[,3] > 0)

# composition_size <- matrix(data = 0, nrow = length(composition), ncol = 1)
# for (i in 1:length(composition)){
#   composition_size[i] <- length(composition[[i]])
# }

mutans_g1_set_composition <- find_set_list_composition(mutans_falcon_g1_sets, mutans_falcon_g0_sets)
composition_size <- rowSums(mutans_g1_set_composition)
composing_sets <- which(composition_size > 1)
sets_of_interest <- intersect(composing_sets, gene_map)
gene_idxs_of_special_interest <- which(gene_map %in% sets_of_interest)


## FIND SYNTHETIC LETHAL PAIRS

mutans_obj <- get_mutans_model_w_obj()

#singleGeneDels <- oneGeneDel(mutans_obj, geneList = mutans_obj@allGenes)
#doubleGeneDels <- doubleGeneDel(mutans_obj, allComb = TRUE)

#effect <- which(dblGeneDels@hasEffect)
pairs <- doubleGeneDels@dels

print('lethal single dels:')
lethal_single_dels <- which(near(0, single_gene_ko_max_flux))
print(lethal_single_dels)
print('lethal double dels:')
lethal_double_dels <- which(near(0, double_gene_ko_max_flux))
print(lethal_double_dels)

captured_pairs <- matrix(data = FALSE, nrow = length(lethal_double_dels), ncol = 1)
for (i in lethal_double_dels){
  set1 <- get_set_idx(pairs[i, 1], clean_mutans_g1_set)
  set2 <- get_set_idx(pairs[i, 2], clean_mutans_g1_set)
  #if (length(set1) > 0 & length(set2) > 0){
  if (set1 == set2){captured_pairs[i] <- TRUE}
  #}
}

print('synth lethality:')
print(length(which(captured_pairs)))
print('out of')
print(length(lethal_double_dels))

# check overlap between model simulation and observed data
potential_synth_lethals <- pairs[lethal_double_dels,]
potential_synth_lethals <- unique(c(potential_synth_lethals[,1], potential_synth_lethals[,2]))

## test lethality in falcon model

pairs <- doubleGeneDels@dels
new_pairs <- matrix(nrow = nrow(pairs), ncol = ncol(pairs))
for (i in 1:nrow(pairs)){
  for (j in 1:ncol(pairs)){
    new_pairs[i,j] <- paste('Ex_a_', pairs[i,j], sep = '')
  }
}

print('lethal single dels:')
lethal_single_dels <- which(near(0, single_gene_ko_max_flux))
print(lethal_single_dels)
print('lethal double dels:')
lethal_double_dels <- which(near(0, double_gene_ko_max_flux))
print(lethal_double_dels)

lethality <- matrix(data = FALSE, nrow = length(lethal_double_dels), ncol = 1)
for (i in 1:length(lethal_double_dels)){
  lethality[i] <- GRB_maximize(mutans_obj_falcon_model, 477, new_pairs[lethal_double_dels[i],])
}

print('synth lethality:')
print(length(which(captured_pairs)))
print('out of')
print(length(lethal_double_dels))

## look for synth lethals in g0 sets

g0_gene_set_location <- matrix(data = FALSE, nrow = length(lethal_double_dels), ncol = 1)
for (i in 1:length(lethal_double_dels)){
  location1 <- get_set_idx(new_pairs[lethal_double_dels[i],1], mutans_falcon_g0_sets)
  location2 <- get_set_idx(new_pairs[lethal_double_dels[i],2], mutans_falcon_g0_sets)
  if (is.null(location1) | is.null(location2)){next}
  g0_gene_set_location[i] <- (location1 == location2)
}

print(which(g0_gene_set_location))

## repeat
g0_gene_set_location <- matrix(data = FALSE, nrow = length(lethal_double_dels), ncol = 1)
for (i in 1:length(lethal_double_dels)){
  location1 <- get_set_idx(pairs[lethal_double_dels[i],1], clean_mutans_g0_sets)
  location2 <- get_set_idx(pairs[lethal_double_dels[i],2], clean_mutans_g0_sets)
  if (is.null(location1) | is.null(location2)){next}
  g0_gene_set_location[i] <- (location1 == location2)
}

print(which(g0_gene_set_location))


## look for synth lethals in g1 sets

g1_gene_set_location <- matrix(data = FALSE, nrow = length(lethal_double_dels), ncol = 1)
for (i in 1:length(lethal_double_dels)){
  location1 <- get_set_idx(new_pairs[lethal_double_dels[i],1], mutans_falcon_g1_sets)
  location2 <- get_set_idx(new_pairs[lethal_double_dels[i],2], mutans_falcon_g1_sets)
  #if (is.null(location1) | is.null(location2)){next}
  g1_gene_set_location[i] <- (location1 == location2)
}

print(which(g1_gene_set_location))

captured_pairs <- which(g1_gene_set_location)

## repeat
g1_gene_set_location <- matrix(data = FALSE, nrow = length(lethal_double_dels), ncol = 1)
for (i in 1:length(lethal_double_dels)){
  location1 <- get_set_idx(pairs[lethal_double_dels[i],1], clean_mutans_g1_set)
  location2 <- get_set_idx(pairs[lethal_double_dels[i],2], clean_mutans_g1_set)
  #if (is.null(location1) | is.null(location2)){next}
  g1_gene_set_location[i] <- (location1 == location2)
}

print(which(g1_gene_set_location))


# metabolite distances

# get mets and set up active

library(readr)
mutans_model_met <- read_delim("~/GitHub/PathwayMining/data/mutans_model/mutans_model_met.csv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
# View(mutans_model_met)

inactive_mets <- c('H+','CTP','Adenine','BIOT','H2O','Acetyl-CoA','ocdca','Adenosine','ATP','UMP','H2O2','Mn2+','ADP','trdox','dADP','Ca2+','Phosphate','trdrd','dGDP','Cl-','PPi','O2','dTDP','Co2+','CO2','Isopentenyldiphosphate','dCTP','H2S','NADP','UDP-N-acetylglucosamine','dCDP','dTMP','NADPH','GTP','dUDP','Na+','CoA','S-Adenosyl-L-methionine','Heme','Cd2+','ACP','UTP','Malonyl-CoA','dCMP','NAD','GDP','Zn2+','H2S2O3','NADH','Niacin','Cu2+','Hg2+','CMP','dATP','H2CO3','Pb','Malonyl-acyl-carrierprotein-','GMP','Mg','Biomass','NH3','K+','dUMP','dGMP','AMP','dGTP','dUTP','O2-','UDP','Fe2+','fe3')
inactive_met_idxs <- c()

for (met in inactive_mets){
  idx <- which(mutans_falcon@met_name == met)
  if (length(idx) < 1){next}
  inactive_met_idxs <- c(idx, inactive_met_idxs)
}

a_pairs <- matrix(nrow = nrow(pairs), ncol = ncol(pairs)) # lethal pairs represented by corresponding metabolites (a_SMU...)
for (i in 1:nrow(pairs)){
  for (j in 1:ncol(pairs)){
    a_pairs[i,j] <- paste('a_', pairs[i,j], sep = '')
  }
}
lethal_a_pairs <- a_pairs[lethal_double_dels,]

active_mets <- matrix(data = TRUE, nrow = 1, ncol = length(mutans_falcon@met_id)) # need to substitute falcon model
# colnames(active_mets) <- mutans_falcon@met_id
active_mets[inactive_met_idxs] <- FALSE

for (i in 1:nrow(lethal_a_pairs)){
  if (i %in% captured_pairs){next}
  gene_1 <- which(mutans_falcon@met_id == lethal_a_pairs[i,1])
  gene_2 <- which(mutans_falcon@met_id == lethal_a_pairs[i,2])
  
  output <- get_path_mtx_between_reactions(mutans_falcon@S, gene_1, gene_2, active_mets = active_mets)
  path <- trace_path_mtx_between_reactions(output, gene_2, gene_1)
  print(length(path))
}

for (i in captured_pairs){
  gene_1 <- which(mutans_falcon@met_id == lethal_a_pairs[i,1])
  gene_2 <- which(mutans_falcon@met_id == lethal_a_pairs[i,2])
  # print(gene_1)
  # print(gene_2)
  output <- get_path_mtx_between_reactions(mutans_falcon@S, gene_1, gene_2, active_mets = active_mets)
  path <- trace_path_mtx_between_reactions(output, gene_2, gene_1)
  print(length(path))
}
