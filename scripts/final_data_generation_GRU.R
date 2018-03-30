source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/gene_tools.R')

# LOAD MODELS
load('~/GitHub/PathwayMining/data/mutans_model/mutans_model.RData')
mutans_sybil <- mutans
load('~/GitHub/PathwayMining/data/yeast_model/Maranas_model/maranas_model_lipid_exch.RData')
yeast_sybil <- yeast_model

mutans <- GRB_mutans_model()
yeast <- GRB_yeast_model()

# GENERATE FALCON MODELS
mutans_falcon <- GRB_mutans_falcon_model()
yeast_falcon <- GRB_yeast_falcon_model()

# MODEL STATISTICS
mutans_n <- mutans$get_sizes()$NumVars
mutans_vars <- mutans$get_names()$VarName
yeast_n <- yeast$get_sizes()$NumVars
yeast_vars <- yeast$get_names()$VarName

mutans_gpr <- generate_gpr(mutans_sybil)
yeast_gpr <- generate_gpr(yeast_sybil)

save(mutans_gpr, file = 'final_paper_data/mutans_gpr.RData')
save(yeast_gpr, file = 'final_paper_data/yeast_gpr.RData')

# R0 SETS
mutans_r0_sets <- GRB_generate_set_list(mutans)
yeast_r0_sets <- GRB_generate_set_list(yeast)

save(mutans_r0_sets, file = 'final_paper_data/mutans_r0_sets.RData')
save(yeast_r0_sets, file = 'final_paper_data/yeast_r0_sets.RData')

# G(R0) SETS
mutans_gr0_sets <- gene_set_from_rxn_set(mutans_gpr, mutans_r0_sets)
yeast_gr0_sets <- gene_set_from_rxn_set(yeast_gpr, yeast_r0_sets)

save(mutans_r0_sets, file = 'final_paper_data/mutans_gr0_sets.RData')
save(yeast_r0_sets, file = 'final_paper_data/yeast_gr0_sets.RData')

# G0 SETS
mutans_falcon_r0_sets <- GRB_generate_set_list(mutans_falcon)
yeast_falcon_r0_sets <- GRB_generate_set_list(yeast_falcon)

save(mutans_r0_sets, file = 'final_paper_data/mutans_g0_sets.RData')
save(yeast_r0_sets, file = 'final_paper_data/yeast_g0_sets.RData')

# R1 SETS
mutans_coupling_array <- GRB_generate_set_lists_array(mutans, 1:mutans_n, compare_known_r0_sets = TRUE, optimize_suppr=TRUE)
mutans_r1_matrix <- coupling_matrix_from_array(mutans_coupling_array)
mutans_r1_matrix <- (mutans_r1_matrix > 0)
mutans_r1_sets <- list(get_list_of_sets(return_couples(mutans_r1_matrix)))

save(mutans_r1_sets, file = 'final_paper_data/mutans_r1_sets.RData')

yeast_coupling_array <- GRB_generate_set_lists_array(yeast, 1:yeast_n, compare_known_r0_sets = TRUE, optimize_suppr=TRUE)
yeast_r1_matrix <- coupling_matrix_from_array(yeast_coupling_array)
yeast_r1_matrix <- (yeast_r1_matrix > 0)
yeast_r1_sets <- list(get_list_of_sets(return_couples(yeast_r1_matrix)))

save(yeast_r1_sets, file = 'final_paper_data/yeast_r1_sets.RData')

# G(R1) SETS
mutans_gr1_sets <- gene_set_from_rxn_set(mutans_gpr, mutans_r1_sets)
yeast_gr1_sets <- gene_set_from_rxn_set(yeast_gpr, yeast_r1_sets)

save(mutans_gr1_sets, file = 'final_paper_data/mutans_gr1_sets.RData')
save(yeast_gr1_sets, file = 'final_paper_data/yeast_gr1_sets.RData')

# G1 SETS

# SET STATISTICS
mutans_r0_set_size <- get_size_list(mutans_r0_sets)
mutans_r0_set_size_dist <- get_size_distribution(mutans_r0_sets)
yeast_r0_set_size <- get_size_list(yeast_r0_sets)
yeast_r0_set_size_dist <- get_size_distribution(yeast_r0_sets)

mutans_g0_set_size <- get_size_list(mutans_g0_sets)
mutans_g0_set_size_dist <- get_size_distribution(mutans_g0_sets)
yeast_g0_set_size <- get_size_list(yeast_g0_sets)
yeast_g0_set_size_dist <- get_size_distribution(yeast_g0_sets)

mutans_gr0_set_size <- get_size_list(mutans_gr0_sets)
mutans_gr0_set_size_dist <- get_size_distribution(mutans_gr0_sets)
yeast_gr0_set_size <- get_size_list(yeast_gr0_sets)
yeast_gr0_set_size_dist <- get_size_distribution(yeast_gr0_sets)

mutans_r1_set_size <- get_size_list(mutans_r1_sets)
mutans_r1_set_size_dist <- get_size_distribution(mutans_r1_sets)
yeast_r1_set_size <- get_size_list(yeast_r1_sets)
yeast_r1_set_size_dist <- get_size_distribution(yeast_r1_sets)

mutans_gr1_set_size <- get_size_list(mutans_gr1_sets)
mutans_gr1_set_size_dist <- get_size_distribution(mutans_gr1_sets)
yeast_gr1_set_size <- get_size_list(yeast_gr1_sets)
yeast_gr1_set_size_dist <- get_size_distribution(yeast_gr1_sets)

mutans_g1_set_size <- get_size_list(mutans_g1_sets)
mutans_g1_set_size_dist <- get_size_distribution(mutans_g1_sets)

size_data <- matrix(nrow = 2, ncol = 12)
rownames(size_data) <- c('mutans', 'yeast')
colnames(size_data) <- c('r0_size', 'r0_size_dist', 'g0_size', 'g0_size_dist', 'gr0_size', 'gr0_size_dist', 'r1_size', 'r1_size_dist', 
                         'gr1_size', 'gr1_size_dist', 'g1_size', 'g1_size_dist')
size_data[1,1] <- list(mutans_r0_set_size)
size_data[1,2] <- list(mutans_r0_set_size_dist)
size_data[1,3] <- list(mutans_g0_set_size)
size_data[1,4] <- list(mutans_g0_set_size_dist)
size_data[1,5] <- list(mutans_gr0_set_size)
size_data[1,6] <- list(mutans_gr0_set_size_dist)
size_data[1,7] <- list(mutans_r1_set_size)
size_data[1,8] <- list(mutans_r1_set_size_dist)
size_data[1,9] <- list(mutans_gr1_set_size)
size_data[1,10] <- list(mutans_gr1_set_size_dist)
size_data[1,11] <- list(mutans_g1_set_size)
size_data[1,12] <- list(mutans_g1_set_size_dist)

size_data[2,1] <- list(yeast_r0_set_size)
size_data[2,2] <- list(yeast_r0_set_size_dist)
size_data[2,3] <- list(yeast_g0_set_size)
size_data[2,4] <- list(yeast_g0_set_size_dist)
size_data[2,5] <- list(yeast_gr0_set_size)
size_data[2,6] <- list(yeast_gr0_set_size_dist)
size_data[2,7] <- list(yeast_r1_set_size)
size_data[2,8] <- list(yeast_r1_set_size_dist)
size_data[2,9] <- list(yeast_gr1_set_size)
size_data[2,10] <- list(yeast_gr1_set_size_dist)
# size_data[1,11] <- list(yeast_g1_set_size)
# size_data[1,12] <- list(yeast_g1_set_size_dist)