# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
#source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/falcon_tools.R')

## ECOLI MODEL

ptm <- proc.time() # timing start

ecoli <- GRB_ecoli_model()
n <- ecoli$get_sizes()$NumVars
vars <- ecoli$get_names()$VarName

ecoli_r0_set_list <- GRB_generate_set_list(ecoli)

#proc.time() - ptm # timing end

ecoli <- GRB_ecoli_model()
ecoli_set_lists <- GRB_generate_set_lists(ecoli, 1:n, compare_known_r0_sets = TRUE, optimize_suppr=TRUE) #GRB_generate_set_lists(ecoli, ecoli_r0_set_list, 1:n)

ecoli_coupling_array <- GRB_generate_set_lists_array(ecoli, 1:n, compare_known_r0_sets = TRUE, optimize_suppr=TRUE)
ecoli_r1_matrix <- coupling_matrix_from_array(ecoli_coupling_array)
ecoli_r1_matrix <- (ecoli_r1_matrix > 0)
ecoli_r1_sets <- list(get_list_of_sets(return_couples(ecoli_r1_matrix)))

ecoli_composition_set_full <- return_composition_sets(ecoli_r0_set_list, ecoli_set_lists, ecoli)
ecoli_composition_set <- ecoli_composition_set_full$composition

check_composition_error(ecoli_composition_set_full, ecoli_r0_set_list)

print(ecoli_composition_set_full$error)
print('deletion error check')

ecoli_deletion_list <- check_set_list_for_deletion(vars, ecoli_set_lists)
check_deletion_error(ecoli_deletion_list, ecoli_r0_set_list)

save(ecoli_r0_set_list, ecoli_set_lists, ecoli_composition_set_full, ecoli_deletion_list, file = "ecoli_run_data.RData")

#proc.time() - ptm # timing end

ecoli_r0_pairs <- return_pairs_from_set_list(ecoli_r0_set_list)
ecoli_new_r1_pairs <- new_pairs_from_composition(ecoli_r0_set_list, ecoli_composition_set)
ecoli_r1_pairs <- append_pair_lists(ecoli_r0_pairs, ecoli_new_r1_pairs)
ecoli_r0_set_list <- ecoli_r0_set_list
ecoli_r1_set_list <- get_list_of_sets(ecoli_r1_pairs)
save(ecoli_r0_pairs, ecoli_new_r1_pairs, ecoli_r1_pairs, ecoli_r0_set_list, ecoli_r1_set_list, file = 'ecoli_pairs_sets.RData')

print('FIN')


## MUTANS MODEL

ptm <- proc.time() # timing start

mutans <- GRB_mutans_model()
n <- mutans$get_sizes()$NumVars
vars <- mutans$get_names()$VarName

mutans_r0_set_list <- GRB_generate_set_list(mutans)

#proc.time() - ptm # timing end

mutans_r0_coupling_mtx <- flux_coupling_raptor(mutans)$coupled
mutans_r0_sets <- list(get_list_of_sets(return_couples(mutans_r0_coupling_mtx)))


mutans <- GRB_mutans_model()
ptm <- proc.time() # timing start
mutans_coupling_array <- GRB_generate_set_lists_array(mutans, compare_known_r0_sets = TRUE, optimize_suppr = TRUE)
mutans_r1_matrix <- coupling_matrix_from_array(mutans_coupling_array)
mutans_r1_matrix <- (mutans_r1_matrix > 0)
mutans_r1_sets <- get_list_of_sets(return_couples(mutans_r1_matrix))
proc.time() - ptm
save(mutans_r1_matrix, file = 'data/mutans_model/mutans_r1_matrix.RData')

mutans_set_lists <- GRB_generate_set_lists(mutans, mutans_r0_set_list, 1:n)

mutans_composition_set_full <- return_composition_sets(mutans_r0_set_list, mutans_set_lists, mutans)
mutans_composition_set <- mutans_composition_set_full$composition

check_composition_error(mutans_composition_set_full, mutans_r0_set_list)

print(mutans_composition_set_full$error)
print('deletion error check')

mutans_deletion_list <- check_set_list_for_deletion(vars, mutans_set_lists)
check_deletion_error(mutans_deletion_list, mutans_r0_set_list)

save(mutans_r0_set_list, mutans_set_lists, mutans_composition_set_full, mutans_deletion_list, file = "mutans_run_data.RData")

proc.time() - ptm # timing end

mutans_r0_pairs <- return_pairs_from_set_list(mutans_r0_set_list)
mutans_new_r1_pairs <- new_pairs_from_composition(mutans_r0_set_list, mutans_composition_set)
mutans_r1_pairs <- append_pair_lists(mutans_r0_pairs, mutans_new_r1_pairs)
mutans_r0_set_list <- mutans_r0_set_list
mutans_r1_set_list <- get_list_of_sets(mutans_r1_pairs)
save(mutans_r0_pairs, mutans_new_r1_pairs, mutans_r1_pairs, mutans_r0_set_list, mutans_r1_set_list, file = 'mutans_pairs_sets.RData')

print('FIN')


## YEAST MODEL

ptm <- proc.time() # timing start

yeast <- GRB_yeast_model()
n <- yeast$get_sizes()$NumVars
vars <- yeast$get_names()$VarName

yeast_r0_set_list <- GRB_generate_set_list(yeast)

proc.time() - ptm # timing end


yeast_r0_coupling_mtx <- flux_coupling_raptor(yeast)$coupled
yeast_r0_sets <- list(get_list_of_sets(return_couples(yeast_r0_coupling_mtx)))


yeast <- GRB_yeast_model()
ptm <- proc.time() # timing start
yeast_coupling_array <- GRB_generate_set_lists_array(yeast, compare_known_r0_sets = TRUE, optimize_suppr = TRUE)
yeast_r1_matrix <- coupling_matrix_from_array(yeast_coupling_array)
yeast_r1_matrix <- (yeast_r1_matrix > 0)
yeast_r1_sets <- list(get_list_of_sets(return_couples(yeast_r1_matrix)))
proc.time() - ptm
save(yeast_r1_matrix, file = 'data/yeast_model/yeast_r1_matrix.RData')


yeast <- GRB_yeast_model()
yeast_set_lists <- GRB_generate_set_lists(yeast, yeast_r0_set_list, 1:n)

yeast_composition_set_full <- return_composition_sets(yeast_r0_set_list, yeast_set_lists, yeast)
yeast_composition_set <- yeast_composition_set_full$composition

check_composition_error(yeast_composition_set_full, yeast_r0_set_list)

print(yeast_composition_set_full$error)
print('deletion error check')

yeast_deletion_list <- check_set_list_for_deletion(vars, yeast_set_lists)
check_deletion_error(yeast_deletion_list, yeast_r0_set_list)

save(yeast_r0_set_list, yeast_set_lists, yeast_composition_set_full, yeast_deletion_list, file = "yeast_run_data.RData")

#proc.time() - ptm # timing end

yeast_r0_pairs <- return_pairs_from_set_list(yeast_r0_set_list)
yeast_new_r1_pairs <- new_pairs_from_composition(yeast_r0_set_list, yeast_composition_set)
yeast_r1_pairs <- append_pair_lists(yeast_r0_pairs, yeast_new_r1_pairs)
yeast_r0_set_list <- yeast_r0_set_list
yeast_r1_set_list <- get_list_of_sets(yeast_r1_pairs)
save(yeast_r0_pairs, yeast_new_r1_pairs, yeast_r1_pairs, yeast_r0_set_list, yeast_r1_set_list, file = 'yeast_pairs_sets.RData')

rm(list = ls())
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/gene_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/data_tools.R')
library(parallel)
library(tictoc)
library(sybil)

load('~/GitHub/PathwayMining/data/yeast_model/Maranas_model/maranas_model_lipid_exch.RData')

yeast <- GRB_yeast_model()
n <- yeast$get_sizes()$NumVars
vars <- yeast$get_names()$VarName

coupling_vector <- read_coupling_csv('~/GitHub/PathwayMining/scripts/yeast_r1_coupling.csv')
avoid_rxns <- coupling_vector$completed_idxs

ptm <- proc.time()
tic()
coupling_vector_list <- GRB_generate_set_lists_cluster(yeast, compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, file_output = '~/GitHub/PathwayMining/scripts/yeast_r1_coupling.csv')
save(coupling_vector_list, file = '~/GitHub/PathwayMining/scripts/yeast_g1_coupling_vector.RData')
toc()
print(proc.time() - ptm)

ptm <- proc.time()
tic()
yeast_g1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
yeast_g1_sets <- get_list_of_sets_from_mtx(yeast_g1_coupling_matrix)
toc()
print(proc.time() - ptm)
save(yeast_g1_sets, file = '~/GitHub/PathwayMining/scripts/final_paper_data/yeast_r1_sets.RData')


print('FIN')


# PAO MODEL
rm(list = ls())
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/gene_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/data_tools.R')
library(parallel)
library(tictoc)
library(sybil)
pao <- GRB_pao_model()
n <- pao$get_sizes()$NumVars
vars <- pao$get_names()$VarName

output_file <- '~/GitHub/PathwayMining/scripts/pao_r1_coupling.csv'

coupling_vector <- read_coupling_csv(output_file)
avoid_rxns <- coupling_vector$completed_idxs

ptm <- proc.time()
tic()
coupling_vector_list <- GRB_generate_set_lists_cluster(pao, compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, file_output = output_file)
save(coupling_vector_list, file = 'pao_coupling_vector.RData')
toc()
print(proc.time() - ptm)

ptm <- proc.time()
tic()
pao_r1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
pao_r1_sets <- get_list_of_sets_from_mtx(pao_r1_coupling_matrix)
toc()
print(proc.time() - ptm)
save(pao_r1_sets, file = '~/GitHub/PathwayMining/scripts/final_paper_data/pao_r1_sets.RData')


pao_r0_coupling_mtx <- flux_coupling_raptor(pao)$coupled
r0_sets <- get_list_of_sets_from_mtx(pao_r0_coupling_mtx)
pao_falcon <- GRB_pao_falcon_model()
n <- pao_falcon$get_sizes()$NumVars
vars <- pao_falcon$get_names()$VarName

pao_g0_coupling_mtx <- flux_coupling_raptor(pao_falcon, reaction_indexes = 1:2699)$coupled
g0_sets <- get_list_of_sets_from_mtx(pao_g0_coupling_mtx)

check_all_sets_for_containing(r0_sets, g0_sets)

load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')
pao_falcon_model <- GRB_generate_falcon_model(pao_model)
n <- pao_falcon_model$get_sizes()$NumVars
vars <- pao_falcon_model$get_names()$VarName

non_gene_assc_rxns <- which(pao_model@genes == "")
gene_indexes <- grep('Ex_a', pao_falcon_model$get_names()$VarName)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)

pao_g0_coupling_mtx_2 <- flux_coupling_raptor(pao_falcon_model, reaction_indexes = suppr_indexes)$coupled
g0_sets_2 <- get_list_of_sets_from_mtx(pao_g0_coupling_mtx)

check_all_sets_for_containing(r0_sets, g0_sets_2)
check_all_sets_for_containing(g0_sets, g0_sets_2)
