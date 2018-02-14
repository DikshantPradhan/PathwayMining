# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
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

proc.time() - ptm # timing end

ecoli <- GRB_ecoli_model()
ecoli_set_lists <- GRB_generate_set_lists(ecoli, ecoli_r0_set_list, 1:n)

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

mutans <- GRB_mutans_model()
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

print('FIN')
