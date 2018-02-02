# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/falcon_tools.R')

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
