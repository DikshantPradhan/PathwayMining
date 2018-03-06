source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/data_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')

# test_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled
#
# save(test_mtx, file = 'test_mtx.RData')
#
# mutans_test_matrix <- isolate_gene_matrix(test_mtx)
# clean_mutans_test_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_test_matrix)))[[1]])
#
# mutans_falcon <- GRB_mutans_falcon_model()
# mutans_falcon_coupling_array_test <- GRB_generate_set_lists_array(mutans_falcon, suppression_idxs = c(1), reaction_indexes = reaction_indexes, compare_known_r0_sets = TRUE, optimize_suppr = TRUE)
# mutans_test_matrix <- isolate_gene_matrix(mutans_falcon_coupling_array_test[,,1])
# clean_mutans_test_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_test_matrix)))[[1]])
# length(unlist(clean_mutans_test_set))

load('~/GitHub/PathwayMining/data/mutans_model/mutans_model_w_obj.RData')
load('~/GitHub/PathwayMining/data/mutans_model/mutans_gene_dels.RData')
# model name: 'mutans_obj'
sybil_mutans_obj <- mutans_obj
mutans_obj <- as_GRBmodel(mutans_obj)
mutans_obj$show_output(FALSE)

single_gene_ko_flux_dels <- singleGeneDels@fluxdels
double_gene_ko_flux_dels <- doubleGeneDels@fluxdels

single_gene_ko_max_flux <- matrix(nrow = length(single_gene_ko_flux_dels), ncol = 1)
double_gene_ko_max_flux <- matrix(nrow = length(double_gene_ko_flux_dels), ncol = 1)

print(paste('single gene kos', length(single_gene_ko_max_flux)))
for (i in 1:length(single_gene_ko_flux_dels)){
  single_gene_ko_max_flux[i] <- GRB_maximize(mutans_obj, 477, single_gene_ko_flux_dels[[i]])
}
print(paste('double gene kos', length(double_gene_ko_max_flux)))
for (i in 1:length(double_gene_ko_flux_dels)){
  double_gene_ko_max_flux[i] <- GRB_maximize(mutans_obj, 477, double_gene_ko_flux_dels[[i]])
}

print('lethal single dels:')
lethal_single_dels <- which(near(0, single_gene_ko_max_flux))
print(lethal_single_dels)
print('lethal double dels:')
lethal_double_dels <- which(near(0, double_gene_ko_max_flux))
print(lethal_double_dels)

#save(double_gene_ko_max_flux, single_gene_ko_max_flux, file = '~/GitHub/PathwayMining/data/mutans_model/deletion/mutans_gene_del_flux.RData')

r1_set_containment <- matrix(FALSE, nrow = length(double_gene_ko_flux_dels), ncol = 1)
for (i in 1:length(double_gene_ko_flux_dels)){
  r1_set_containment[i] <- check_sets_for_containing(double_gene_ko_flux_dels[[i]], mutans_r1_sets)
}

lethal_flux_dels <- double_gene_ko_flux_dels[lethal_double_dels]
r1_lethal_set_containment <- matrix(FALSE, nrow = length(lethal_flux_dels), ncol = 1)
for (i in 1:length(lethal_flux_dels)){
  r1_set_containment[i] <- check_sets_for_containing(lethal_flux_dels[[i]], mutans_r1_sets)
}
