test_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled

save(test_mtx, file = 'test_mtx.RData')

mutans_test_matrix <- isolate_gene_matrix(test_mtx)
clean_mutans_test_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_test_matrix)))[[1]])

mutans_falcon <- GRB_mutans_falcon_model()
mutans_falcon_coupling_array_test <- GRB_generate_set_lists_array(mutans_falcon, suppression_idxs = c(1), reaction_indexes = reaction_indexes, compare_known_r0_sets = TRUE, optimize_suppr = TRUE)
mutans_test_matrix <- isolate_gene_matrix(mutans_falcon_coupling_array_test[,,1])
clean_mutans_test_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_test_matrix)))[[1]])
length(unlist(clean_mutans_test_set))
