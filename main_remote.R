source('grb_tools.R')

test <- coupling_matrix_from_array(mutans_falcon_coupling_array)
test <- (test > 0)

#mutans_falcon_g1_sets <- list(get_list_of_sets(return_couples(mutans_falcon_g1_matrix)))
mutans_g1_matrix <- isolate_gene_matrix(test)
clean_mutans_g1_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_g1_matrix)))[[1]])


#mutans_test_g1_matrix <- isolate_gene_matrix(mutans_falcon_coupling_array[,,100])
#clean_mutans_g1_set_test <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_test_g1_matrix)))[[1]])
