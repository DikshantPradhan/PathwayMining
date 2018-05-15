rm(list = ls())
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/gene_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
library(parallel)
library(tictoc)

model_data_generation <- function(sybil_model, grb_model, model_name,
                                  r0 = TRUE, gr0 = TRUE, g0 = TRUE, r1 = TRUE, gr1 = TRUE, g1 = TRUE, gpr_save = FALSE, cores = 1){

  # setup
  falcon_model <- GRB_generate_falcon_model(sybil_model)
  n <- grb_model$get_sizes()$NumVars
  vars <- grb_model$get_names()$VarName

  falcon_n <- falcon_model$get_sizes()$NumVars
  falcon_vars <- falcon_model$get_names()$VarName
  non_vars <- grep('conversion', falcon_vars)

  non_gene_assc_rxns <- which(sybil_model@genes == "")
  gene_indexes <- grep('Ex_a', falcon_model$get_names()$VarName)
  falcon_rxn_idxs <- c(non_gene_assc_rxns, gene_indexes)

  falcon_full_rxn_idxs <- 1:(non_vars[1]-1)

  falcon_gene_rxn_idxs <- grep('Ex_a_', falcon_vars)
  if (length(falcon_rxn_idxs) == 0){
    falcon_rxn_idxs <- 1:(falcon_gene_rxn_idxs[length(falcon_gene_rxn_idxs)])
  }

  gpr <- generate_gpr(sybil_model)
  if (gpr_save){
    save(gpr, file = paste('final_paper_data/', model_name, '_gpr.RData', sep = ''))
  }

  # R0 sets
  if (r0){
    print('r0 sets...')
    print(model_name)
    ptm <- proc.time()
    tic()
    r0_coupling_mtx <- flux_coupling_raptor(grb_model)$coupled
    toc()
    print(proc.time() - ptm)
    ptm <- proc.time()
    tic()
    r0_sets <- get_list_of_sets_from_mtx(r0_coupling_mtx) #get_list_of_sets(return_couples(r0_coupling_mtx))
    toc()
    print(proc.time() - ptm)

    save(r0_sets, file = paste('final_paper_data/', model_name, '_r0_sets.RData', sep = ''))

    # GR0 sets
    if (gr0){
      print('gr0 sets...')
      gr0_sets <- gene_set_from_rxn_set(gpr, r0_sets)
      save(gr0_sets, file = paste('final_paper_data/', model_name, '_gr0_sets.RData', sep = ''))
    }
  }

  # G0 sets
  if (g0){
    print('g0 sets...')
    print(model_name)

    ptm <- proc.time()
    tic()
    g0_coupling_mtx <- flux_coupling_raptor(falcon_model, reaction_indexes = falcon_full_rxn_idxs)$coupled
    toc()
    print(proc.time() - ptm)
    save(g0_coupling_mtx, file = paste('final_paper_data/', model_name, '_g0_coupling_mtx.RData', sep = ''))
    ptm <- proc.time()
    tic()
    #g0_sets <- get_list_of_sets(return_couples(g0_coupling_mtx)) #GRB_generate_set_list(mutans_falcon)
    g0_sets <- get_list_of_sets_from_mtx(g0_coupling_mtx)
    toc()
    print(proc.time() - ptm)

    save(g0_sets, file = paste('final_paper_data/', model_name, '_g0_sets.RData', sep = ''))
  }

  # R1 Sets
  if (r1){
    print('r1 sets...')
    print(model_name)

    ptm <- proc.time()
    tic()
    coupling_vector_list <- GRB_generate_set_lists_cluster(grb_model, 1:n, compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = cores)
    toc()
    print(proc.time() - ptm)
    save(coupling_vector_list, file = paste('final_paper_data/', model_name, '_r1_coupling_vector_list.RData', sep = ''))
    ptm <- proc.time()
    tic()
    r1_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector_list, n, vars) #coupling_matrix_from_array(mutans_coupling_array)
    r1_sets <- get_list_of_sets_from_mtx(r1_matrix)
    toc()
    print(proc.time() - ptm)

    save(r1_sets, file = paste('final_paper_data/', model_name, '_r1_sets.RData', sep = ''))

    # GR1 Sets
    if (gr1){
      print('gr1 sets...')
      gr1_sets <- gene_set_from_rxn_set(gpr, r1_sets)

      save(gr1_sets, file = paste('final_paper_data/', model_name, '_gr1_sets.RData', sep = ''))
    }
  }

  # G1 Sets
  if (g1){
    non_gene_assc_rxns <- which(sybil_model@genes == "")
    gene_indexes <- grep('Ex_a', falcon_vars)
    suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)

    print('g1 sets...')
    print(model_name)
    ptm <- proc.time()
    tic()
    coupling_vector_list <- GRB_generate_set_lists_cluster(falcon_model, suppression_idxs = suppr_indexes, reaction_indexes = suppr_indexes,
                                                           compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = cores)
    toc()
    print(proc.time() - ptm)
    save(coupling_vector_list, file = paste('final_paper_data/', model_name, '_g1_coupling_vector_list.RData', sep = ''))
    ptm <- proc.time()
    tic()
    g1_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector_list, falcon_n, falcon_vars)
    g1_sets <- get_list_of_sets_from_mtx(g1_matrix)
    toc()
    print(proc.time() - ptm)

    save(g1_sets, file = paste('final_paper_data/', model_name, '_g1_sets.RData', sep = ''))
  }

  # sets size distributions
}
