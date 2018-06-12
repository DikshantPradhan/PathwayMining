# PACT Analysis

# initialization
source('~/GitHub/PathwayMining/data_tools.R')
pao_g1_coupling <- read_coupling_csv('~/Documents/jensn lab/PAO_model/pao_g1_coupling.csv')
pao_g1_coupling_list <- pao_g1_coupling$coupling_vector

vars <- pao_falcon@react_id

full_pao_coupling_mtx <- coupling_matrix_from_coupling_vector_list(pao_g1_coupling_list, n_react = length(vars), vars = vars)
fullish_pao_coupling_mtx <- full_ish_coupling_matrix_from_coupling_vector_list(pao_g1_coupling_list, n_react = length(vars), vars = vars)

uncoupled_mtx <- identify_intermediate_uncoupled(full_pao_coupling_mtx, fullish_pao_coupling_mtx)