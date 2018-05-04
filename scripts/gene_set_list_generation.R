# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
# source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/raptor_coupling.R')
source('~/GitHub/PathwayMining/grb_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')
source('~/GitHub/PathwayMining/falcon_tools.R')

ptm <- proc.time()

# # test composition of set_lists (make sure that blocking any reaction in an og_set results in the same g1 set)

## ECOLI MODEL

ecoli_falcon <- GRB_ecoli_falcon_model()
n <- ecoli_falcon$get_sizes()$NumVars
vars <- ecoli_falcon$get_names()$VarName

reaction_indexes <- c()
reaction_indexes <- grep('Ex_a', vars)

print(paste('num genes:', length(reaction_indexes)))

ecoli_falcon_og_set_list <- GRB_generate_set_list(ecoli_falcon, reaction_indexes = reaction_indexes)

ecoli_falcon <- GRB_ecoli_falcon_model()
ecoli_falcon_set_lists <- GRB_generate_set_lists(ecoli_falcon, ecoli_falcon_og_set_list, 1:n, reaction_indexes)

ecoli_falcon_composition_set_full <- return_composition_sets(ecoli_falcon_og_set_list, ecoli_falcon_set_lists, ecoli_falcon)

ecoli_falcon_composition_set <- ecoli_falcon_composition_set_full$composition

print('composition error check')

# check for errors in compositions of g1 sets (intermediate) (each blockage within the same set should have the same effects)
for (i in 1:length(ecoli_falcon_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_falcon_og_set_list[[i]])
  if (length(ecoli_falcon_og_set_list[[i]]) > 1){
    for (j in 1:(length(ecoli_falcon_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, ecoli_falcon_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, ecoli_falcon_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      #print(paste(ecoli_falcon_composition_set[[rxn1]], ecoli_falcon_composition_set[[rxn2]]))
      if (rxn1 > length(ecoli_falcon_composition_set) | rxn2 > length(ecoli_falcon_composition_set)){next}

      if (length(ecoli_falcon_composition_set[[rxn1]]) != length(ecoli_falcon_composition_set[[rxn2]])){
        print("composition error")
        print(c(ecoli_falcon_og_set_list[[i]][[j]], ";", ecoli_falcon_og_set_list[[i]][[j+1]]))
        print(c(ecoli_falcon_composition_set[[rxn1]], ";", ecoli_falcon_composition_set[[rxn2]]))
      }
      #print(j)
      #print(ecoli_falcon_composition_set[[GRB_get_rxn_idx(ecoli_falcon, j)]])
    }
  }
}

print(ecoli_falcon_composition_set_full$error)

print('deletion error check')

ecoli_falcon_deletion_list <- check_set_list_for_deletion(vars, ecoli_falcon_set_lists)

# check for inconsistencies in deletions of sets (each blockage within the same set should have the same effects)
for (i in 1:length(ecoli_falcon_og_set_list)){ # print sets joined by each deletion
  #print(ecoli_falcon_og_set_list[[i]])
  if (length(ecoli_falcon_og_set_list[[i]]) > 1){
    for (j in 1:(length(ecoli_falcon_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, ecoli_falcon_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, ecoli_falcon_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      if (length(ecoli_falcon_deletion_list[[rxn1]]) != length(ecoli_falcon_deletion_list[[rxn2]])){
        print("deletion error")
        print(c(ecoli_falcon_og_set_list[[i]][[j]], ";", ecoli_falcon_og_set_list[[i]][[j+1]]))
        print(c(ecoli_falcon_deletion_list[[rxn1]], ";", ecoli_falcon_deletion_list[[rxn2]]))
      }
      #print(j)
      #print(ecoli_falcon_composition_set[[GRB_get_rxn_idx(ecoli_falcon, j)]])
    }
  }
}

save(ecoli_falcon_og_set_list, ecoli_falcon_set_lists, ecoli_falcon_composition_set_full, ecoli_falcon_deletion_list, file = "ecoli_falcon_run_data.RData")

proc.time() - ptm

print('specific error')

#ct <- 0
#blocked <- c()
#for (i in 1:length(ecoli_composition_set_full$error)){
#  if (nchar(ecoli_composition_set_full$error[i]) > 8){
#    ct <- ct + 1
#    #blocked <- c(blocked, ecoli_vars[i])
#    print(ecoli_composition_set_full$error[i])
#    #print(ecoli_composition_set_full$error[i])
#  }
#}

print(ct)

ecoli_falcon_g0_pairs <- return_pairs_from_set_list(ecoli_falcon_og_set_list)
ecoli_falcon_new_g1_pairs <- new_pairs_from_composition(ecoli_falcon_og_set_list, ecoli_falcon_composition_set)
ecoli_falcon_g1_pairs <- append_pair_lists(ecoli_falcon_g0_pairs, ecoli_falcon_new_g1_pairs)

ecoli_falcon_g0_set_list <- ecoli_falcon_og_set_list
ecoli_falcon_g1_set_list <- get_list_of_sets(ecoli_falcon_g1_pairs)


## MUTANS MODEL ~ 6 HRS

load('~/GitHub/PathwayMining/data/mutans_model/mutans_model.RData')

sybil_mutans <- mutans
non_gene_assc_rxns <- which(sybil_mutans@genes == "")

mutans_falcon <- GRB_mutans_falcon_model()
n <- mutans_falcon$get_sizes()$NumVars
vars <- mutans_falcon$get_names()$VarName

gene_indexes <- c()
gene_indexes <- grep('Ex_a', vars)
suppr_indexes <- c(gene_indexes, non_gene_assc_rxns)
reaction_indexes <- 1:1115

mutans_g0_coupling_mtx <- flux_coupling_raptor(mutans_falcon, reaction_indexes = reaction_indexes)$coupled
mutans_falcon_g0_sets <- get_list_of_sets(return_couples(mutans_g0_coupling_mtx))
mutans_g0_matrix <- isolate_gene_matrix(fill_coupling_matrix(mutans_g0_coupling_mtx))
clean_mutans_g0_sets <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_g0_matrix)))[[1]])
save(mutans_g0_coupling_mtx, file = 'data/mutans_model/mutans_g0_coupling_mtx.RData')

print(paste('num genes:', length(reaction_indexes)))
#ptm <- proc.time() # timing start
#mutans_falcon_og_set_list <- GRB_generate_set_list(mutans_falcon, reaction_indexes = reaction_indexes)

#proc.time() - ptm # timing end
ptm <- proc.time() # timing start
mutans_falcon <- GRB_mutans_falcon_model()
mutans_falcon_coupling_array <- GRB_generate_set_lists_array(mutans_falcon, reaction_indexes = reaction_indexes, compare_known_r0_sets = TRUE, optimize_suppr = TRUE)
mutans_falcon_g1_matrix <- coupling_matrix_from_array(mutans_falcon_coupling_array)
mutans_falcon_g1_matrix <- (mutans_falcon_g1_matrix > 0)
proc.time() - ptm
save(mutans_falcon_coupling_array, file = '~/GitHub/PathwayMining/data/mutans_model/mutans_falcon_coupling_array.RData')
save(mutans_falcon_g1_matrix, file = '~/GitHub/PathwayMining/data/mutans_model/mutans_falcon_g1_matrix.RData')

mutans_falcon_g1_sets <- get_list_of_sets(return_couples(mutans_falcon_g1_matrix))
mutans_g1_matrix <- isolate_gene_matrix(fill_coupling_matrix(mutans_falcon_g1_matrix))
clean_mutans_g1_set <- clean_rxn_names_in_set(list(get_list_of_sets(return_couples(mutans_g1_matrix)))[[1]])

mutans_falcon_set_lists <- GRB_generate_set_lists(mutans_falcon, mutans_falcon_og_set_list, suppr_indexes, reaction_indexes)
mutans_falcon_composition_set_full <- return_composition_sets(mutans_falcon_og_set_list, mutans_falcon_set_lists, mutans_falcon)

mutans_falcon_composition_set <- mutans_falcon_composition_set_full$composition

print('composition error check')

# check for errors in compositions of g1 sets (intermediate) (each blockage within the same set should have the same effects)
for (i in 1:length(mutans_falcon_og_set_list)){ # print sets joined by each deletion
  #print(mutans_falcon_og_set_list[[i]])
  if (length(mutans_falcon_og_set_list[[i]]) > 1){
    for (j in 1:(length(mutans_falcon_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, mutans_falcon_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, mutans_falcon_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      #print(paste(mutans_falcon_composition_set[[rxn1]], mutans_falcon_composition_set[[rxn2]]))
      if (rxn1 > length(mutans_falcon_composition_set) | rxn2 > length(mutans_falcon_composition_set)){next}

      if (length(mutans_falcon_composition_set[[rxn1]]) != length(mutans_falcon_composition_set[[rxn2]])){
        print("composition error")
        print(c(mutans_falcon_og_set_list[[i]][[j]], ";", mutans_falcon_og_set_list[[i]][[j+1]]))
        print(c(mutans_falcon_composition_set[[rxn1]], ";", mutans_falcon_composition_set[[rxn2]]))
      }
      #print(j)
      #print(mutans_falcon_composition_set[[GRB_get_rxn_idx(mutans_falcon, j)]])
    }
  }
}

print(mutans_falcon_composition_set_full$error)

print('deletion error check')

mutans_falcon_deletion_list <- check_set_list_for_deletion(vars, mutans_falcon_set_lists)

# check for inconsistencies in deletions of sets (each blockage within the same set should have the same effects)
for (i in 1:length(mutans_falcon_og_set_list)){ # print sets joined by each deletion
  #print(mutans_falcon_og_set_list[[i]])
  if (length(mutans_falcon_og_set_list[[i]]) > 1){
    for (j in 1:(length(mutans_falcon_og_set_list[[i]])-1)){
      rxn1 <- get_rxn_idx(vars, mutans_falcon_og_set_list[[i]][[j]])
      rxn2 <- get_rxn_idx(vars, mutans_falcon_og_set_list[[i]][[j+1]])
      #print(c(rxn1, rxn2))
      if (length(mutans_falcon_deletion_list[[rxn1]]) != length(mutans_falcon_deletion_list[[rxn2]])){
        print("deletion error")
        print(c(mutans_falcon_og_set_list[[i]][[j]], ";", mutans_falcon_og_set_list[[i]][[j+1]]))
        print(c(mutans_falcon_deletion_list[[rxn1]], ";", mutans_falcon_deletion_list[[rxn2]]))
      }
      #print(j)
      #print(mutans_falcon_composition_set[[GRB_get_rxn_idx(mutans_falcon, j)]])
    }
  }
}

save(mutans_falcon_og_set_list, mutans_falcon_set_lists, mutans_falcon_composition_set_full, mutans_falcon_deletion_list, file = "mutans_falcon_run_data.RData")

proc.time() - ptm

print('specific error')

#ct <- 0
#blocked <- c()
#for (i in 1:length(mutans_composition_set_full$error)){
#  if (nchar(mutans_composition_set_full$error[i]) > 8){
#    ct <- ct + 1
#    #blocked <- c(blocked, mutans_vars[i])
#    print(mutans_composition_set_full$error[i])
#    #print(mutans_composition_set_full$error[i])
#  }
#}

#print(ct)

mutans_falcon_g0_pairs <- return_pairs_from_set_list(mutans_falcon_og_set_list)
mutans_falcon_new_g1_pairs <- new_pairs_from_composition(mutans_falcon_og_set_list, mutans_falcon_composition_set)
mutans_falcon_g1_pairs <- append_pair_lists(mutans_falcon_g0_pairs, mutans_falcon_new_g1_pairs)

mutans_falcon_g0_set_list <- mutans_falcon_og_set_list
mutans_falcon_g1_set_list <- get_list_of_sets(mutans_falcon_g1_pairs)

save(mutans_falcon_g0_pairs, mutans_falcon_new_g1_pairs, mutans_falcon_g1_pairs, mutans_falcon_g0_set_list, mutans_falcon_g1_set_list, file = 'mutans_falcon_pairs_sets.RData')

# id g1 set compositions

g0_sets <- clean_mutans_g0_sets
g1_sets <- clean_mutans_g1_set

g1_compositions <- matrix(data = 0, nrow = length(g1_sets), ncol = length(g0_sets))
for (i in 1:length(g1_sets)){
  composition <- find_composing_sets(g1_sets[[i]], g0_sets)
  for (j in composition){
    g1_compositions[i, j] <- 1
  }
}

g1_comp_num <- rowSums(g1_compositions)
print(which(g1_comp_num > 1))

## YEAST MODEL
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
yeast_falcon_model <- GRB_generate_falcon_model(yeast_model)

n <- yeast_falcon_model$get_sizes()$NumVars
vars <- yeast_falcon_model$get_names()$VarName

non_gene_assc_rxns <- which(yeast_model@genes == "")
gene_indexes <- grep('Ex_a', yeast_falcon_model$get_names()$VarName)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)

coupling_vector <- read_coupling_csv('~/GitHub/PathwayMining/scripts/yeast_g1_coupling.csv')
avoid_rxns <- coupling_vector$completed_idxs

ptm <- proc.time()
tic()
coupling_vector_list <- GRB_generate_set_lists_cluster(yeast_falcon_model, suppression_idxs = suppr_indexes, reaction_indexes = suppr_indexes, compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, file_output = '~/GitHub/PathwayMining/scripts/yeast_g1_coupling.csv')
save(coupling_vector_list, file = '~/GitHub/PathwayMining/scripts/yeast_g1_coupling_vector.RData')
toc()
print(proc.time() - ptm)

ptm <- proc.time()
tic()
yeast_g1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
yeast_g1_sets <- get_list_of_sets_from_mtx(yeast_g1_coupling_matrix)
toc()
print(proc.time() - ptm)
save(yeast_g1_sets, file = '~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g1_sets.RData')


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
load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')

output_file <- '~/GitHub/PathwayMining/scripts/pao_g1_coupling.csv'

pao_falcon_model <- GRB_generate_falcon_model(pao_model)

n <- pao_falcon_model$get_sizes()$NumVars
vars <- pao_falcon_model$get_names()$VarName

non_gene_assc_rxns <- which(pao_model@genes == "")
gene_indexes <- grep('Ex_a', pao_falcon_model$get_names()$VarName)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)

#coupling_vector <- read_coupling_csv(output_file)
avoid_rxns <- c() #coupling_vector$completed_idxs

ptm <- proc.time()
tic()
coupling_vector_list <- GRB_generate_set_lists_cluster(pao_falcon_model, suppression_idxs = suppr_indexes, reaction_indexes = suppr_indexes, compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, file_output = output_file)
save(coupling_vector_list, file = 'pao_coupling_vector.RData')
toc()
print(proc.time() - ptm)

ptm <- proc.time()
tic()
pao_g1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
pao_g1_sets <- get_list_of_sets_from_mtx(pao_g1_coupling_matrix)
toc()
print(proc.time() - ptm)
save(pao_g1_sets, file = '~/GitHub/PathwayMining/scripts/final_paper_data/pao_g1_sets.RData')


r0_coupling_mtx <- flux_coupling_raptor(pao_falcon_model, reaction_indexes = suppr_indexes)$coupled
r0_coupling_mtx <- fill_coupling_matrix(r0_coupling_mtx)

suppr_vector <- Matrix(data = FALSE, nrow = 1, ncol = n, sparse = TRUE)
suppr_vector[suppr_indexes] <- TRUE
i <- 1
while (i <= n){
  if (suppr_vector[i]){ # if tagged to be suppressed
    set_idx <- which(r0_coupling_mtx[,i])[1] # which is first reaction (row) i is coupled to
    if (!is.na(set_idx)){
      rxn_idxs <- which(r0_coupling_mtx[set_idx,]) # other reactions in set
      # only suppress first reaction in set since, theoretically, suppressing any should have the same effect
      suppr_vector[rxn_idxs] <- FALSE
      suppr_vector[rxn_idxs[1]] <- TRUE
      }
    else {
      suppr_vector[i] <- FALSE
    }
  }
  i <- i+1
}

init_suppr <- length(which(suppr_vector))


coupling_vector <- read_coupling_csv('~/GitHub/PathwayMining/scripts/pao_coupling.csv')
avoid_rxns <- coupling_vector$completed_idxs

suppr_vector[avoid_rxns] <- FALSE

iter <- 1
while (length(avoid_rxns) < init_suppr){
  print(paste(length(avoid_rxns), '/', 712, "...", sep = ''))

  suppr_vector[avoid_rxns] <- FALSE

  ptm <- proc.time()
  tic()
  coupling_vector_list <- GRB_generate_set_lists_cluster(pao_falcon_model, suppression_idxs = which(suppr_vector), reaction_indexes = suppr_indexes, compare_known_r0_sets = TRUE, optimize_suppr=FALSE, cores = 6, avoid_idxs = avoid_rxns, file_output = '~/GitHub/PathwayMining/scripts/pao_coupling.csv')
  save(coupling_vector_list, file = 'pao_coupling_vector.RData')
  toc()
  print(proc.time() - ptm)
  print(paste('FIN', iter))
  iter <- iter + 1
  coupling_vector <- read_coupling_csv('~/GitHub/PathwayMining/scripts/pao_coupling.csv')
  avoid_rxns <- coupling_vector$completed_idxs
}

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
load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')
pao_falcon_model <- GRB_generate_falcon_model(pao_model)

non_gene_assc_rxns <- which(pao_model@genes == "")
gene_indexes <- grep('Ex_a', pao_falcon_model$get_names()$VarName)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)
avoid_rxns <- coupling_vector$completed_idxs
ptm <- proc.time()
tic()
coupling_vector_list <- GRB_generate_set_lists_cluster(pao_falcon_model, suppression_idxs = suppr_indexes, reaction_indexes = suppr_indexes, compare_known_r0_sets = TRUE, optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, file_output = '~/GitHub/PathwayMining/scripts/pao_coupling.csv')
save(coupling_vector_list, file = 'pao_coupling_vector.RData')
toc()
print(proc.time() - ptm)
print('FIN2')
