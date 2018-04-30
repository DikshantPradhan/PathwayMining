load("~/GitHub/PathwayMining/gi_matrix.RData")
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/gene_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g1_sets.RData")
load("~/GitHub/PathwayMining/scripts/final_paper_data/yeast_g0_sets.RData")

e_matrix <- gi_matrix$e

clean_row_names <- function(names){
  for (i in 1:length(names)){
    names[i] <- strsplit(names[i], split = "_")[[1]][1]
  }
  return(names)
}

rownames(e_matrix) <- clean_row_names(rownames(e_matrix))
colnames(e_matrix) <- clean_row_names(colnames(e_matrix))

g0 <- clean_rxn_names_in_set(g0_sets)
g1 <- clean_rxn_names_in_set(yeast_g1_sets)

g1_interact <- get_num_interactions(g1, enrichment_mtx, 0)
print(g1_interact)
print(get_num_interactions_in_set(g1[[986]], e_matrix, 0))
