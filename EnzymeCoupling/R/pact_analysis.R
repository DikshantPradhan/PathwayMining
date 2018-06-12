# PACT Analysis

# initialization
source('~/GitHub/PathwayMining/data_tools.R')
pao_g1_coupling <- read_coupling_csv('~/GitHub/PathwayMining/scripts/pao/pao_g1_coupling.csv')
pao_g1_coupling_list <- pao_g1_coupling$coupling_vector

# pao_falcon <- generate_falcon_model(pao_model)
vars <- pao_falcon$get_names()$VarName #pao_falcon@react_id

full_pao_coupling_mtx <- coupling_matrix_from_coupling_vector_list(pao_g1_coupling_list, n_react = length(vars), vars = vars)
fullish_pao_coupling_mtx <- full_ish_coupling_matrix_from_coupling_vector_list(pao_g1_coupling_list, n_react = length(vars), vars = vars)

uncoupled_mtx <- identify_intermediate_uncoupled(full_pao_coupling_mtx, fullish_pao_coupling_mtx, n_react = length(vars))

uncoupled_pairs <- c()
for (i in 1:nrow(uncoupled_mtx)){
  for (j in which(uncoupled_mtx[i,])){
    rxn_1 <- rownames(uncoupled_mtx)[i]
    rxn_2 <- colnames(uncoupled_mtx)[j]
    
    uncoupled_pairs <- c(uncoupled_pairs, list(c(rownames(uncoupled_mtx)[i], colnames(uncoupled_mtx)[j])))
    # if (grepl('PA', rxn_1) & grepl('PA', rxn_2)){
    #   uncoupled_pairs <- c(uncoupled_pairs, list(c(rownames(uncoupled_mtx)[i], colnames(uncoupled_mtx)[j])))
    # }
  }
  # for (j in 1:ncol(uncoupled_mtx)){
  #   if (uncoupled_mtx[i,j]){
  #     uncoupled_pairs <- c(uncoupled_pairs, list(c(rownames(uncoupled_mtx)[i], colnames(uncoupled_mtx)[j])))
  #   }
  # }
}

uncoupled_pairs <- clean_rxn_names_in_set(uncoupled_pairs)
uncoupled_genes <- unique(unlist(uncoupled_pairs))

pao_data <- pao1#pa14[which(pa14$condition == unique(pa14$condition)[1]),]

# if using pa14 dataset, swap pa14 genes for pao1 homologs
# load("~/GitHub/PathwayMining/EnzymeCoupling/data/pao1_key.Rdata")
# pao_data$gene <- pao1_key[pao_data$gene]

g0_set_df_og <- gene_set_dataframe(g0_sets)
g1_set_df_og <- gene_set_dataframe(g1_sets)

genes_in_model <- unique(unlist(g0_sets))

df_genes <- pao_data$gene[which(pao_data$sig_df)]
df_genes <- df_genes[which(df_genes %in% genes_in_model)]
de_genes <- pao_data$gene[which(pao_data$sig_de)]
de_genes <- de_genes[which(de_genes %in% genes_in_model)]
dfde_genes <- pao_data$gene[which(pao_data$sig_de & pao_data$sig_df)]
dfde_genes <- dfde_genes[which(dfde_genes %in% genes_in_model)]
num_df <- length(df_genes)
num_de <- length(de_genes)
num_dfde <- length(dfde_genes)

coupled_genes <- setdiff(genes_in_model, uncoupled_genes)

print(length(df_genes))
print(length(de_genes))
print(length(dfde_genes))

print(length(uncoupled_genes))
print(length(coupled_genes))

print(length(intersect(de_genes, uncoupled_genes)))
print(length(intersect(de_genes, coupled_genes)))

print(length(intersect(df_genes, uncoupled_genes)))
print(length(intersect(df_genes, coupled_genes)))

print(length(intersect(dfde_genes, uncoupled_genes)))
print(length(intersect(dfde_genes, coupled_genes)))

uncoupled_g0_set_idxs <- c()
for (gene in uncoupled_genes){
  uncoupled_g0_set_idxs <- c(uncoupled_g0_set_idxs, get_set_idx(gene, g0_df$clean.sets))
}
uncoupled_g0_set_idxs <- unique(uncoupled_g0_set_idxs)

coupled_g0_set_idxs <- setdiff(1:nrow(g0_df), uncoupled_g0_set_idxs)
# for (gene in coupled_genes){
#   coupled_g0_set_idxs <- c(coupled_g0_set_idxs, get_set_idx(gene, g0_df$clean.sets))
# }
# coupled_g0_set_idxs <- unique(coupled_g0_set_idxs)

length(intersect(which(g0_df$df_frac == 1), coupled_g0_set_idxs))
length(intersect(which(g0_df$df_frac == 1), uncoupled_g0_set_idxs))
length(intersect(which(g0_df$de_frac == 1), uncoupled_g0_set_idxs))
length(intersect(which(g0_df$de_frac == 1), coupled_g0_set_idxs))
length(intersect(which(g0_df$dfde_frac == 1), coupled_g0_set_idxs))
length(intersect(which(g0_df$dfde_frac == 1), uncoupled_g0_set_idxs))