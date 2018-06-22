library(ggplot2)

gene_data_frame <- function(model, gr0, gr1){
  genes <- model@allGenes
  
  # df <- data.frame(genes)
  
  gpr_idxs <- lapply(genes, function(x){get_set_idx(x, model@genes)})
  gpr_promiscuity <- vapply(gpr_idxs, function(x){length(x)}, c(1))
  gr0_idxs <- lapply(genes, function(x){get_set_idx(x, gr0)})
  gr0_promiscuity <- vapply(gr0_idxs, function(x){length(x)}, c(1))
  gr1_idxs <- lapply(genes, function(x){get_set_idx(x, gr1)})
  gr1_promiscuity <- vapply(gr1_idxs, function(x){length(x)}, c(1))
  
  # df$gpr_promiscuity <- gpr_promiscuity
  # df$gr0_promiscuity <- gr0_promiscuity
  # df$gr1_promiscuity <- gr1_promiscuity
  # 
  # return(df)
  # list(genes = genes, gpr_idxs = gpr_idxs, gpr_promiscuity = gpr_promiscuity,
  #      gr0_idxs = gr0_idxs, gr0_promiscuity = gr0_promiscuity, gr1_idxs = gr1_idxs, gr1_promiscuity = gr1_promiscuity)
  # cbind(genes = genes, gpr_promiscuity = gpr_promiscuity, gr0_promiscuity = gr0_promiscuity, gr1_promiscuity = gr1_promiscuity)
  data.frame(genes = genes, gpr_promiscuity = gpr_promiscuity, gr0_promiscuity = gr0_promiscuity, gr1_promiscuity = gr1_promiscuity)
}


load("~/GitHub/PathwayMining/data/final_figures_data/pao_gr0_gr1_sets.RData")
pao_gene_data_frame <- gene_data_frame(pao_model, pao_gr0_sets, pao_gr1_sets)
pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gpr_promiscuity)),]
# pao_gene_data_frame <- pao_gene_data_frame[1:n_genes,]
pao_gene_data_frame$numbers <- 1:nrow(pao_gene_data_frame)

# pao_gene_data_frame <- data.frame(pao_gene_data_frame)
# pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gpr_promiscuity)),]

n_genes <- length(unique(unlist(pao_gr0_sets)))

ggplot(data=pao_gene_data_frame, aes(x=numbers, y=gpr_promiscuity)) +
  geom_bar(stat="identity") + ggtitle('Gene Promiscuity between Reactions') + labs(y='Number of Reactions', x ='Gene') #+
  # geom_bar(aes(fill = cut))+
  # theme(#axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())

pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gr0_promiscuity)),]
pao_gene_data_frame <- pao_gene_data_frame[1:n_genes,]
pao_gene_data_frame$numbers <- 1:n_genes
ggplot(data=pao_gene_data_frame, aes(x=numbers, y=gr0_promiscuity)) +
  geom_bar(stat="identity") + ggtitle('Gene Promiscuity between GR0 Sets') + labs(y='Number of Sets', x ='Gene')

pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gr1_promiscuity)),]
pao_gene_data_frame <- pao_gene_data_frame[1:n_genes,]
pao_gene_data_frame$numbers <- 1:n_genes
ggplot(data=pao_gene_data_frame, aes(x=numbers, y=gr1_promiscuity)) +
  geom_bar(stat="identity") + ggtitle('Gene Promiscuity between GR1 Sets') + labs(y='Number of Sets', x ='Gene')


# ggplot(data=pao_gene_data_frame, aes(x=reorder(genes, -gr1_promiscuity), y=gr1_promiscuity)) +
#   geom_bar(stat="identity") + ggtitle('Gene Promiscuity between GR1 Sets') + labs(y='Number of Sets', x ='Gene') + 
#   theme(#axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank())

load("~/GitHub/PathwayMining/scripts/pao/pao_g0_g1_data_frames.RData")

gr0_df <- c()# data.frame()
gr0_df$size <- vapply(pao_gr0_sets, function(x){length(x)}, c(1))
gr0_df <- data.frame(gr0_df)
gr0_df <- gr0_df[which(gr0_df$size > 0),]
gr0_df <- gr0_df[rev(order(gr0_df$size)),]
gr0_df$numbers <- 1:nrow(gr0_df)
ggplot(data=gr0_df, aes(x=numbers, y=size)) +
  geom_bar(stat="identity") + ggtitle('GR0 Set Size Distribution') + labs(y='Number of Genes in Set', x ='Set')

g0_df$X..genes <- unlist(g0_df$X..genes)
g0_df <- g0_df[which(g0_df$X..genes > 0),]
g0_df <- g0_df[rev(order(g0_df$X..genes)),]
g0_df$numbers <- 1:nrow(g0_df)
# g0_sets <- g0_df$clean.sets
ggplot(data=g0_df, aes(x=numbers, y=X..genes)) +
  geom_bar(stat="identity") + ggtitle('G0 Set Size Distribution') + labs(y='Number of Genes in Set', x ='Set')

