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

pao_gene_data_frame <- gene_data_frame(pao_model, gr0_sets, gr1_sets)
# pao_gene_data_frame <- data.frame(pao_gene_data_frame)
# pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gpr_promiscuity)),]
ggplot(data=pao_gene_data_frame, aes(x=reorder(genes, -gpr_promiscuity), y=gpr_promiscuity)) +
  geom_bar(stat="identity")

ggplot(data=pao_gene_data_frame, aes(x=reorder(genes, -gr0_promiscuity), y=gr0_promiscuity)) +
  geom_bar(stat="identity")

ggplot(data=pao_gene_data_frame, aes(x=reorder(genes, -gr1_promiscuity), y=gr1_promiscuity)) +
  geom_bar(stat="identity")
