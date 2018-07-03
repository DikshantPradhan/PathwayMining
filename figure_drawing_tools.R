library(ggplot2)
load("~/GitHub/PathwayMining/data/pao_model/pao_model.RData")

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

qplot(pao_gene_data_frame$gpr_promiscuity, geom="histogram", binwidth = 1) +
  xlab('Number of Associated Reactions') + ylab('Number of Genes') +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

qplot(pao_gene_data_frame$gr0_promiscuity[which(pao_gene_data_frame$gr0_promiscuity > 0)], geom="histogram", binwidth = 1, xlim = c(0,15)) +
  xlab('Number of Associated Sets') + ylab('Number of Genes') +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

# ggplot(data=pao_gene_data_frame, aes(pao_gene_data_frame$gpr_promiscuity)) + 
#   geom_histogram(stat="identity")# + labs(y='Number of Reactions', x ='Gene') #Gene Promiscuity between Reactions
#   # geom_bar(aes(fill = cut))+
#   # theme(#axis.title.x=element_blank(),
#   #       axis.text.x=element_blank(),
#   #       axis.ticks.x=element_blank())

pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gr0_promiscuity)),]
pao_gene_data_frame <- pao_gene_data_frame[1:n_genes,]
pao_gene_data_frame$numbers <- 1:n_genes
ggplot(data=pao_gene_data_frame, aes(x=numbers, y=gr0_promiscuity)) +
  geom_bar(stat="identity") + labs(y='Number of Sets', x ='Gene') # Gene Promiscuity between GR0 Sets

pao_gene_data_frame <- pao_gene_data_frame[rev(order(pao_gene_data_frame$gr1_promiscuity)),]
pao_gene_data_frame <- pao_gene_data_frame[1:n_genes,]
pao_gene_data_frame$numbers <- 1:n_genes
ggplot(data=pao_gene_data_frame, aes(x=numbers, y=gr1_promiscuity)) +
  geom_bar(stat="identity") + labs(y='Number of Sets', x ='Gene') # Gene Promiscuity between GR1 Sets


# ggplot(data=pao_gene_data_frame, aes(x=reorder(genes, -gr1_promiscuity), y=gr1_promiscuity)) +
#   geom_bar(stat="identity") + ggtitle('Gene Promiscuity between GR1 Sets') + labs(y='Number of Sets', x ='Gene') + 
#   theme(#axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank())

# load("~/GitHub/PathwayMining/scripts/pao/pao_g0_g1_data_frames.RData")
load("~/GitHub/PathwayMining/EnzymeCoupling/data/pao_g0_g1_data_frames.RData")

gr0_df <- 1:length(pao_gr0_sets) #c()# data.frame()
gr0_df$size <- vapply(pao_gr0_sets, function(x){length(x)}, c(1))
gr0_df <- data.frame(gr0_df)
gr0_df <- gr0_df[which(gr0_df$size > 0),]
gr0_df <- gr0_df[rev(order(gr0_df$size)),]
gr0_df$numbers <- 1:nrow(gr0_df)
ggplot(data=gr0_df, aes(x=numbers, y=size)) +
  geom_bar(stat="identity") + labs(y='Number of Genes in Set', x ='Set') # GR0 Set Size Distribution

gr1_df <- 1:length(pao_gr1_sets) #c()# data.frame()
gr1_df$size <- vapply(pao_gr1_sets, function(x){length(x)}, c(1))
gr1_df <- data.frame(gr1_df)
gr1_df <- gr1_df[which(gr1_df$size > 0),]
gr1_df <- gr1_df[rev(order(gr1_df$size)),]
gr1_df$numbers <- 1:nrow(gr1_df)

g0_df$X..genes <- unlist(g0_df$X..genes)
g0_df <- g0_df[which(g0_df$X..genes > 0),]
g0_df <- g0_df[rev(order(g0_df$X..genes)),]
g0_df$numbers <- 1:nrow(g0_df)
# g0_sets <- g0_df$clean.sets
ggplot(data=g0_df, aes(x=numbers, y=X..genes)) +
  geom_bar(stat="identity") + labs(y='Number of Genes in Set', x ='Set') # G0 Set Size Distribution

ggplot(data=set_df, aes(y=g0_sets)) +
  geom_bar(stat="identity") + labs(y='Number of Genes in Set', x ='Set') # G0 Set Size Distribution

gr0_set_sizes <- gr0_df$size[which(gr0_df$size > 0)]
g0_set_sizes <- g0_df$X..genes[which(g0_df$X..genes > 0)]

gr1_set_sizes <- gr1_df$size[which(gr1_df$size > 0)]
g1_set_sizes <- unlist(g1_df$X..genes[which(g1_df$X..genes > 0)])


# set_df <- data.frame(sets = c(g0_set_sizes, gr0_set_sizes), category = c(rep('G0 Sets', length(g0_set_sizes)), rep('GR0 Sets', length(gr0_set_sizes))))
set_df <- data.frame(sets = c(g1_set_sizes, gr1_set_sizes), category = c(rep('G1 Sets', length(g1_set_sizes)), rep('GR1 Sets', length(gr1_set_sizes))))
ggplot(set_df, aes(x = sets, fill = category)) + xlim(c(8,125)) + #geom_histogram(alpha = 0.3) + xlab('Set Size') + ylab('Count') + 
  geom_histogram(position = "dodge", alpha = .8) + 
  # scale_fill_manual(values=c("blue", "red")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
  #scale_y_log10()# + scale_fill_manual(values=c("grey20", "grey60")) 
# ggplot(set_df, aes(sets)) + 
#   geom_histogram(data = g0_sets, fill = "red", alpha = 0.2) + 
#   geom_histogram(data = gr0_sets, fill = "blue", alpha = 0.2)
