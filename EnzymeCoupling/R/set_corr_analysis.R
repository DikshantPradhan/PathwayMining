library(gplots)
library(RColorBrewer)
set_cor <- function(C, set, avg = TRUE){
  if (length(set) < 2){
    return(NA)
  }

  cor_list <- c()

  for (i in 1:(length(set)-1)){
    for (j in (i+1):length(set)){
      # print(paste(set[i],set[j]))
      if (set[i] %in% rownames(C) & set[j] %in% rownames(C)){
        cor_list <- c(cor_list, C[set[i],set[j]])
      }
    }
  }

  if (!avg){
    return(cor_list)
  }

  return(mean(cor_list))
}

g0_gr0_compar <- function(C, g0_set, gr0_set){
  new_gr0 <- gr0_set[!gr0_set %in% g0_set]
  # print(new_gr0)
  genes <- c(g0_set, new_gr0)

  return(C[genes,genes])
}

g0_gr0_compar_wrapper <- function(idx_1, idx_2){
  return(g0_gr0_compar(C, g0[[idx_1]], gr0[[idx_2]]))
}

g0_gr0_sets <- function(g0_idx){
  print(g0[[g0_idx]])

  first_rxn <- g0[[g0_idx]][1]
  gr0_idx <- get_set_idx(first_rxn, gr0)[1]
  print(gr0[[gr0_idx]])

  return(c(g0_idx, gr0_idx))
}

g0_avg_set_cor <- vapply(1:length(g0), function(x){return(set_cor(C, g0[[x]]))}, FUN.VALUE = c(1))
# avoid <- which(is.na(avg_set_cor_list))
# avg_set_cor_list[-c(avoid)]
plot(1:length(g0_avg_set_cor), g0_avg_set_cor, main = 'Average Correlations Within G0 Sets', xlab = 'Set Index', ylab = 'Correlation')
hist(g0_avg_set_cor, breaks = 20, main = 'Average Correlations Within G0 Sets', xlab = 'Correlation', ylab = 'Number of Sets')

gr0_avg_set_cor <- vapply(1:length(gr0), function(x){return(set_cor(C, gr0[[x]]))}, FUN.VALUE = c(1))
plot(1:length(gr0_avg_set_cor), gr0_avg_set_cor)
hist(gr0_avg_set_cor, breaks = 20, main = 'Average Correlations Within GR0 Sets', xlab = 'Correlation', ylab = 'Number of Sets')
# compar <- g0_gr0_compar_wrapper(393, 104)
# #g0_gr0_compar(C, g0[[898]], gr0[[451]])
# compar2 <- compar[nrow(compar):1,]
# heatmap(compar2,Rowv=NA, Colv=NA, col=brewer.pal(9,"Blues"))

draw_heatmap <- function(g0_idx){
  idxs <- g0_gr0_sets(g0_idx)
  compar <- g0_gr0_compar_wrapper(idxs[1], idxs[2])
  #g0_gr0_compar(C, g0[[898]], gr0[[451]])
  # compar2 <- compar[nrow(compar):1,]
  heatmap.2(compar,Rowv=FALSE,Colv=FALSE, col=brewer.pal(9,"Blues"), symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')
  # heatmap(compar2,Rowv=NA, Colv=NA, col=brewer.pal(9,"Blues"))
}

# interesting idxs: 393, 395, 405, 382, 384, 410
# : 547, 403, 645, 849, 890, 896
# : 628?
draw_heatmap(384)
# compar_df <- data.frame(compar)
# ggplot(data.frame(compar)) + geom_tile(aes(fill = rescale), + colour = "white") + scale_fill_gradient(low = "white", +high = "steelblue")
# print(g0_gr0_compar(C, g0[[898]], gr0[[451]]))
which(g0_avg_set_cor < 0.6 & g0_avg_set_cor > 0.5)
