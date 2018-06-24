library(gplots)
library(RColorBrewer)
source('~/GitHub/PathwayMining/set_tools.R')
set_cor <- function(C, set, avg = TRUE){
  if (length(set) < 2){
    return(0)
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

g1_gr1_compar <- function(C, g1_set, gr1_set){
  new_gr1 <- gr1_set[!gr1_set %in% g1_set]
  # print(new_gr0)
  genes <- c(g1_set, new_gr1)

  return(C[genes,genes])
}

g1_gr1_compar_wrapper <- function(idx_1, idx_2){
  return(g0_gr0_compar(C, g1[[idx_1]], gr1[[idx_2]]))
}

g1_gr1_sets <- function(g1_idx){
  print(g1[[g1_idx]])

  first_rxn <- g1[[g1_idx]][1]
  gr1_idx <- get_set_idx(first_rxn, gr1)[1]
  print(gr1[[gr1_idx]])

  return(c(g1_idx, gr1_idx))
}


plot_correlation_and_quantiles <- function(df, title = '', x = '', y = ''){
  df <- df[rev(order(df$avg_set_cor)),]
  df <- df[which(df$avg_set_cor != 0),]
  len <- 1:nrow(df)
  # plot(len, df$avg_set_cor, type = 'l', main = title, xlab = x, ylab = y, col = 'red')
  plot(0, type = 'n', xlim = c(0, nrow(df)), ylim = c(-0.5, 1), , main = title, xlab = x, ylab = y)
  for (i in len){
    rect(i-0.5, df$quantiles[[i]][2], i-0.5, df$quantiles[[i]][4], col = 'blue', border = 'blue')
  }
  lines(len, df$avg_set_cor, type = 'l', col = 'red')
}

set_analysis <- function(sets){
  avg_set_cor <- vapply(1:length(sets), function(x){return(set_cor(C, sets[[x]]))}, FUN.VALUE = c(1))
  avg_set_cor <- avg_set_cor
  set_quantile <- lapply(1:length(sets), function(x){return(quantile(set_cor(C, sets[[x]], avg = FALSE)))})
  set_analysis <- data.frame(avg_set_cor)
  set_analysis$quantiles <- set_quantile

  return(set_analysis)
}

g0_avg_set_cor <- vapply(1:length(g0), function(x){return(set_cor(C, g0[[x]]))}, FUN.VALUE = c(1))
gr0_avg_set_cor <- vapply(1:length(gr0), function(x){return(set_cor(C, gr0[[x]]))}, FUN.VALUE = c(1))
g1_avg_set_cor <- vapply(1:length(g1), function(x){return(set_cor(C, g1[[x]]))}, FUN.VALUE = c(1))
gr1_avg_set_cor <- vapply(1:length(gr1), function(x){return(set_cor(C, gr1[[x]]))}, FUN.VALUE = c(1))
# avg_set_cor <- g0_avg_set_cor
# g0_set_quantile <- lapply(1:length(g0), function(x){return(quantile(set_cor(C, g0[[x]], avg = FALSE)))})
# g0_set_analysis <- data.frame(avg_set_cor)
# g0_set_analysis$quantiles <- g0_set_quantile
#
# plot_correlation_and_quantiles(g0_set_analysis, 'G0 Average Set Correlation and Quartiles', 'Set', 'Correlation')

g0_set_analysis <- set_analysis(g0)
plot_correlation_and_quantiles(g0_set_analysis, 'G0 Average Set Correlation and Quartiles', 'Set', 'Correlation')

gr0_set_analysis <- set_analysis(gr0)
plot_correlation_and_quantiles(gr0_set_analysis, 'GR0 Average Set Correlation and Quartiles', 'Set', 'Correlation')

g1_set_analysis <- set_analysis(g1)
plot_correlation_and_quantiles(g1_set_analysis, 'G1 Average Set Correlation and Quartiles', 'Set', 'Correlation')

gr1_set_analysis <- set_analysis(gr1)
plot_correlation_and_quantiles(gr1_set_analysis, 'GR1 Average Set Correlation and Quartiles', 'Set', 'Correlation')

# avoid <- which(is.na(avg_set_cor_list))
# avg_set_cor_list[-c(avoid)]
# plot(1:length(g0_avg_set_cor), g0_avg_set_cor, main = 'Average Correlations Within G0 Sets', xlab = 'Set Index', ylab = 'Correlation')
# hist(g0_avg_set_cor, breaks = 20, main = 'Average Correlations Within G0 Sets', xlab = 'Correlation', ylab = 'Number of Sets')

# gr0_avg_set_cor <- vapply(1:length(gr0), function(x){return(set_cor(C, gr0[[x]]))}, FUN.VALUE = c(1))
# avg_set_cor <- gr0_avg_set_cor
# gr0_set_quantile <- lapply(1:length(gr0), function(x){return(quantile(set_cor(C, gr0[[x]], avg = FALSE)))})
# gr0_set_analysis <- data.frame(avg_set_cor)
# gr0_set_analysis$quantiles <- gr0_set_quantile
#
# plot_correlation_and_quantiles(gr0_set_analysis, 'GR0 Average Set Correlation and Quartiles', 'Set', 'Correlation')


# plot(1:length(gr0_avg_set_cor), gr0_avg_set_cor)
# hist(gr0_avg_set_cor, breaks = 20, main = 'Average Correlations Within GR0 Sets', xlab = 'Correlation', ylab = 'Number of Sets')
# compar <- g0_gr0_compar_wrapper(393, 104)
# #g0_gr0_compar(C, g0[[898]], gr0[[451]])
# compar2 <- compar[nrow(compar):1,]
# heatmap(compar2,Rowv=NA, Colv=NA, col=brewer.pal(9,"Blues"))

draw_heatmap_g0 <- function(g0_idx){
  idxs <- g0_gr0_sets(g0_idx)
  compar <- g0_gr0_compar_wrapper(idxs[1], idxs[2])
  cols <- rep('black', nrow(compar))
  cols[row.names(compar) %in% g0[[idxs[1]]]] <- 'red'
  # colors = c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100))
  #
  # my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 300)

  #g0_gr0_compar(C, g0[[898]], gr0[[451]])
  # compar2 <- compar[nrow(compar):1,]
  heatmap.2(compar,Rowv=FALSE,Colv=FALSE, colRow = cols, colCol = cols, col=rev(bluered(256)), symbreaks = TRUE, symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')
  # heatmap.2(compar,Rowv=FALSE,Colv=FALSE, col=brewer.pal(9,"Blues"), symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')
  # heatmap(compar2,Rowv=NA, Colv=NA, col=brewer.pal(9,"Blues"))
}

draw_heatmap_g1 <- function(g1_idx){
  idxs <- g1_gr1_sets(g1_idx)
  compar <- g1_gr1_compar_wrapper(idxs[1], idxs[2])
  # print(paste('compar:', rownames(compar)))
  col_list <- c('brown3', 'cadetblue3','chartreuse3','chocolate3','burlywood3','darkgoldenrod4','aquamarine3','darkorchid3','deeppink','lightgoldenrod4',
                'midnightblue','mistyrose3','lightsteelblue2','lightsalmon','yellowgreen','tan1','springgreen3','snow3','yellow3','slateblue','sienna3',
                'violetred3','seagreen4','royalblue2','palevioletred','paleturquoise3','seagreen3','plum4')
  row_cols <- rep('white', nrow(compar))
  col_cols <- rep('white', nrow(compar))
  for (i in 1:length(rownames(compar))){
    # print(rownames(compar)[i])
    g1_set <- g1[[get_set_idx(rownames(compar)[i], g1)]]
    # print(row.names(compar) %in% g1_set)
    # col_idx <- i%%length(col_list)+1
    # row_cols[row.names(compar) %in% g1_set] <- col_list[col_idx]
    if (length(which(g1_set %in% row.names(compar))) > 1){
      # print(i)
      col_idx <- i%%length(col_list)+1
      row_cols[row.names(compar) %in% g1_set] <- col_list[col_idx]
      # col_list <- col_list[-col_idx]
      # print(col_list)
    }
  }
  for (i in 1:length(rownames(compar))){
    # print(rownames(compar)[i])
    g0_set <- g0[[get_set_idx(rownames(compar)[i], g0)]]
    # print(row.names(compar) %in% g1_set)
    # col_idx <- i%%length(col_list)+1
    # col_cols[row.names(compar) %in% g0_set] <- col_list[col_idx]
    if (length(which(g0_set %in% row.names(compar))) > 1){
      # print(i)

      col_idx <- i%%length(col_list)+1
      col_cols[row.names(compar) %in% g0_set] <- col_list[col_idx]
      # col_list <- col_list[-col_idx]
      # print(col_list)
    }
  }
  # print(cols)

  # colors = c(seq(-1,-0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100))

  # my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 300)

  #g0_gr0_compar(C, g0[[898]], gr0[[451]])
  # compar2 <- compar[nrow(compar):1,]
  heatmap.2(compar,Rowv=FALSE,Colv=FALSE, RowSideColors = row_cols,ColSideColors = col_cols, col=rev(bluered(256)), symbreaks = TRUE, symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')
  # heatmap.2(compar,Rowv=FALSE,Colv=FALSE, colRow = cols, colCol = cols, col=rev(bluered(256)), symbreaks = TRUE, symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')
  # heatmap.2(compar,Rowv=FALSE,Colv=FALSE, col=brewer.pal(9,"Blues"), symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')
  # heatmap(compar2,Rowv=NA, Colv=NA, col=brewer.pal(9,"Blues"))
}

draw_heatmap_g1_99 <- function(){
  idxs <- g1_gr1_sets(99)
  compar <- g1_gr1_compar_wrapper(idxs[1], idxs[2])
  cols <- rep('white', nrow(compar))
  cols[row.names(compar) %in% g1[[idxs[1]]]] <- 'red'
  cols[row.names(compar) %in% g1[[142]]] <- 'green'
  # side_colors <-

  heatmap.2(compar,Rowv=FALSE,Colv=FALSE, RowSideColors = cols, col=rev(bluered(256)), symbreaks = TRUE, symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none')

}
draw_heatmap_g1_99()
# interesting idxs g0: 5, 211, 91, 18, 19, 20, 37, 42, 55, 60, 44
# interesting idxs g1: 39, 99?. 110 (same as g0 set 211), 352, 375?, 30, 38?, 108, 405?, 17, 28, 31, 95

draw_heatmap_g0(384)
# compar_df <- data.frame(compar)
# ggplot(data.frame(compar)) + geom_tile(aes(fill = rescale), + colour = "white") + scale_fill_gradient(low = "white", +high = "steelblue")
# print(g0_gr0_compar(C, g0[[898]], gr0[[451]]))
which(g0_avg_set_cor < 0.6 & g0_avg_set_cor > 0.5)

which(g1_avg_set_cor < 0.6 & g1_avg_set_cor > 0.5)
