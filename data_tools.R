## data tools

# combine two tables
combine_tables <- function(table_1, table_2, iter_col_1 = 1, iter_col_2 = 1){ # iterate down table 1, specify row 
  aligned_mtx <- matrix(data = NA, nrow = nrow(table_1), ncol = ncol(table_2))
  
  for (i in 1:nrow(table_1)){
    it <- table_1[i, iter_col_1]
    pos <- which(table_2[, iter_col_2] == it)
    if (length(pos) > 0){
      aligned_mtx[i,] <- table_2[pos,]
    }
  }
  
  return(cbind(table_1, aligned_mtx[,-c(iter_col_2)]))
}

