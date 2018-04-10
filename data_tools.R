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

read_coupling_csv <- function(filename, coupling_vector = c(), completed_idxs = c()){
  lines <- readLines(filename)
  n_lines <- length(lines)

  for (line in lines){
    nums <- as.numeric(strsplit(line, split = ',')[[1]])
    idx <- nums[1]
    if (idx %in% completed_idxs){next}
    coupling_idxs <- nums[2:length(nums)]

    completed_idxs <- c(completed_idxs, idx)
    coupling_vector[[idx]] <- coupling_idxs
  }


  list(
    coupling_vector = coupling_vector,
    completed_idxs = completed_idxs
  )
}

coupling_vector <- read_coupling_csv('~/GitHub/PathwayMining/scripts/pao_coupling.csv')
