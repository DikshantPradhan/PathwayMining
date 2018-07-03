library(ggplot2)
library(dplyr)

gene_set_dataframe <- function(sets){
  # clean_sets <- clean_rxn_names_in_set(sets)
  num_genes <- lapply(sets, length)
  
  frame <- cbind(sets, sets, num_genes)
  colnames(frame) <- c('sets', 'clean sets', '# genes')
  return(frame)
}

map_genes_to_sets <- function(genes, set){
  map <- matrix(data = 0, nrow = length(genes), ncol = 1)
  for (i in 1:length(map)){
    pos <- grep(genes[i], set)
    if (length(pos) < 1){next}
    # print(pos)
    map[i] <- pos
  }
  return(map)
}

binary_set_assignment <- function(sets, assignment){
  
  # new_sets <- c()
  
  for (i in 1:length(sets)){
    # new_set <- c()
    if (is.null(sets[[i]])){next()}
    for (j in 1:length(sets[[i]])){
      if (sets[[i]][j] %in% assignment){
        sets[[i]][j] <- 1
        # new_set[k] <- clean_ex_a(set_list[[i]][j])
        # k <- k+1
      }
      else {
        sets[[i]][j] <- 0
        # next
      }
    }
    # if (length(new_set) > 0){
    #   new_sets[i] <- list(new_set)
    # }
  }
  
  return(sets)
}

binary_set_histogram <- function(binary_sets){
  assignment_fraction <- matrix(data = -1, nrow = length(binary_sets), ncol = 1)
  for (i in 1:length(binary_sets)){
    assignment_fraction[i] <- length(which(binary_sets[[i]] == 1))/length(binary_sets[[i]])
  }
  return(assignment_fraction)
}

double_histogram <- function(data1, data2, title){
  dat1 = data.frame(x = as.numeric(data1), group="observed")
  dat2 = data.frame(x = as.numeric(data2), group="bootstrapped")
  dat = rbind(dat1, dat2)
  ggplot(dat, aes(x, fill=group, colour=group)) +
    geom_histogram(aes(y=..density..), breaks=seq(0,1,0.05), alpha=0.6, 
                   position="identity", lwd=0.2) +
    ggtitle(title)
}

gene_set_expanded_dataframe <- function(set_df, sets, df_genes, de_genes, interest_genes, pure_df_genes, pure_de_genes, dfde_genes){
  df <- data.frame(set_df)
  df$df <- binary_set_assignment(sets, df_genes)
  df$de <- binary_set_assignment(sets, de_genes)
  df$interest_genes <- binary_set_assignment(sets, interest_genes)
  df$pure_df <- binary_set_assignment(sets, pure_df_genes)
  df$pure_de <- binary_set_assignment(sets, pure_de_genes)
  df$dfde <- binary_set_assignment(sets, dfde_genes)
  df$df_frac <- binary_set_histogram(df$df)
  df$de_frac <- binary_set_histogram(df$de)
  df$pure_df_frac <- binary_set_histogram(df$pure_df)
  df$pure_de_frac <- binary_set_histogram(df$pure_de)
  df$dfde_frac <- binary_set_histogram(df$dfde)
  df$interest_frac <- binary_set_histogram(df$interest_genes)
  
  return(df)
}


plot_density <- function(mtx, xlimits=c(0, 1), ylimits=c(0, 200), start = 2){
  n <- ncol(mtx)
  plot(1, type="n", xlab="", ylab="", xlim=xlimits, ylim=ylimits)
  # freqs <- table(as.numeric(mtx[,2]))
  # freqs_x <- as.numeric(names(freqs))
  # plot(freqs_x, freqs)
  for (i in start:n){
    freqs <- table(as.numeric(mtx[,i]))
    freqs_x <- as.numeric(names(freqs))
    lines(freqs_x, freqs, col = 'grey')
    # line(density(mtx[,i]))
  }
}

plot_set_characteristic_hypothesis <- function(observed, simulated, start = 2, end = 1001, title = '', color = 'lightskyblue'){
  obs_frac <- length(which(observed == 0 | observed == 1))/length(observed)
  sim_fracs <- c()
  
  for (i in start:end){
    sim <- length(which(simulated[,i] == 0 | simulated[,i] == 1))/nrow(simulated)
    # print(i)
    sim_fracs <- c(sim_fracs, sim)
  }
  sim_fracs <- data.frame(sim_fracs)
  colnames(sim_fracs) <- 'frac'
  # print(length(sim_fracs))
  # hist(sim_fracs, breaks = 50)
  # line(obs_frac)
  ggplot(sim_fracs, aes(x=frac)) + #ggtitle(title) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.01,
                   colour="black", fill=color) +
    geom_density() + scale_x_continuous(limits = c(0,1)) +
    xlab('Complete Set Fraction') + ylab('Density') +
    geom_vline(xintercept = obs_frac, colour = 'red', linetype="dashed") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
    # guide_legend('test')
    # geom_histogram(binwidth=.01, colour="black", fill="white")
  # ggplot(sim_fracs, aes(x = frac)) + geom_density() +
  # qplot(sim_fracs, geom = 'histogram', bins = 50, main = title, xlim = c(0,1.1)) + geom_density() +
  # qplot(sim_fracs, main = title, xlim = c(0,1.1)) + geom_density() +
    # geom_vline(xintercept = obs_frac, colour = 'red') + xlab('Complete Set Fraction') + ylab('count')
}

plot_genes_barplot <- function(sets, pure_df, pure_de, dfde_list, genes_of_interest, percent = FALSE){
  genes <- unique(unlist(sets))
  n_genes <- length(genes)
  n_sets <- length(sets)
  
  set_idxs <- matrix(data = 'none', nrow = n_genes, ncol = 1)
  category <- matrix(data = 'none', nrow = n_genes, ncol = 1)
  
  for (i in 1:n_genes){
    idx <- get_set_idx(genes[i], sets)
    # print(idx)
    set_idxs[i] <- idx
    set_pos <- which(sets[[idx]] == genes[i])
    # print(set_pos)
    df <- pure_df[[idx]]
    de <- pure_de[[idx]]
    dfde <- dfde_list[[idx]]
    
    # print(dfde)
    
    if (as.numeric(df[set_pos]) == 1){category[i] <- 'df'}
    if (as.numeric(de[set_pos]) == 1){category[i] <- 'de'}
    if (as.numeric(dfde[set_pos]) == 1){category[i] <- 'dfde'}
    
    # else {
    #     if (de[set_pos] == 1){
    #       category[i] == 'de'
    #     } else {
    #         if (dfde[set_pos] == 1){category[i] == 'dfde'}
    #     }
    #   }
  }
  data <- data.frame(genes = genes, set_idx = set_idxs, category = category)
  # print(length(which(data$category != 'none')))
  # print(length(which(category != 'none')))
  print(head(data))
  set_idxs <- rep(1:n_sets, times=1, each=4)
  cats <- c('df','dfde','de','none')
  set_category <- rep(cats, times = n_sets, each = 1)#matrix(data = 'na', nrow = n_sets*4, ncol = 4)
  
  print(head(set_idxs))
  print(head(set_category))
  
  value <- matrix(data = 0, nrow = 4*n_sets, ncol = 1)
  for (i in 1:length(set_idxs)){
    value[i] <- length(which(data$set_idx == as.numeric(set_idxs[i]) & data$category == set_category[i]))
  }
  
  bar_plot <- data.frame(set_idx = set_idxs, category = set_category, value = value)
  bar_plot$category <- factor(bar_plot$category, levels = rev(cats), ordered = TRUE) 
  
  # print(factor(bar_plot$category, levels=c("df", "dfde", "de", "none")))#bar_plot[order(bar_plot$category, leve),]
  # print(head(bar_plot))
  
  if (percent){
    # bar_plot$category <- factor(bar_plot$category, levels = rev(cats), ordered = TRUE)
    ggplot() + geom_bar(data=bar_plot, aes(fill=category, y=value, x=set_idx), stat="identity", position = 'fill') +
      labs(y='percent') +
      coord_flip() + ggtitle('Gene Distribution Percentages in Sets')
  }
  else {
    bar_plot_none <- bar_plot
    bar_plot_none$value[which(bar_plot$category != 'none')] <- 0
    bar_plot_none$value <- -1*bar_plot_none$value
    bar_plot$value[which(bar_plot$category == 'none')] <- 0
    
    ggplot() + geom_bar(data=bar_plot, aes(fill=category, y=value, x=set_idx), stat="identity") + labs(y = 'frequency') +
      geom_bar(data=bar_plot_none, aes(fill=category, y=value, x=set_idx), stat="identity") +
      coord_flip() + ggtitle('Gene Distributions in Sets')
  }
}

# shuffle whole sets and assign df genes and dde genes to complete sets and observe overlap to assign dfde
sample_dfde_sets <- function(sets, n_df, n_de){
  
  n_sets <- length(sets)
  all_genes <- unique(unlist(sets))
  n_genes <- length(unique(unlist(sets)))
  
  sample_sets <- function(sets, n){
    # shuffle sets
    set_order <- sample(1:n_sets, size = n_sets)
    sets <- sets[set_order]
    
    assignment <- rep(FALSE, n_genes)
    names(assignment) <- all_genes
    
    # set assignments
    remaining <- n
    set_idx <- 1
    while (remaining > 0){
      genes <- sets[[set_idx]]
      
      if (length(genes) > remaining){
        genes <- sample(genes, size = remaining)
        assignment[genes] <- TRUE
        remaining_df <- 0
      }
      
      assignment[genes] <- TRUE
      remaining <- remaining - length(genes)
      set_idx <- set_idx+1
    }
    
    return(assignment)
  }
  
  df_assignment <- sample_sets(sets, n_df)
  de_assignment <- sample_sets(sets, n_de)
  
  dfde_genes <- names(df_assignment[which(df_assignment & de_assignment)])
  
  return(dfde_genes)
}

plot_dfde_sampling <- function(sets, num_df, num_de, num_dfde){
  # g0_dfde_obs <- g0_df$dfde[which(g0_set_df[,3] > 1)]
  # 
  # binary_df_g1_sets <- binary_set_assignment(g1_sets, df_genes)
  # binary_df_g1_assign <- binary_set_histogram(binary_df_g1_sets)
  # binary_de_g1_sets <- binary_set_assignment(g1_sets, de_genes)
  # binary_de_g1_assign <- binary_set_histogram(binary_de_g1_sets)
  # g1_set_df <- cbind(g1_set_df, binary_df_g1_assign, binary_de_g1_assign)
  # colnames(g1_set_df) <- c('sets', 'clean sets', '# genes', 'df frac', 'de frac')
  # g1_df_obs <- g1_set_df[which(g1_set_df[,3] > 1),4]
  # g1_de_obs <- g1_set_df[which(g1_set_df[,3] > 1),5]
  # 
  # num_g0_sets <- length(which(g0_set_df[,3] > 1))
  # num_g1_sets <- length(which(g1_set_df[,3] > 1))
  # 
  n = 1000
  # g0_dfde_vec <- matrix(data = 0, nrow = num_g0_sets, ncol = 1)
  
  num_sampled_dfde <- c()
  for (i in 1:n){
    dfde_genes <- sample_dfde_sets(sets, num_df, num_de)
    num_sampled_dfde <- c(num_sampled_dfde, length(dfde_genes))
    # binary_dfde_g0_sets <- binary_set_assignment(g0_sets, dfde_genes)
    # binary_dfde_g0_assign <- binary_set_histogram(binary_dfde_g0_sets)
    # g0_set_df <- cbind(g0_set_df_og, binary_dfde_g0_assign)
    # g0_dfde_vec <- cbind(g0_df_vec, g0_set_df[which(g0_set_df[,3] > 1),4])
    
    
  }
  # print(num_sampled_dfde)
  
  qplot(num_sampled_dfde, geom = 'histogram', bins = 50, main = 'dFdE gene bootstrap') + 
    geom_vline(xintercept = num_dfde, colour = 'red') + xlab('number of dfde genes') + ylab('frequency')
  
  
  # plot_density(g0_dfde_vec)
  # freqs <- table(as.numeric(g0_dfde_obs))
  # freqs_x <- as.numeric(names(freqs))
  # lines(freqs_x, freqs, col = 'red')
  # title(main="g0 dFdE set composition", xlab="% of genes", ylab="frequency")
  
}

# shuffle df and de fractions and observe overlap to assign dfde
bootstrap_dfde <- function(df_frac_list, de_frac_list, set_size_list, fracs = TRUE){
  if (length(df_frac_list) != length(de_frac_list)){return()}
  n_sets <- length(df_frac_list)
  
  # df_binary_sets <- empty_binary_sets
  # de_binary_sets <- empty_binary_sets
  
  dfde_sample <- matrix(data = 0, ncol = 1, nrow = n_sets)
  
  for (i in 1:n_sets){
    df_frac <- df_frac_list[i]
    de_frac <- de_frac_list[i]
    
    set_size <- set_size_list[[i]]
    # print(set_size)
    
    df_sampled_idxs <- sample(1:set_size, size = floor(set_size*df_frac))
    de_sampled_idxs <- sample(1:set_size, size = floor(set_size*de_frac))
    
    # df_binary_sets[[i]][df_sampled_idxs] <- 1
    # de_binary_sets[[i]][de_sampled_idxs] <- 1
    if (fracs){
      dfde_sample[i] <- length(intersect(df_sampled_idxs, de_sampled_idxs))/set_size
    }
    else {
      dfde_sample[i] <- length(intersect(df_sampled_idxs, de_sampled_idxs))
    }
  }
  # print(dim(dfde_frac_list))
  return(dfde_sample)
}

plot_dfde_bootstrap <- function(df, num_dfde){
  reduced_df <- df[which(df$X..genes > 1),]
  n_sets <- nrow(df)

  n = 1000
  dfde_vec <- matrix(data = 0, nrow = n_sets, ncol = 1)
  # num_sampled_dfde <- c()
  
  for (i in 1:n){
    df_frac_sampled <- sample(df$df_frac, size = n_sets)
    de_frac_sampled <- sample(df$de_frac, size = n_sets)
    # print(paste(length(g0_df_frac_sampled), length(g0_de_frac_sampled)))
    dfde_sample <- bootstrap_dfde(df_frac_sampled, de_frac_sampled, df$X..genes, fracs = TRUE)
    # print(dim(dfde_frac_list))
    # num_sampled_dfde <- c(num_sampled_dfde, sum(dfde_sample))
    dfde_vec <- cbind(dfde_vec, dfde_sample)
  }
  
  # qplot(num_sampled_dfde, geom = 'histogram', bins = 50, main = 'dFdE gene bootstrap') + 
  #   geom_vline(xintercept = num_dfde, colour = 'red') + xlab('number of dfde genes') + ylab('frequency')
  
  plot_density(dfde_vec, ylimits = c(0,110))
  freqs <- table(as.numeric(reduced_df$dfde_frac))
  freqs_x <- as.numeric(names(freqs))
  lines(freqs_x, freqs, col = 'red')
  title(main="dFdE set composition", xlab="% of genes with fitness defect & expression change", ylab="frequency")
  
}

ratio_vector <- function(categories, set_idxs){
  idxs <- 1:length(unique(set_idxs))
  ratios <- matrix(data = 0, nrow = length(unique(set_idxs)), ncol = 1)
  for (i in idxs){
    sets <- which(set_idxs == i)
    total <- length(set)
    
    interest <- length(grep('none', categories[sets]))
    ratio <- interest/total
    ratios[i] <- ratio
  }
  return(ratios)
}

de_df_mixing <- function(categories, g1_sets){
  de_genes <- grep('full', categories)
  mixed_genes <- matrix(data = FALSE, nrow = length(de_genes), ncol = 1)
  
  for (i in de_genes){
    g1_set <- g1_sets[i]
    mixed <- ('full_df' %in% categories[which(g1_sets == g1_set)]) | ('full_dfde' %in% categories[which(g1_sets == g1_set)])
    if (mixed){mixed_genes[i] <- TRUE}
  }
  
  return(mixed_genes)
}

g1_composition_sampling <- function(g0_categories, g1_set_comp){
  plot_idxs <- which(g1_set_comp > 1)
  g0_categories <- g0_categories[sample(x = 1:length(g0_categories), size = length(g0_categories))]
  value <- rep(x = 1, size = length(g0_categories))
  g1_set_idxs <- rep(1:length(g1_set_comp), times = g1_set_comp)
  return(ratio_vector(g0_categories, g1_set_idxs))
  
  # g1_composition_df <- data.frame(g1_set = g1_set_idxs, category = g0_categories, value = value)
  # g1_composition_df$category <- factor(g1_composition_df$category, levels = c(1,0), ordered = TRUE)
  # ggplot() + geom_bar(data=g1_composition_df[which(g1_set_idxs %in% plot_idxs),], aes(fill=category, y=value, x=g1_set), stat="identity") + labs(y = '# of G0 sets') +
  #   coord_flip() + ggtitle('G1 Set Composition')
  
  # de_mixed <- de_df_mixing(g0_categories, g1_set_idxs)
  # return(length(which(de_mixed)))
  # g1_ratios <- matrix(data = 0, nrow = length(plot_idxs), ncol = 1)
  # for (i in plot_idxs){
  #   set_idx <- which(g1_set_idxs==i)
  #   df <- as.integer(length(which(g0_categories[set_idx] == 'full_df')) >= 1)
  #   de <- as.integer(length(which(g0_categories[set_idx] == 'full_de')) >= 1)
  #   # dfde <- as.integer(length(which(g0_categories[set_idx] == 'full_dfde')) >= 1)
  #   total <- df+de
  #   # if (df == 0 | de == 0){next}
  #   if (total > 1){g1_ratios[i] <- 1}
  #   # ratio <- df/de
  #   # g1_ratios[i] <- 1 #ratio
  # }
  # return(sum(g1_ratios))
  # return(sum(g0_categories[plot_idxs])/sum(g0_categories))
}

g1_ratio_sampling <- function(g0_categories, g1_set_comp, n, original_g1_idxs){
  plot_idxs <- which(g1_set_comp > 1)
  # og_sum <- sum(g0_categories[plot_idxs])/sum(g0_categories)
  
  ratios <- ratios <- matrix(data = 0, nrow = length(unique(original_g1_idxs)), ncol = 1)  #c()
  for (i in 1:n){
    ratios <- cbind(ratios, g1_composition_sampling(g0_categories, g1_set_comp))
  }
  
  # og_de_mixed <- length(which(de_df_mixing(g0_categories, original_g1_idxs)))
  og_ratios <- ratio_vector(g0_categories, original_g1_idxs)
  plot_density(ratios, ylimits = c(0,110))
  freqs <- table(as.numeric(og_ratios))
  freqs_x <- as.numeric(names(freqs))
  lines(freqs_x, freqs, col = 'red')
  title(main="g0 dFdE set composition", xlab="% of genes with fitness defect & expression change", ylab="frequency")
  
  
  # print()
  # return(g1_ratios)
  # # print(g1_ratios)
  # sampled_ratios <- table(ratios)
  # observed_ratios <- table(g1_ratios)
  # 
  # hist(ratios)
  # qplot(ratios, geom = 'histogram', bins = 50, main = '') +
  #   geom_vline(xintercept = og_de_mixed, colour = 'red') + xlab('') + ylab('frequency')
  
  # qplot(ratios, geom = 'histogram', bins = 50, main = 'set composition ratios') +
  # #   qplot(g1_ratios, geom = 'histogram', bins = 50) +
  #   xlab('fraction of g0 sets in g1 sets which are full') + ylab('frequency') #+ scale_y_continuous(formatter="percent")

  # qplot(g1_ratios, geom = 'histogram', bins = 50, main = 'set composition ratios') +
  #   #   qplot(g1_ratios, geom = 'histogram', bins = 50) +
  #   xlab('fraction of g0 sets in g1 sets which are full') + ylab('frequency') #+ scale_y_continuous(formatter="percent")
}

g1_architecture_measurement <- function(coupling_mtx, df_frac, de_frac){
  edge_list <- matrix(nrow = length(which(coupling_mtx)), ncol = 3)
  set_idxs <- as.numeric(rownames(coupling_mtx))
  
  edge_idx <- 1
  for (i in 1:nrow(coupling_mtx)){
    for (j in which(coupling_mtx[i,])){
      set_idx_1 <- set_idxs[i]
      set_idx_2 <- set_idxs[j]
      df_diff <- abs(df_frac[set_idx_1] - df_frac[set_idx_2])
      de_diff <- abs(de_frac[set_idx_1] - de_frac[set_idx_2])
      euclidian_diff <- sqrt((df_frac[set_idx_1] - df_frac[set_idx_2])^2 + (de_frac[set_idx_1] - de_frac[set_idx_2])^2)
      
      # g1_set <- get_set_idx(set_idxs[i], g1_composition)
      # print(g1_set)
      if (is.nan(df_diff)){
        df_diff <- -0.1
        de_diff <- -0.1
        euclidian_diff <- -0.1
      }
      
      edge_list[edge_idx,1] <- df_diff
      edge_list[edge_idx,2] <- de_diff
      edge_list[edge_idx,3] <- euclidian_diff
      # edge_list[edge_idx,4] <- g1_set
      
      edge_idx <- edge_idx + 1
    }
  }
  
  return(edge_list)
}

g1_architecture_measurement_binary <- function(coupling_mtx, df_frac, de_frac){
  edge_list <- matrix(nrow = length(which(coupling_mtx)), ncol = 3)
  set_idxs <- as.numeric(rownames(coupling_mtx))
  
  edge_idx <- 1
  for (i in 1:nrow(coupling_mtx)){
    for (j in which(coupling_mtx[i,])){
      set_idx_1 <- set_idxs[i]
      set_idx_2 <- set_idxs[j]
      df_diff <- df_frac[set_idx_1]&df_frac[set_idx_2]
      de_diff <- de_frac[set_idx_1]&de_frac[set_idx_2]
      euclidian_diff <-(df_frac[set_idx_1]&de_frac[set_idx_2])|(df_frac[set_idx_2]&de_frac[set_idx_1]) #((df_frac[set_idx_1]&!de_frac[set_idx_1])&(de_frac[set_idx_2]&!df_frac[set_idx_2]))|((df_frac[set_idx_2]&!de_frac[set_idx_2])&(de_frac[set_idx_1]&!df_frac[set_idx_1]))
        #(df_frac[set_idx_1]&de_frac[set_idx_2])&(df_frac[set_idx_2]&de_frac[set_idx_1])
      
      # if ((i < 5) & (j < 5)){
      #   print(paste(set_idx_1, set_idx_2, df_diff, de_diff))
      # }
      
      edge_list[edge_idx,1] <- df_diff
      edge_list[edge_idx,2] <- de_diff
      edge_list[edge_idx,3] <- euclidian_diff
      
      edge_idx <- edge_idx + 1
    }
  }
  
  return(edge_list)
}

g1_architecture_bootstrap <- function(coupling_mtx, df_frac, de_frac, n = 1000){
  sets <- rownames(coupling_mtx)
  
  boot_sets <- sample(sets, size = length(sets))
  rownames(coupling_mtx) <- boot_sets
  colnames(coupling_mtx) <- boot_sets
  
  edge_list <- g1_architecture_measurement_binary(coupling_mtx, df_frac, de_frac)
  
  for (i in 2:n){
    boot_sets <- sample(sets, size = length(sets))
    rownames(coupling_mtx) <- boot_sets
    colnames(coupling_mtx) <- boot_sets
    
    edge_list <- cbind(edge_list, g1_architecture_measurement_binary(coupling_mtx, df_frac, de_frac))
  }
  
  return(edge_list)
}

binary_edge_histogram <- function(edge_lists, og_edges, xlimits = c(0,1), title = '', color = 'lightskyblue'){
  
  bootstrapped <- c()
  for (i in 1:ncol(edge_lists)){
    bootstrapped <- c(bootstrapped, length(which(edge_lists[,i]))/nrow(edge_lists))
  }
  
  observed <- length(which(og_edges))/length(og_edges)
  # print(observed)
  bootstrapped <- data.frame(bootstrapped)
  colnames(bootstrapped) <- 'frac'
  ggplot(bootstrapped, aes(x=frac)) + ggtitle(title) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.01,
                   colour="black", fill=color) +
    geom_density() + scale_x_continuous(limits = c(0,1)) +
    geom_vline(xintercept = observed, colour = 'red', linetype="dashed") + xlab('Fraction') + ylab('Count') +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
  
  # qplot(bootstrapped, geom = 'histogram', bins = 50, xlim = xlimits, colour="black", fill=color) + 
  #   geom_vline(xintercept = observed, colour = 'red', linetype="dashed") + xlab('Fraction') + ylab('count')
}

node_degree <- function(coupling_mtx, interest_frac){
  set_idxs <- as.numeric(rownames(coupling_mtx))
  
  degree_frac <- matrix(data = 0, nrow = 2, ncol = length(set_idxs))
  
  for (i in 1:nrow(coupling_mtx)){
    degree <- length(which(coupling_mtx[,i])) + length(which(coupling_mtx[i,]))
    frac <- interest_frac[set_idxs[i]]
    if (is.na(frac)){
      frac <- -0.1
    }
    
    degree_frac[1,i] <- degree
    degree_frac[2,i] <- frac
  }
  
  return(degree_frac)
}
