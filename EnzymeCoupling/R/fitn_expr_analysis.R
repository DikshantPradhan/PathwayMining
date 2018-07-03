# g1_sets <- pao_g1_sets
g0_sets <- g0_df$clean.sets #clean_rxn_names_in_set(g0_sets)
g1_sets <- g1_df$clean.sets #clean_rxn_names_in_set(g1_sets)

# set data to be used:
pao_data <- pao1#pa14[which(pa14$condition == unique(pa14$condition)[1]),]

# if using pa14 dataset, swap pa14 genes for pao1 homologs
load("~/GitHub/PathwayMining/EnzymeCoupling/data/pao1_key.Rdata")
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

g0_set_df <- gene_set_dataframe(g0_sets)
g1_set_df <- gene_set_dataframe(g1_sets)

binary_df_g0_sets <- binary_set_assignment(g0_sets, df_genes)
binary_df_g0_assign <- binary_set_histogram(binary_df_g0_sets)
binary_de_g0_sets <- binary_set_assignment(g0_sets, de_genes)
binary_de_g0_assign <- binary_set_histogram(binary_de_g0_sets)
g0_set_df <- cbind(g0_set_df, binary_df_g0_assign, binary_de_g0_assign)
colnames(g0_set_df) <- c('sets', 'clean sets', '# genes', 'df frac', 'de frac')
g0_df_obs <- g0_set_df[which(g0_set_df[,3] > 1),4]
g0_de_obs <- g0_set_df[which(g0_set_df[,3] > 1),5]

binary_df_g1_sets <- binary_set_assignment(g1_sets, df_genes)
binary_df_g1_assign <- binary_set_histogram(binary_df_g1_sets)
binary_de_g1_sets <- binary_set_assignment(g1_sets, de_genes)
binary_de_g1_assign <- binary_set_histogram(binary_de_g1_sets)
g1_set_df <- cbind(g1_set_df, binary_df_g1_assign, binary_de_g1_assign)
colnames(g1_set_df) <- c('sets', 'clean sets', '# genes', 'df frac', 'de frac')
g1_df_obs <- g1_set_df[which(g1_set_df[,3] > 1),4]
g1_de_obs <- g1_set_df[which(g1_set_df[,3] > 1),5]

num_g0_sets <- length(which(g0_set_df[,3] > 1))
num_g1_sets <- length(which(g1_set_df[,3] > 1))

n = 1000
g0_df_vec <- matrix(data = 0, nrow = num_g0_sets, ncol = 1)
g0_de_vec <- matrix(data = 0, nrow = num_g0_sets, ncol = 1)
g1_df_vec <- matrix(data = 0, nrow = num_g1_sets, ncol = 1)
g1_de_vec <- matrix(data = 0, nrow = num_g1_sets, ncol = 1)

print(num_de)
print(num_df)

for (i in 1:n){
  df_genes <- sample(genes_in_model, size = num_df, replace = FALSE)
  de_genes <- sample(genes_in_model, size = num_de, replace = FALSE)

  binary_df_g0_sets <- binary_set_assignment(g0_sets, df_genes)
  binary_df_g0_assign <- binary_set_histogram(binary_df_g0_sets)
  binary_de_g0_sets <- binary_set_assignment(g0_sets, de_genes)
  binary_de_g0_assign <- binary_set_histogram(binary_de_g0_sets)
  g0_set_df <- cbind(g0_set_df_og, binary_df_g0_assign, binary_de_g0_assign)
  colnames(g0_set_df) <- c('sets', 'clean sets', '# genes', 'df frac', 'de frac')
  # print(dim(g0_set_df))
  # sets_of_interest <- g0_set_df[which(g0_set_df[,3] > 1),4]
  g0_df_vec <- cbind(g0_df_vec, g0_set_df[which(g0_set_df[,3] > 1),4])
  # g0_df_sets <- c(g0_df_sets, sets_of_interest)
  # sets_of_interest <- matrix(data = sets_of_interest, nrow = 1, ncol = length(sets_of_interest))
  # hist(as.numeric(sets_of_interest), breaks=20, freq = FALSE)
  # sets_of_interest <- g0_set_df[which(g0_set_df[,3] > 1),5]
  g0_de_vec <- cbind(g0_de_vec, g0_set_df[which(g0_set_df[,3] > 1),5])
  # g0_de_sets <- c(g0_de_sets, sets_of_interest)
  # sets_of_interest <- matrix(data = sets_of_interest, nrow = 1, ncol = length(sets_of_interest))
  # hist(as.numeric(sets_of_interest), breaks=20, freq = FALSE)

  binary_df_g1_sets <- binary_set_assignment(g1_sets, df_genes)
  binary_df_g1_assign <- binary_set_histogram(binary_df_g1_sets)
  binary_de_g1_sets <- binary_set_assignment(g1_sets, de_genes)
  binary_de_g1_assign <- binary_set_histogram(binary_de_g1_sets)
  g1_set_df <- cbind(g1_set_df_og, binary_df_g1_assign, binary_de_g1_assign)
  colnames(g1_set_df) <- c('sets', 'clean sets', '# genes', 'df frac', 'de frac')
  # print(dim(g1_set_df))
  # sets_of_interest <- g1_set_df[which(g1_set_df[,3] > 1),4]
  # g1_df_sets <- c(g1_df_sets, sets_of_interest)
  g1_df_vec <- cbind(g1_df_vec, g1_set_df[which(g1_set_df[,3] > 1),4])
  # sets_of_interest <- matrix(data = sets_of_interest, nrow = 1, ncol = length(sets_of_interest))
  # hist(as.numeric(sets_of_interest), breaks=20, freq = FALSE)
  # sets_of_interest <- g1_set_df[which(g1_set_df[,3] > 1),5]
  # g1_de_sets <- c(g1_de_sets, sets_of_interest)
  g1_de_vec <- cbind(g1_de_vec, g1_set_df[which(g1_set_df[,3] > 1),5])
  # sets_of_interest <- matrix(data = sets_of_interest, nrow = 1, ncol = length(sets_of_interest))
  # hist(as.numeric(sets_of_interest), breaks=20, freq = FALSE)
}

plot_density(g0_df_vec, ylimits = c(0,90))
freqs <- table(as.numeric(g0_df_obs))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'blue')
title(main="Differential Fitness Composition in G0 Sets", xlab="% of Genes in Set with Differential Fitness", ylab="Frequency")
legend(0.45, 80, legend = c('bootstrapped', 'observed'), col = c('grey', 'blue'), lwd = 1)

plot_density(g0_de_vec, ylimits = c(0,60))
freqs <- table(as.numeric(g0_de_obs))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'chartreuse4')
title(main="Differential Expression Composition in G0 Sets", xlab="% of Genes in Set with Differential Expression", ylab="Frequency")
legend(0.45, 60, legend = c('bootstrapped', 'observed'), col = c('grey', 'chartreuse4'), lwd = 1)

plot_density(g1_df_vec, ylimits = c(0,80))
freqs <- table(as.numeric(g1_df_obs))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'blue')
title(main="Differential Fitness Composition in G1 Sets", xlab="% of Genes in Set with Differential Fitness", ylab="Frequency")
legend(0.45, 80, legend = c('bootstrapped', 'observed'), col = c('grey', 'blue'), lwd = 1)

plot_density(g1_de_vec, ylimits = c(0,60))
freqs <- table(as.numeric(g1_de_obs))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'chartreuse4')
title(main="Differential Expression Composition in G1 Sets", xlab="% of Genes in Set with Differential Expression", ylab="Frequency")
legend(0.45, 60, legend = c('bootstrapped', 'observed'), col = c('grey', 'chartreuse4'), lwd = 1)

pure_df_genes <- pao_data$gene[which(pao_data$sig_df & !pao_data$sig_de)]
pure_df_genes <- df_genes[which(pure_df_genes %in% genes_in_model)]
df_genes <- pao_data$gene[which(pao_data$sig_df)]
df_genes <- df_genes[which(df_genes %in% genes_in_model)]
pure_de_genes <- pao_data$gene[which(pao_data$sig_de & !pao_data$sig_df)]
pure_de_genes <- df_genes[which(pure_de_genes %in% genes_in_model)]
de_genes <- pao_data$gene[which(pao_data$sig_de)]
de_genes <- de_genes[which(de_genes %in% genes_in_model)]
dfde_genes <- pao_data$gene[which(pao_data$sig_de & pao_data$sig_df)]
dfde_genes <- dfde_genes[which(dfde_genes %in% genes_in_model)]
num_df <- length(df_genes)
num_de <- length(de_genes)
num_dfde <- length(dfde_genes)
interest_genes <- unique(c(df_genes, de_genes))

g0_df <- gene_set_expanded_dataframe(g0_set_df_og, g0_sets, df_genes, de_genes, interest_genes, pure_df_genes, pure_de_genes, dfde_genes)
g1_df <- gene_set_expanded_dataframe(g1_set_df_og, g1_sets, df_genes, de_genes, interest_genes, pure_df_genes, pure_de_genes, dfde_genes)

# g0_df <- data.frame(g0_set_df)
# g0_df$df <- binary_set_assignment(g0_sets, df_genes)
# g0_df$de <- binary_set_assignment(g0_sets, de_genes)
# g0_df$interest_genes <- binary_set_assignment(g0_sets, interest_genes)
# g0_df$pure_df <- binary_set_assignment(g0_sets, pure_df_genes)
# g0_df$pure_de <- binary_set_assignment(g0_sets, pure_de_genes)
# g0_df$dfde <- binary_set_assignment(g0_sets, dfde_genes)
# g0_df$df_frac <- binary_set_histogram(g0_df$df)
# g0_df$de_frac <- binary_set_histogram(g0_df$de)
# g0_df$pure_df_frac <- binary_set_histogram(g0_df$pure_df)
# g0_df$pure_de_frac <- binary_set_histogram(g0_df$pure_de)
# g0_df$dfde_frac <- binary_set_histogram(g0_df$dfde)
# g0_df$interest_frac <- binary_set_histogram(g0_df$interest_genes)

# set <- data.frame(as.numeric(g0_df_obs), as.numeric(g0_de_obs))
# 
# ggplot(g0_df, aes(x=df_frac, y=de_frac)) + geom_point(aes(size=dfde_frac), alpha = 0.1) + ggtitle('g0 Sets')
# ggplot(g1_df, aes(x=df_frac, y=de_frac)) + geom_point(aes(size=dfde_frac), alpha = 0.1) + ggtitle('g1 Sets')

print(df_genes)
print(de_genes)
print(dfde_genes)

plot_set_characteristic_hypothesis(g0_df_obs, g0_df_vec, color = 'lightskyblue') #Distribution of Differential Fitness Genes in G0 Sets
plot_set_characteristic_hypothesis(g0_de_obs, g0_de_vec, color = 'chartreuse3') #Distribution of Differential Expression Genes in G0 Sets
plot_set_characteristic_hypothesis(g1_df_obs, g1_df_vec, color = 'lightskyblue') #Distribution of Differential Fitness Genes in G1 Sets
plot_set_characteristic_hypothesis(g1_de_obs, g1_de_vec, color = 'chartreuse3') #Distribution of Differential Expression Genes in G1 Sets

# qplot(g1_df$dfde_frac[which(g1_df$X..genes > 1)], geom = 'histogram', xlab = 'composition', ylab = 'count', main = 'g1 dedf genes in sets', bins = 20)

# DFDE analysis: count

bootstrapped_dfde_count <- c()
num_genes <- length(genes_in_model)
for (i in 1:1000){
  dfde_binary <- matrix(data = 0, nrow = 1, ncol = num_genes)
  df <- sample(1:num_genes, size = num_df, replace = FALSE)
  de <- sample(1:num_genes, size = num_de, replace = FALSE)
  dfde_binary[df] <- dfde_binary[df] + 1
  dfde_binary[de] <- dfde_binary[de] + 1

  bootstrapped_dfde_count <- c(bootstrapped_dfde_count, length(which(dfde_binary == 2)))
}
bootstrapped_dfde_count <- data.frame(bootstrapped_dfde_count)
colnames(bootstrapped_dfde_count) <- 'ct'
ggplot(bootstrapped_dfde_count, aes(x=ct)) + ggtitle('Number of Genes with Differential Fitness and Expression') + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=2,
                 colour="black", fill='mediumpurple1') +
  geom_density() + scale_x_continuous(limits = c(0,125)) +
  xlab('Number of Genes') + ylab('Density') +
  geom_vline(xintercept = num_dfde, colour = 'red', linetype="dashed")
# qplot(bootstrapped_dfde_count, geom = 'histogram', bins = 50, main = 'Bootstrapped DFDE gene couont', xlim = c(0,150)) + geom_vline(xintercept = num_dfde, colour = 'red') + xlab('number of dfde genes') + ylab('count')

# DFDE analysis: sets

reduced_g0_df <- g0_df[which(g0_df$X..genes > 1),]
reduced_g1_df <- g1_df[which(g1_df$X..genes > 1),]

n_g0_sets <- nrow(reduced_g0_df)
n_g1_sets <- nrow(reduced_g1_df)


n = 1000
g0_dfde_vec <- matrix(data = 0, nrow = n_g0_sets, ncol = 1)
g1_dfde_vec <- matrix(data = 0, nrow = n_g1_sets, ncol = 1)

for (i in 1:n){
  g0_df_frac_sampled <- sample(reduced_g0_df$df_frac, size = n_g0_sets)
  g0_de_frac_sampled <- sample(reduced_g0_df$de_frac, size = n_g0_sets)
  # print(paste(length(g0_df_frac_sampled), length(g0_de_frac_sampled)))
  dfde_frac_list <- bootstrap_df_de_fracs(g0_df_frac_sampled, g0_de_frac_sampled, reduced_g0_df$X..genes)
  # print(dim(dfde_frac_list))
  g0_dfde_vec <- cbind(g0_dfde_vec, dfde_frac_list)

  g1_df_frac_sampled <- sample(reduced_g1_df$df_frac, size = n_g1_sets)
  g1_de_frac_sampled <- sample(reduced_g1_df$de_frac, size = n_g1_sets)
  dfde_sample <- bootstrap_df_de_fracs(g1_df_frac_sampled, g1_de_frac_sampled, reduced_g1_df$X..genes)
  g1_dfde_vec <- cbind(g1_dfde_vec, dfde_sample)
}
plot_density(g0_dfde_vec, ylimits = c(0,110))
freqs <- table(as.numeric(reduced_g0_df$dfde_frac))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="g0 dFdE set composition", xlab="% of genes with fitness defect & expression change", ylab="frequency")

plot_density(g1_dfde_vec, ylimits = c(0,110))
freqs <- table(as.numeric(reduced_g1_df$dfde_frac))
freqs_x <- as.numeric(names(freqs))
lines(freqs_x, freqs, col = 'red')
title(main="g1 dFdE set composition", xlab="% of genes with fitness defect & expression change", ylab="frequency")

# DFDE Resampling

# gene distributions in sets

interest_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$interest_genes[[x]]==1))}, c(1))
non_interest_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$interest_genes[[x]]==0))}, c(1))
df_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$pure_df[[x]]==1))}, c(1))
de_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$pure_de[[x]]==1))}, c(1))
reduced_g0_df <- reduced_g0_df[order(-interest_count,-df_count,de_count,non_interest_count),]
plot_genes_barplot(reduced_g0_df$sets, reduced_g0_df$pure_df, reduced_g0_df$pure_de, reduced_g0_df$dfde, reduced_g0_df$interest_genes)

interest_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$interest_genes[[x]]==1))}, c(1))
non_interest_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$interest_genes[[x]]==0))}, c(1))
df_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$pure_df[[x]]==1))}, c(1))
de_count <- vapply(1:nrow(reduced_g0_df), function(x){length(which(reduced_g0_df$pure_de[[x]]==1))}, c(1))
reduced_g0_df <- reduced_g0_df[order(-reduced_g0_df$interest_frac, -reduced_g0_df$pure_df_frac, reduced_g0_df$pure_de_frac),]
plot_genes_barplot(reduced_g0_df$sets, reduced_g0_df$pure_df, reduced_g0_df$pure_de, reduced_g0_df$dfde, reduced_g0_df$interest_genes, percent = TRUE)

interest_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$interest_genes[[x]]==1))}, c(1))
non_interest_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$interest_genes[[x]]==0))}, c(1))
df_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$pure_df[[x]]==1))}, c(1))
de_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$pure_de[[x]]==1))}, c(1))
reduced_g1_df <- reduced_g1_df[order(-interest_count,-df_count,de_count,non_interest_count),]
plot_genes_barplot(reduced_g1_df$sets, reduced_g1_df$pure_df, reduced_g1_df$pure_de, reduced_g1_df$dfde, reduced_g1_df$interest_genes)

interest_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$interest_genes[[x]]==1))}, c(1))
non_interest_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$interest_genes[[x]]==0))}, c(1))
df_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$pure_df[[x]]==1))}, c(1))
de_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$pure_de[[x]]==1))}, c(1))
reduced_g1_df <- reduced_g1_df[order(-reduced_g1_df$interest_frac, -reduced_g1_df$pure_df_frac, reduced_g1_df$pure_de_frac),]
plot_genes_barplot(reduced_g1_df$sets, reduced_g1_df$pure_df, reduced_g1_df$pure_de, reduced_g1_df$dfde, reduced_g1_df$interest_genes, percent = TRUE)

# analyze df de groupings in g1 sets 
# og_g1_sets <- g1_sets
cats <- c('full_df','full_dfde','full_de','full','partial','none')
# g1_sets <- og_g1_sets
g0_categories <- matrix(data = 'na', nrow = length(g0_sets), ncol = 1)
g1_set_idxs <- matrix(data = 0, nrow = length(g0_sets), ncol = 1)
value <- matrix(data = 0, nrow = length(g0_sets), ncol = 1)
for (i in 1:length(g0_sets)){
  cat <- 'none'
  if (g0_df$X..genes[i] < 1){next}
  if (g0_df$interest_frac[i] != 0){
    cat <- 'partial'
  }
  if (g0_df$interest_frac[i] == 1){
    cat <- 'full'
  }
  if (g0_df$df_frac[i] == 1){
    cat <- 'full_df'
  }
  if (g0_df$de_frac[i] == 1){
    cat <- 'full_de'
  }
  if (g0_df$dfde_frac[i] == 1){
    cat <- 'full_dfde'
  }
  
  g0_categories[i] <- cat
  
  g1_idx <- get_containing_set_idx(g0_sets[[i]], g1_sets)
  g1_set_idxs[i] <- g1_idx
  value[i] <- length(g0_sets[[i]])
}
remove <- which(value == 0 | g1_set_idxs == 0)
g0_categories <- g0_categories[-remove]
g1_set_idxs <- g1_set_idxs[-remove]
value <- value[-remove]
g1_composition_df <- data.frame(g1_set = g1_set_idxs, category = g0_categories, value = value)

comp_count <- vapply(1:length(g1_sets), function(x){length(which(g1_composition_df$g1_set == x))}, c(1))
g1_sets <- g1_sets[order(-comp_count)]
# g1_sets <- g1_sets[1:40]
# non_interest_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$interest_genes[[x]]==0))}, c(1))
# df_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$pure_df[[x]]==1))}, c(1))
# de_count <- vapply(1:nrow(reduced_g1_df), function(x){length(which(reduced_g1_df$pure_de[[x]]==1))}, c(1))
# reduced_g1_df <- reduced_g1_df[order(-reduced_g1_df$interest_frac, -reduced_g1_df$pure_df_frac, reduced_g1_df$pure_de_frac),]
g1_composition_df$category <- factor(g1_composition_df$category, levels = rev(cats), ordered = TRUE)
# g1_composition_df <- g1_composition_df[order(comp_count),]

g1_composition_df$value <- 1

ggplot() + geom_bar(data=g1_composition_df, aes(fill=category, y=value, x=g1_set), stat="identity") + labs(y = '# of G0 sets') +
  coord_flip() + ggtitle('G1 Set Composition')

binary_g0_categories <- matrix(data = 0, nrow = length(g0_categories), ncol = 1)
cats <- c('full_df','full_dfde','full_de','partial','none')
binary_g0_categories[grep('full', g0_categories)] <- 1
binary_g0_categories[grep(cats[3], g0_categories)] <- 1

g1_composition_sampling(binary_g0_categories, comp_count)

g1_ratio_sampling(g0_categories, comp_count, 1000, g1_set_idxs)
