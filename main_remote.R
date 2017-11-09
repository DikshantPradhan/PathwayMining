
r0_size_class_2 <- c()

# randomize order of sets classified by size
for (i in 1:length(r0_size_class)){
  #print(i)
  #r0_size_class_2[[i]] <- sample(r0_size_class[[i]], size = length(r0_size_class[[i]]), replace = FALSE)
  #if (i == 13 | i == 21){
  #  print(r0_size_class[[i]])
  #  print(r0_size_class_2[[i]])
  #}
  #if (length(unique(r0_size_class_2[[i]])) != length(unique(r0_size_class[[i]]))){print('incorrect length')}
  if (length(r0_size_class[[i]]) == 1){
    r0_size_class_2[[i]] <- r0_size_class[[i]]
  }
  else {
    r0_size_class_2[[i]] <- sample(r0_size_class[[i]], size = length(r0_size_class[[i]]), replace = FALSE)
  }
}

total <- 0
for (i in 1:length(r0_size_class_2)){
  for (j in r0_size_class_2[[i]]){
    total <- total + length(r0_set_list[[j]])
  }
}

print(total)

for (i in 1:length(r0_size_class_2)){
  for (j in r0_size_class_2[[i]]){
    if (length(r0_set_list[[j]]) != i){
      print('wrong size')
      print(paste(i, ' ', j))
      print(r0_set_list[[j]])
    }
  }
}
