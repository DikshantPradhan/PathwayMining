# S <- model@S

g <- make_empty_graph()
g <- g + vertices(model@met_id, color = "red") + vertices(model@react_id, color = "green")

for (i in 1:ncol(S)){
  for (j in which(S[,i] < 0)){
    g <- g + edge(model@met_id[j], model@react_id[i], weight = avg_sample[i])
  }
  
  for (j in which(S[,i] > 0)){
    g <- g + edge(model@react_id[i], model@met_id[j], weight = avg_sample[i])
  }
}

#plot(g)
plot(g, rescale = FALSE, ylim=c(-15,15),xlim=c(-15,15), asp = 0)
