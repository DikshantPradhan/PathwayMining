library(igraph)
S <- model@S

g <- make_empty_graph()
# g <- g + vertices(model@met_id, color = "red")
g <- g + vertices(model@react_id, color = "green")


# 
# for (i in 1:ncol(S)){
#   for (j in which(S[,i] < 0)){
#     g <- g + edge(model@met_id[j], model@react_id[i], weight = avg_sample[i])
#   }
#   
#   for (j in which(S[,i] > 0)){
#     g <- g + edge(model@react_id[i], model@met_id[j], weight = avg_sample[i])
#   }
# }
# 
# #plot(g)
# plot(g, rescale = FALSE, ylim=c(-15,15),xlim=c(-15,15), asp = 0)

graph_rxn_sets <- function(set_list, edge_color = "grey"){
  graph <- make_empty_graph()
  graph <- graph + vertices(unique(unlist(set_list)), color = "green")
  
  graph <- rxn_set_edges(graph, list(unlist(rxn_set)), edge_color)
  graph <- rxn_set_edges(graph, set_list, edge_color = "blue")
  
  return(graph)
}

rxn_set_edges <- function(graph, set_list, edge_color){
  
  for (set in set_list){
    for (rxn in set){
      ds_rxns <- intersect(get_dwnst_paths(get_rxn_idx(rxn)), get_rxn_idx(set))
      us_rxns <- intersect(get_upst_paths(get_rxn_idx(rxn)), get_rxn_idx(set))
      
      for (ds_rxn in ds_rxns){
        graph <- graph + edge(rxn, get_rxn_id_from_idx(ds_rxn), color = edge_color)
      }
      for (us_rxn in us_rxns){
        graph <- graph + edge(get_rxn_id_from_idx(us_rxn), rxn, color = edge_color)
      }
    }
  }
  
  return(graph)
}

# graph_rxn_sets(list(unlist(rxn_set)))
rxn_sets <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$gained_pairs[5]))

# g <- graph_rxn_sets(list(unlist(rxn_set)), edge_color = "blue")
g <- graph_rxn_sets(rxn_set)

plot(g)