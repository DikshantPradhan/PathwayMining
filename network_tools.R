library(igraph)
S <- model@S

og_pairs <- return_couples(flux_coupling_cor(sample_og))
og_rxn_set <- get_list_of_sets(og_pairs)

# g <- make_empty_graph()
# g <- g + vertices(model@met_id, color = "red")
# g <- g + vertices(model@react_id, color = "green")


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

add_edge <- function(graph, v1, v2, color){
  if (!are_adjacent(graph, v1, v2)){
    graph <- graph + edge(v1, v2, color = color)
  }
  return(graph)
}

add_vertex <- function(graph, v, color = "green"){
  if (!(v %in% names(V(graph)))){
    graph <- graph + vertices(v, color = color)
  }
  return(graph)
}

add_set_vertices <- function(graph, set_list, color = "green"){
  for (vertex in unlist(set_list)){
    graph <- add_vertex(graph, vertex)
  }
  
  return(graph)
}

graph_rxn_sets <- function(set_list, edge_color = "grey", show_theorietical = FALSE, direct = TRUE){
  graph <- make_empty_graph()
  graph <- graph + vertices(unique(unlist(set_list)), color = "green")
  
  graph <- rxn_set_edges(graph, set_list, edge_color, direct = direct)
  
  if (show_theorietical){
    graph <- rxn_set_edges(graph, list(unlist(set_list)), edge_color = "grey")
  }
  
  return(graph)
}

rxn_set_edges <- function(graph, set_list, edge_color, direct = TRUE){
  for (vertex in unlist(set_list)){
    graph <- add_vertex(graph, vertex)
  }
  
  for (set in set_list){
    for (rxn in set){
      ds_rxns <- c()
      us_rxns <- c()
      
      if (direct){
        ds_rxns <- intersect(get_rxn_id_from_idx(get_dwnst_rxns(get_rxn_idx(rxn))), set)
        us_rxns <- intersect(get_rxn_id_from_idx(get_upst_rxns(get_rxn_idx(rxn))), set)
      }
      else {
        ds_rxns <- intersect(get_rxn_id_from_idx(get_dwnst_paths(get_rxn_idx(rxn))), set)
        us_rxns <- intersect(get_rxn_id_from_idx(get_upst_paths(get_rxn_idx(rxn))), set)
      }
      
      for (ds_rxn in ds_rxns){
        # graph <- graph + edge(rxn, get_rxn_id_from_idx(ds_rxn), color = edge_color)
        graph <- add_edge(graph, rxn, ds_rxn, color = edge_color)
      }
      for (us_rxn in us_rxns){
        # graph <- graph + edge(get_rxn_id_from_idx(us_rxn), rxn, color = edge_color)
        graph <- add_edge(graph, us_rxn, rxn, color = edge_color)
      }
    }
  }
  
  return(graph)
}

graph_correlation_set <- function(set_list){
  g <- graph_rxn_sets(set_list, edge_color = "blue", show_theorietical = FALSE)
  plot(g)
}

graph_containing_set <- function(rxn_id, set_list, g = make_empty_graph()){
  set_idx <- grep(rxn_id, set_list)
  
  # g <- make_empty_graph()
  # g <- g + vertices(unique(set_list[[set_idx]]), color = "green")
  g <- add_set_vertices(g, unique(set_list[[set_idx]]))
  g <- rxn_set_edges(g, set_list[set_idx], edge_color = "blue")
  
  plot(g)
  
  return(g)
}

compare_containing_sets <- function(rxn_id, og_rxn_set, suppr_rxn_set, direct = TRUE){
  og_set_idx <- grep(rxn_id, og_rxn_set)
  suppr_set_idx <- grep(rxn_id, suppr_rxn_set)
  
  # print(og_set_idx)
  # print(suppr_set_idx)
  
  rxns <- union(og_rxn_set[[og_set_idx]], suppr_rxn_set[[suppr_set_idx]])
  
  # print(suppr_rxn_set[[suppr_set_idx]])
  
  g <- make_empty_graph()
  g <- g + vertices(unique(rxns), color = "green")
  
  g <- rxn_set_edges(g, suppr_rxn_set[suppr_set_idx], edge_color = "blue", direct = direct)
  g <- rxn_set_edges(g, og_rxn_set[og_set_idx], edge_color = "red", direct = direct)
  
  plot(g)
}

compare_model_sets <- function(comparison_num){
  og_pairs <- return_couples(flux_coupling_cor(sample_og))
  og_rxn_set <- get_list_of_sets(og_pairs)
  
  g <- make_empty_graph()
  g <- g + vertices(model@react_id, color = "green")
  g <- rxn_set_edges(g, og_rxn_set, "grey", direct = TRUE)
  
  compar_sample <- sampler(suppressed_model(model, comparison_num)) #sample_suppr[comparison_num, , ]
  suppr_pairs <- return_couples(flux_coupling_cor(compar_sample))
  suppr_rxn_set <- get_list_of_sets(suppr_pairs)
  
  #g <- rxn_set_edges(g, suppr_rxn_set, "blue")
  
  # added_rxns <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$gained_pairs[comparison_num]))
  # lost_rxns <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$lost_pairs[comparison_num]))
  # 
  # g <- rxn_set_edges(g, added_rxns, "green", direct = TRUE)
  # g <- rxn_set_edges(g, lost_rxns, "red", direct = TRUE)
  plot(g)
  
  return(suppr_rxn_set)
}

# graph_rxn_sets(list(unlist(rxn_set)))
# rxn_sets <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$gained_pairs[5]))
# 
# g <- graph_rxn_sets(list(unlist(rxn_set)), edge_color = "blue")
# g <- graph_rxn_sets(rxn_set)
# 
# plot(g)