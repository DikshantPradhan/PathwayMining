library(igraph)

# og_pairs <- return_couples(flux_coupling_cor(sample_og))
# og_rxn_set <- get_list_of_sets(og_pairs)

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
    graph <- add_vertex(graph, vertex, color = color)
  }
  
  return(graph)
}

plot_graph <- function(graph){
  plot(graph,edge.arrow.size=0.3,vertex.label.color = "black",vertex.size=10)
}

graph_rxn_sets <- function(set_list, edge_color = "grey", show_theorietical = FALSE, direct = TRUE, sample = NULL){
  graph <- make_empty_graph()
  # graph <- graph + vertices(unique(unlist(set_list)), color = "green")
  
  graph <- rxn_set_edges(set_list, graph, edge_color, direct = direct)
  
  if (show_theorietical){
    graph <- rxn_set_edges(list(unlist(set_list)), graph, edge_color = "grey", sample = sample)
  }
  plot_graph(graph)
  return(graph)
}

rxn_set_edges <- function(set_list, graph = make_empty_graph(), edge_color = "blue", vertex_color = "green", 
                          direct = TRUE, addition = "", sample = NULL){
  for (vertex in unlist(set_list)){
    graph <- add_vertex(graph, paste(vertex, addition, sep = ""), color = vertex_color)
  }
  
  for (set in set_list){
    for (rxn in set){
      ds_rxns <- c()
      us_rxns <- c()
      
      rxn_idx <- get_rxn_idx(rxn)
      
      if (direct){
        ds_rxns <- intersect(get_rxn_id_from_idx(get_dwnst_rxns(rxn_idx, sample)), set)
        # us_rxns <- intersect(get_rxn_id_from_idx(get_upst_rxns(get_rxn_idx(rxn))), set)
      }
      else {
        ds_rxns <- intersect(get_rxn_id_from_idx(get_dwnst_paths(rxn_idx, sample)), set)
        # us_rxns <- intersect(get_rxn_id_from_idx(get_upst_paths(get_rxn_idx(rxn))), set)
      }
      
      for (ds_rxn in ds_rxns){
        # graph <- graph + edge(rxn, get_rxn_id_from_idx(ds_rxn), color = edge_color)
        graph <- add_edge(graph, paste(rxn, addition, sep = ""), paste(ds_rxn, addition, sep = ""), color = edge_color)
      }
      # for (us_rxn in us_rxns){
      #   # graph <- graph + edge(get_rxn_id_from_idx(us_rxn), rxn, color = edge_color)
      #   graph <- add_edge(graph, paste(us_rxn, addition, sep = ""), paste(rxn, addition, sep = ""), color = edge_color)
      # }
    }
  }
  
  return(graph)
}

graph_correlation_set <- function(set_list){
  g <- graph_rxn_sets(set_list, edge_color = "blue", show_theorietical = FALSE)
  plot_graph(g)
}

graph_containing_set <- function(rxn_id, set_list, g = make_empty_graph(), plot = TRUE, sample = NULL){
  set_idx <- get_set_idx(rxn_id, set_list) #grep(rxn_id, set_list)
  
  # g <- make_empty_graph()
  # g <- g + vertices(unique(set_list[[set_idx]]), color = "green")
  # print(paste(set_idx, ", ", rxn_id, ", "))
  
  if (length(set_idx) == 0){
    print(paste(rxn_id, " not found in set"))
    return(g)
  }
  
  g <- add_set_vertices(g, unique(set_list[[set_idx]]))
  g <- rxn_set_edges(set_list[set_idx], g, edge_color = "blue", sample = sample)
  
  if (plot){
    plot_graph(g)
  }
  return(g)
}

graph_multiple_sets <- function(rxn_ids, set_lists, graph = make_empty_graph(), plot = TRUE){
  
  for (i in 1:length(rxn_ids)){
    # print(rxn_ids[i])
    for (j in 1:length(set_lists)){
      graph <- graph_containing_set(rxn_ids[i], set_lists[[j]], graph, plot = FALSE)
    }
  }
  
  if (plot){
    plot_graph(graph)
  }
  return(graph)
}

compare_containing_sets <- function(rxn_id, og_rxn_set, suppr_rxn_set, direct = TRUE){
  og_set_idx <- get_set_idx(rxn_id, og_rxn_set) #grep(rxn_id, og_rxn_set)
  suppr_set_idx <- get_set_idx(rxn_id, suppr_rxn_set) #grep(rxn_id, suppr_rxn_set)
  
  rxns <- union(og_rxn_set[[og_set_idx]], suppr_rxn_set[[suppr_set_idx]])
  
  g <- make_empty_graph()
  g <- g + vertices(unique(rxns), color = "green")
  
  g <- rxn_set_edges(suppr_rxn_set[suppr_set_idx], g, edge_color = "blue", direct = direct)
  g <- rxn_set_edges(og_rxn_set[og_set_idx], g, edge_color = "red", direct = direct)
  
  plot_graph(g)
}

compare_model_sets <- function(comparison_num){
  og_pairs <- return_couples(flux_coupling_cor(sample_og))
  og_rxn_set <- get_list_of_sets(og_pairs)
  
  g <- make_empty_graph()
  g <- g + vertices(model@react_id, color = "green")
  g <- rxn_set_edges(og_rxn_set, g, "grey", direct = TRUE)
  
  compar_sample <- sampler(suppressed_model(model, comparison_num)) #sample_suppr[comparison_num, , ]
  suppr_pairs <- return_couples(flux_coupling_cor(compar_sample))
  suppr_rxn_set <- get_list_of_sets(suppr_pairs)
  
  #g <- rxn_set_edges(g, suppr_rxn_set, "blue")
  
  # added_rxns <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$gained_pairs[comparison_num]))
  # lost_rxns <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$lost_pairs[comparison_num]))
  # 
  # g <- rxn_set_edges(g, added_rxns, "green", direct = TRUE)
  # g <- rxn_set_edges(g, lost_rxns, "red", direct = TRUE)
  plot_graph(g)
  
  return(suppr_rxn_set)
}

compare_degen_sets_to_og <- function(rxn_id, og_rxn_set, suppr_rxn_set, addition = "", graph = make_empty_graph(), plot = TRUE){
  
  og_set_idx <- get_set_idx(rxn_id, og_rxn_set)
  suppr_set_idx <- get_set_idx(rxn_id, suppr_rxn_set)
  
  graph <- rxn_set_edges(og_rxn_set[og_set_idx], graph, edge_color = "black")
  # graph <- rxn_set_edges(graph, og_rxn_set[og_set_idx], edge_color = "grey", direct = FALSE)
  
  # print(suppr_set_idx)
  
  if (length(og_set_idx) > 0 & length(suppr_set_idx) > 0){
    if (!setequal(og_rxn_set[[og_set_idx]], suppr_rxn_set[[suppr_set_idx]])){
      # print("unequal set: ")
      # print(og_rxn_set[[og_set_idx]])
      # print(suppr_rxn_set[[suppr_set_idx]])
      graph <- rxn_set_edges(suppr_rxn_set[suppr_set_idx], graph, edge_color = "blue", vertex_color = "white", addition = addition)
      # graph <- rxn_set_edges(graph, suppr_rxn_set[suppr_set_idx], edge_color = "grey", vertex_color = "white",
      #                        addition = addition, direct = FALSE)
      graph <- add_edge(graph, rxn_id, paste(rxn_id, addition, sep = ""), color = "red")
    }
  }
  
  if (plot){
    plot_graph(graph)
  }
  return(graph)
}

compare_multiple_degen_sets <- function(rxn_ids, og_rxn_set, suppr_rxn_sets, graph = make_empty_graph(), plot = TRUE){
  
  # used_sets <- c()
  
  for (i in 1:length(rxn_ids)){
    # print(rxn_ids[i])
    for (j in 1:length(suppr_rxn_sets)){ # addition = paste("_", j, sep = ""),
      graph <- compare_degen_sets_to_og(rxn_ids[i], og_rxn_set, suppr_rxn_sets[[j]], addition = "_", graph = graph, plot = FALSE) 
        #graph_containing_set(rxn_ids[i], set_lists[[j]], graph, plot = FALSE)
      # used_sets <- c(used_sets, list())
    }
  }
  
  if (plot){
    plot_graph(graph)
  }
  return(graph)
}

graph_redundancies <- function(rxns){
  graph <- make_empty_graph()
  for (rxn in rxns){
    graph <- add_vertex(graph = graph, rxn, color = "yellow")
  }
  graph <- rxn_set_edges(set_list = list(rxns), graph = graph, sample = NULL)
  graph <- rxn_set_edges(set_list = list(total_union(find_all_sets_for_rxns(rxns, set_lists))), graph = graph, sample = NULL)
  plot_graph(graph)
  return(graph)
}

# graph_rxn_sets(list(unlist(rxn_set)))
# rxn_sets <- get_list_of_sets(convert_pair_strings_to_vector(media_cond$gained_pairs[5]))
# 
# g <- graph_rxn_sets(list(unlist(rxn_set)), edge_color = "blue")
# g <- graph_rxn_sets(rxn_set)
# 
# plot(g)