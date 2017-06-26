library(igraph)

S <- model@S

# g <- make_empty_graph()
# g <- g + vertices(model@met_id, color = "red") + vertices(model@react_id, color = "green")
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

graph_rxn_sets <- function(set_list){
  g <- make_empty_graph()
  g <- g + vertices(unique(unlist(set_list)), color = "green")
  
  for (set in set_list){
    for (rxn in set){
      
      # ds_path <- get_dwnst_paths(get_rxn_idx(rxn))
      
      ds_rxns <- intersect(get_dwnst_paths(get_rxn_idx(rxn)), get_rxn_idx(set))
      print(ds_rxns)
      
      us_rxns <- intersect(get_upst_paths(get_rxn_idx(rxn)), get_rxn_idx(set))
      print(us_rxns)
      
      for (ds_rxn in ds_rxns){
        g <- g + edge(rxn, get_rxn_id_from_idx(ds_rxn))
      }
      
      for (us_rxn in us_rxns){
        g <- g + edge(get_rxn_id_from_idx(us_rxn), rxn)
      }
    }
    #print(set)
  }
  
  plot(g)
}
