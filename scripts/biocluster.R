## script for biocluster

cluster_func <- function(){
  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes,
                                              compare_mtx = compare_known_r0_sets, known_set_mtx = r0_coupling_mtx)$coupled
  
  return(coupling_mtx)
}

lapply()