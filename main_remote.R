for (rxn_idx in 1:length(yeast_open_mod@react_id)){
  if (!(yeast_open_mod@react_id[rxn_idx] %in% unlist(yeast_falcon_test))){next}
  for (gene in yeast_open_mod@genes[[rxn_idx]]){
    if (gene == ''){next}
    if (!(gene %in% unlist(yeast_falcon_test))){
      print(paste('error', rxn_idx, gene))
    }
  }
}

for (rxn_idx in 1:length(mutans@react_id)){
  if (!(mutans@react_id[rxn_idx] %in% unlist(mutans_falcon_test))){next}
  for (gene in mutans@genes[[rxn_idx]]){
    if (gene == ''){next}
    if (!(gene %in% unlist(mutans_falcon_test))){
      print(paste('error', rxn_idx, mutans@react_id[rxn_idx], gene))
    }
  }
}
