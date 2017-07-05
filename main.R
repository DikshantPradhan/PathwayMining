redundant_id <- c()

for (i in 1:95){
  print(model@react_id[i])
  if (length(get_set_idx(model@react_id[i], og_set_list)) != 0){
    og_set <- og_set_list[[get_set_idx(model@react_id[i], og_set_list)]]
    if (length(og_set) > 1){
      redundant_id <- c(redundant_id, model@react_id[i])
      print(og_set)
      print(check_set_list_for_containing(og_set, set_lists = set_lists))
    }
  }
}