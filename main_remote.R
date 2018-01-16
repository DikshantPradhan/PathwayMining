# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')

yeast_model <- GRB_yeast_model()
yeast_falcon_model <- GRB_yeast_falcon_model()

yeast_set_list <- GRB_generate_set_list(yeast_model)
yeast_falcon_set_list <- GRB_generate_set_list(yeast_falcon_model)

# for (rxn_idx in 1:length(yeast_open_mod@react_id)){
#   if (!(yeast_open_mod@react_id[rxn_idx] %in% unlist(yeast_falcon_test))){next}
#   for (gene in yeast_open_mod@genes[[rxn_idx]]){
#     if (gene == ''){next}
#     if (!(gene %in% unlist(yeast_falcon_test))){
#       print(paste('error', rxn_idx, gene))
#     }
#   }
# }
# 
# for (rxn_idx in 1:length(mutans@react_id)){
#   if (!(mutans@react_id[rxn_idx] %in% unlist(mutans_falcon_test))){next}
#   for (gene in mutans@genes[[rxn_idx]]){
#     if (gene == ''){next}
#     if (!(gene %in% unlist(mutans_falcon_test))){
#       print(paste('error', rxn_idx, mutans@react_id[rxn_idx], gene))
#     }
#   }
# }

# > mutans_og_set_list <- GRB_generate_set_list(GRB_mutans_model())
# Academic license - for non-commercial use only
# Academic license - for non-commercial use only
# [1] "new rxn list"
# [1] "blocked: 22"
# [1] "blocked: 396"
# [1] "blocked: 489"
# [1] "blocked: 490"
# [1] "blocked: 493"
# [1] "blocked: 499"
# [1] "blocked: 529"
# [1] "blocked: 554"
# [1] "blocked: 555"
# [1] "blocked: 558"
# [1] "blocked: 563"
# [1] "blocked: 607"
# [1] "blocked: 651"
# [1] 1110
# There were 50 or more warnings (use warnings() to see the first 50)
# > mutans_falcon_og_set_list <- GRB_generate_set_list(GRB_mutans_falcon_model())
# Academic license - for non-commercial use only
# Academic license - for non-commercial use only
# [1] "new rxn list"
# [1] "blocked: 22"
# [1] "blocked: 396"
# [1] "blocked: 489"
# [1] "blocked: 490"
# [1] "blocked: 493"
# [1] "blocked: 499"
# [1] "blocked: 529"
# [1] "blocked: 554"
# [1] "blocked: 555"
# [1] "blocked: 558"
# [1] "blocked: 563"
# [1] "blocked: 607"
# [1] "blocked: 651"
# [1] "blocked: 726"
# [1] "blocked: 1001"
# [1] "blocked: 1079"
# [1] "blocked: 1080"
# [1] "blocked: 1081"
# [1] "blocked: 1088"
# [1] "blocked: 1089"
# [1] "blocked: 1090"
# [1] 5394

