# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')

# for (i in 1:length(mutans_model@react_id)){
#   if (mutans_model@lowbnd[i] != mutans_falcon@lowbnd[i]){print(c(i, "wrong lowbnd"))}
#   if (mutans_model@uppbnd[i] != mutans_falcon@uppbnd[i]){print(c(i, "wrong uppbnd"))}
# }

yeast_falcon_compar_react <- function(react_id){
  idx <- which(yeast_model@react_id == react_id)
  f_idx <- which(yeast_falcon@react_id == react_id)
  
  met_idx <- which(yeast_model@S[,idx] != 0)
  f_met_idx <- which(yeast_falcon@S[,f_idx] != 0)
  
  met_coeff <- yeast_model@S[met_idx, idx]
  f_met_coeff <- yeast_falcon@S[f_met_idx, f_idx]
  
  met_id <- yeast_model@met_id[met_idx]
  f_met_id <- yeast_falcon@met_id[f_met_idx]
  
  print(react_id)
  print('Yeast Model')
  print(paste(met_coeff, met_id))
  print(paste(yeast_model@lowbnd[idx], yeast_model@uppbnd[idx]))
  print('Yeast Falcon')
  print(paste(f_met_coeff, f_met_id))
  print(paste(yeast_falcon@lowbnd[f_idx], yeast_falcon@uppbnd[f_idx]))
}

yeast_falcon_compar_met <- function(met_id){
  met_idx <- which(yeast_model@met_id == met_id)
  f_met_idx <- which(yeast_falcon@met_id == met_id)
  
  r_idxs <- which(yeast_model@S[met_idx,] != 0)
  f_r_idxs <- which(yeast_falcon@S[f_met_idx,] != 0)
  
  print(met_id)
  print('Yeast Model')
  print(yeast_model@react_id[r_idxs])
  print('Yeast Falcon')
  print(yeast_falcon@react_id[f_r_idxs])
}

# yeast_falcon_compar_met("s_0097")
# yeast_falcon_compar_met("s_0102")
# yeast_falcon_compar_react("r_0036")
# yeast_falcon_compar_react("r_0084")
# yeast_falcon_compar_react("r_0360")
# yeast_falcon_compar_react("r_1414")
#yeast_falcon_compar_react("EX_s_0199")

model@react_id[3484]
biomass_mets <- which(model@S[,3484] != 0)
for (met in biomass_mets){
  print(paste(model@met_id[met], length(which(model@S[met,] != 0))))
}
