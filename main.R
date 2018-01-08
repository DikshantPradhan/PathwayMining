# source('~/GitHub/PathwayMining/network_tools.R')
source('~/GitHub/PathwayMining/sampling_tools.R')
source('~/GitHub/PathwayMining/fit_comparison_run.R')
source('~/GitHub/PathwayMining/model_tools.R')
source('~/GitHub/PathwayMining/set_tools.R')
source('~/GitHub/PathwayMining/falcon_tools.R')
source('~/GitHub/PathwayMining/load_mod.R')

for (i in 1:length(mutans_model@react_id)){
  if (mutans_model@lowbnd[i] != mutans_falcon@lowbnd[i]){print(c(i, "wrong lowbnd"))}
  if (mutans_model@uppbnd[i] != mutans_falcon@uppbnd[i]){print(c(i, "wrong uppbnd"))}
}