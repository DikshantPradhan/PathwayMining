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

# > unlist(mutans_og_set_list)[mutans_missing]
# [1] "R00130"     "R03035"     "R00190"     "R00259"     "dem00007"  
# [6] "R00115"     "R00094"     "R00428"     "R03067"     "R02237"    
# [11] "R03503"     "R03504"     "R04620"     "R04639"     "R05046"    
# [16] "R05048"     "trans00008" "exc00057"   "dem00005"   "R00497"    
# [21] "R00894"     "R00513"     "R00516"     "R00517"     "R00549"    
# [26] "trans00013" "exc00012"   "R00650"     "R01291"     "R00194"    
# [31] "dem00015"   "R00703"     "dem00020"   "R00955"     "R01092"    
# [36] "R10619"     "R00964"     "R00967"     "R00968"     "R01057"    
# [41] "R01101"     "R01103"     "R01132"     "R01224"     "R07168"    
# [46] "R01229"     "R01547"     "R01548"     "R01549"     "R01560"    
# [51] "R01561"     "R01863"     "R01567"     "trans00028" "exc00037"  
# [56] "R01878"     "R01880"     "R01969"     "R01968"     "R02061"    
# [61] "dem00012"   "R02090"     "R02091"     "R02096"     "R02097"    
# [66] "R02142"     "R02297"     "trans00029" "exc00038"   "R02147"    
# [71] "trans00030" "exc00039"   "R02236"     "R02296"     "trans00027"
# [76] "exc00036"   "R02327"     "R02332"     "R02371"     "R02372"    
# [81] "R02410"     "R02556"     "R02748"     "R02557"     "R02568"    
# [86] "dem00014"   "R03269"     "R03634"     "R03707"     "R10505"    
# [91] "R09078"     "R02948"     "dem00010"   "R05134"     "R04394"    
# [96] "exc00024"   "dem00002"   "R05549"     "R03635"     "R06987"    
# [101] "R00921"     "R01353"     "trans00041" "exc00052"   "R04231"    
# [106] "R03018"     "trans00009" "exc00008"   "R04230"     "R01104"    
# [111] "trans00023" "exc00032"   "R01329"     "R01326"     "trans00025"
# [116] "exc00034"   "R02926"     "R02865"     "trans00026" "exc00035"  
# [121] "R01194"     "trans00024" "exc00033"   "dem00009"   "R11262"    
# [126] "trans00032" "exc00043"   "R11319"     "dem00003"   "R00937"    
# [131] "R00940"     "R02235"     "R02088"     "R10506"     "R10506v2"  
# [136] "R02408"     "trans00045" "exc00056"   "dem00021"   "R02971"    
# [141] "trans00010" "exc00009"   "R04391"     "trans00011" "exc00010"  
# [146] "trans00022" "exc00031"   "trans00033" "exc00044"   "trans00036"
# [151] "exc00047"   "trans00037" "exc00048"   "trans00038" "exc00049"  
# [156] "dem00004"   "dem00018"  

