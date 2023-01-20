ThreeWayInteractionTest_jp = function(data,pheno1,pheno2,pheno3,pheno4){
  # data = copy(result.1_PCSK9)
  # pheno1 = "PCSK9_male_treated"
  # pheno2 = "PCSK9_female_treated"
  # pheno3 = "PCSK9_male_free"
  # pheno4 = "PCSK9_female_free"
  
  # Filter data
  data1 = copy(data)
  data1 = data1[phenotype == pheno1,]
  data2 = copy(data)
  data2 = data2[phenotype == pheno2,]
  data3 = copy(data)
  data3 = data3[phenotype == pheno3,]
  data4 = copy(data)
  data4 = data4[phenotype == pheno4,]
  
  # build information
  result.3 = copy(data1)[,c(1:5)]
  
  result.3[,trait1 := pheno1]
  result.3[,trait1_EAF := data1$EAF ]
  result.3[,trait1_beta := data1$beta ]
  result.3[,trait1_SE := data1$SE ]
  result.3[,trait1_pval := data1$pval ]
  result.3[,trait1_nSamples := data1$nSamples ]
  
  result.3[,trait2 := pheno2]
  result.3[,trait2_EAF := data2$EAF ]
  result.3[,trait2_beta := data2$beta ]
  result.3[,trait2_SE := data2$SE ]
  result.3[,trait2_pval := data2$pval ]
  result.3[,trait2_nSamples := data2$nSamples ]
  
  result.3[,trait3 := pheno3]
  result.3[,trait3_EAF := data3$EAF ]
  result.3[,trait3_beta := data3$beta ]
  result.3[,trait3_SE := data3$SE ]
  result.3[,trait3_pval := data3$pval ]
  result.3[,trait3_nSamples := data3$nSamples ]
  
  result.3[,trait4 := pheno4]
  result.3[,trait4_EAF := data4$EAF ]
  result.3[,trait4_beta := data4$beta ]
  result.3[,trait4_SE := data4$SE ]
  result.3[,trait4_pval := data4$pval ]
  result.3[,trait4_nSamples := data4$nSamples ]
  
  # two way:
  # diff = trait2 - trait1 = females - males = treated - free
  # SE = sqrt(SE1^2 + SE2^2 - r*SE1*SE2)
  
  # three way
  # diff = (trait2 - trait1) - (trait4 - trait3) = (females_treated - males_treated) - (females_free - males_free)
  # SE = sqrt(SE1^2 + SE2^2 + SE3^2 + SE4^2)
  
  result.3[,IA_diff := (trait2_beta - trait1_beta) - (trait4_beta - trait3_beta)]
  result.3[,IA_SE := sqrt(trait1_SE^2 + trait2_SE^2 + trait3_SE^2 + trait4_SE^2) ]
  result.3[,IA_Zscore := IA_diff/IA_SE ]
  result.3[,IA_pval := pnorm(abs(IA_Zscore), lower.tail = F)*2 ]
  
  return(result.3)
}
