TwoWayInteractionTest_jp = function(data,pheno1,pheno2,type,corCol){
  # data = copy(result.2)
  # pheno1 = "females"
  # pheno2 = "males"
  # type = "sexIA"
  # corCol = "cor"
  
  # Filter data by stratum
  data1 = copy(data)
  data1 = data1[grepl(paste0("_",pheno1),phenotype),]
  data2 = copy(data)
  data2 = data2[grepl(paste0("_",pheno2),phenotype),]

  # Filter data by best associated phenotype
  if(type == "sexIA"){
    data1 = data1[statin == bestStatin,]
    data2 = data2[statin == bestStatin,]
  }else{
    data1 = data1[sex == bestSex,]
    data2 = data2[sex == bestSex,]
  }

  # build information
  result.3 = copy(data1)
  myNames = c("markername","chr", "bp_hg19" ,"EA" ,   "OA",   "region" , "candidateGene", "bestPheno",corCol)
  colsOut = setdiff(colnames(result.3), myNames)
  result.3[, get("colsOut") := NULL]
  setnames(result.3,corCol,"cor")
  result.3[,type:= type]
  if(type == "sexIA"){
    result.3[,fix:= data1$statin]
  }else{
    result.3[,fix:= data1$sex]
  }

  result.3[,trait1 := data1$phenotype]
  result.3[,trait1_EAF := data1$EAF ]
  result.3[,trait1_beta := data1$beta ]
  result.3[,trait1_SE := data1$SE ]
  result.3[,trait1_pval := data1$pval ]
  result.3[,trait1_nSamples := data1$nSamples ]
  
  result.3[,trait2 := data2$phenotype]
  result.3[,trait2_EAF := data2$EAF ]
  result.3[,trait2_beta := data2$beta ]
  result.3[,trait2_SE := data2$SE ]
  result.3[,trait2_pval := data2$pval ]
  result.3[,trait2_nSamples := data2$nSamples ]
  
  result.3[,IA_diff := trait2_beta - trait1_beta]
  result.3[,IA_SE := sqrt(trait1_SE^2 + trait2_SE^2 - 2*cor*trait1_SE*trait2_SE) ]
  result.3[,IA_Zscore := IA_diff/IA_SE ]
  result.3[,IA_pval := pnorm(abs(IA_Zscore), lower.tail = F)*2 ]
  
  return(result.3)
}
