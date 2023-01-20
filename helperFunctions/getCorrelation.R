getCorrelation = function(data,pheno1,pheno2,cor_method = "pearson"){
  # data = copy(result.1)
  # pheno1 = "PCSK9_male_treated"
  # pheno2 = "PCSK9_female_treated"
  
  # Filter data for phenotype
  data1 = copy(data)
  data1 = data1[phenotype == pheno1,]
  data2 = copy(data)
  data2 = data2[phenotype == pheno2,]
  
  # Filter data for valid SNPs
  data1 = data1[invalidAssocs  == F,]
  data2 = data2[invalidAssocs  == F,]
  
  # Filter data for same markernames
  data1 = data1[markername %in% data2$markername,]
  data2 = data2[markername %in% data1$markername,]
  
  # Test pearson correlation
  cor = cor.test(data1$beta, data2$beta,method = cor_method)
  res = data.table(pheno1 = pheno1,
                   pheno2 = pheno2,
                   correlation = cor$estimate,
                   pvalue = cor$p.value,
                   method = cor_method,
                   N_SNPs = dim(data1)[1])
  # Return result
  return(res)
}
