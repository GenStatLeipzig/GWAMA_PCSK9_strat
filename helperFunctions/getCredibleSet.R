getCredibleSet<-function(filename, type, CHR, region_start, region_end, NR, SNP, pheno, candidateGene){
  # filename = IndepSignals$input_CS[10]
  # type = "GWAMA"
  # CHR = IndepSignals$Chr[10]
  # region_start = IndepSignals$region_start[10]
  # region_end = IndepSignals$region_end[10]
  # NR = IndepSignals$num[10]
  # SNP = IndepSignals$rsID[10]
  # pheno = IndepSignals$pheno[10]
  # candidateGene = IndepSignals$candidateGene[10]
  
  # Step 1: load data
  data = fread(filename)
  
  # Step 2: Harmonize column names (SNP, chr, pos, beta, se)
  if(type == "cond"){
    myNames1 = c("SNP", "Chr","bp","bC","bC_se","pC")
    data = data[!is.na(bC),]
  }else{
    data<-data[invalidAssocs==F,]
    myNames1 = c("markername","chr", "bp_hg19","beta","SE","pval")
  }
  stopifnot(sum(is.element(colnames(data),myNames1))==6)
  colsOut<-setdiff(colnames(data),myNames1)
  data[,get("colsOut"):=NULL]
  setcolorder(data,myNames1)
  myNames2<-c("SNP","chr","pos","beta","se","pval")
  names(data)<-myNames2
  
  # Step 3: Filter for region
  data = data[chr==CHR & (pos<=region_end & pos>=region_start),]
  
  # Step 4: Get prios
  quant_0.975 = quantile(x = data$beta, probs = 0.975, na.rm = T)
  quant_0.025 = quantile(x = data$beta, probs = 0.025, na.rm = T)
  range_quant = quant_0.975-quant_0.025
  prior = range_quant/(2*1.96)
  message("Using as prior ",prior)
  
  # Step 5: Calculate Bayes factors, Posterior Probability, and cummulative sum
  data[,abf.Wakefield := abf.Wakefield(beta = beta, se = se, priorsd = prior)]
  message("Summary of approximate Bayes factors: ")
  print(summary(data[,abf.Wakefield]))
  sum.abf.r1<-sum(data[,abf.Wakefield], na.rm=T)
  sum.abf.r1
  data[,PostProb:=abf.Wakefield/sum.abf.r1]
  message("Summary of Posterior Probabilities: ")
  print(summary(data[,PostProb]))
  ordering = order(data[,PostProb], decreasing = TRUE)
  data = data[ordering,]
  data[,SumProb:=cumsum(PostProb)]
  data
  
  # Step 6: Create summary
  message(paste0("Phenotype ",pheno, " at locus ",NR," with independent signal ",SNP, " (",type,", candidate gene: ",candidateGene,
                 ")",
                 "\n       Number of SNPs: ", nrow(data), 
                 "\n    Credible Set (99): ", min(which(data[,SumProb]>0.99)), 
                 "\n    Credible Set (95): ", min(which(data[,SumProb]>0.95))))
  
  # Step 7: Create plot
  x99 = min(which(data[,SumProb]>0.99))
  x95 = min(which(data[,SumProb]>0.95))
  
  if(x99<100){
    x_max = 100
    subtext = "only the first 100 SNPs plotted"
  }else{
    x_max = dim(data)[1]
    subtext = "all SNPs of region plotted"
  }
  
  plot(c(1:x_max), data[1:x_max,SumProb], type="l", ylim=c(0,1),
       xlab="position in CS", ylab="sum of Posterior Probability", 
       main=paste0("Phenotype ",pheno, " at locus ",NR," \nwith independent signal ",SNP, 
                   " (",type,", candidate gene: ",candidateGene,")"),
       sub=subtext)
  abline(v=x95,col="blue")
  abline(v=x99,col="red")
  abline(h=0.95,col="blue")
  abline(h=0.99,col="red")
  
  # Step 8: Save complete data set and return filtered data 
  if(dir.exists("../results/04_CredSets/")==F) dir.create("../results/04_CredSets/") 
  outfn = paste0("../results/04_CredSets/",NR,"::",candidateGene,"::",pheno,"::",SNP,".txt")
  fwrite(data,file=outfn,sep = "\t")
  
  data<-data[1:min(which(data[,SumProb]>0.99))]
  return(data)
}
