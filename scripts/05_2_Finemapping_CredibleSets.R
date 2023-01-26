#' ---
#' title: "Credible Sets"
#' subtitle: "PCSK9 GWAMA sex-stratified"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' I want to define the 95% and 99% credible sets for each indendent loci
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/abf.Wakefield.R")
source("../helperFunctions/getCredibleSet.R")

#' # Load result of GCTA ####
#' ***
load("../results/05_GCTA_COJO.RData")
IndepSignals

#' # Get Credible Sets per locus and signal ####
#' ***
dumTab = foreach(i=1:dim(IndepSignals)[1])%do%{
  #i=1
  myRow = IndepSignals[i,]
  if(myRow$multipleSignals==T){
    myType = "cond"
  }else{
    myType = "GWAMA"
  }
  
  credSet = getCredibleSet(filename = myRow$input_CS,
                           type = myType,
                           CHR = myRow$Chr,
                           region_start = myRow$region_start,
                           region_end = myRow$region_end,
                           NR = myRow$num,
                           SNP = myRow$rsID,
                           pheno = myRow$pheno,
                           candidateGene = myRow$candidateGene)
  
  credSet[,phenotype := myRow$pheno]
  credSet[,locus := myRow$num]
  credSet[,leadSNP := myRow$rsID]
  credSet[,candidateGene := myRow$candidateGene]
  credSet[,type := myType]
  credSet[,rank := 1:dim(credSet)[1]]
  credSet[,flag := F]
  credSet[1:min(which(credSet[,SumProb]>0.95)),flag:=T]
  credSet
}

credSet = rbindlist(dumTab)

#' # Some checks and plots ####
#' ***
credSet[,group := paste(locus,leadSNP,sep="::")]
credSet[,table(group)]
credSet

myPlot1 = ggplot(credSet[candidateGene=="PCSK9",], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/05_CredSet_PCSK9_signals.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot1
dev.off()

myPlot2 = ggplot(credSet[candidateGene!="PCSK9",], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/05_CredSet_other_signals.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot2
dev.off()

#' # Save ####
#' ***
names(credSet)
setnames(credSet,"SNP","snp")
save(credSet, file="../results/05_CredSet.RData")

credSet99 = copy(credSet)
credSet99 = credSet99[candidateGene == "PCSK9",]
credSet99 = credSet99[,c(1,17,8,9)] 

credSet95 = copy(credSet)
credSet95 = credSet95[candidateGene == "PCSK9",]
credSet95 = credSet95[flag == T,] 
credSet95 = credSet95[,c(1,17,8,9)] 

write.table(credSet99,file="../results/05_CredSet99_SNPs.txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(credSet95,file="../results/05_CredSet95_SNPs.txt",sep="\t",col.names = T,row.names = F,quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
