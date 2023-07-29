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

#' # Get Credible Sets for PCSK9 locus ####
#' ***
load("../results/03_GCTA_COJO.RData")

ToDoList = data.table(NR = 1:32)

ToDoList[,statistic := list.files(path = "../results/03_GCTA_COJO_cond/",pattern = ".cma.cojo")]
ToDoList[,statistic_path := paste0("../results/03_GCTA_COJO_cond/",statistic)]

ToDoList[,pheno := gsub("_signal.*","",statistic)]
ToDoList[,indepSNP := gsub(".*_signal_","",statistic)]
ToDoList[,indepSNP := gsub(".cma.cojo","",indepSNP)]

dumTab = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = copy(ToDoList)
  myRow = myRow[i,]
  credSet = getCredibleSet(filename = myRow$statistic_path,
                           type = "cond",
                           CHR = 1,
                           region_start = IndepSignals$region_start[1],
                           region_end = IndepSignals$region_end[1],
                           NR = myRow$NR,
                           SNP = myRow$indepSNP,
                           pheno = myRow$pheno,
                           candidateGene = "PCSK9")
  
  credSet[,phenotype := myRow$pheno]
  credSet[,locus := myRow$NR]
  credSet[,leadSNP := myRow$indepSNP]
  credSet[,candidateGene := "PCSK9"]
  credSet[,type := "cond"]
  credSet[,rank := 1:dim(credSet)[1]]
  credSet[,flag := F]
  credSet[1:min(which(credSet[,SumProb]>0.95)),flag:=T]
  credSet
}

credSet_PCSK9 = rbindlist(dumTab)
PCSK9_leadSNPs = unique(credSet_PCSK9$leadSNP)
check1 = credSet_PCSK9[,.N,by=c("phenotype","leadSNP")]
check1[leadSNP == PCSK9_leadSNPs[1]]
check1[leadSNP == PCSK9_leadSNPs[2]]
check1[leadSNP == PCSK9_leadSNPs[3]]
check1[leadSNP == PCSK9_leadSNPs[4]]

check2 = copy(credSet_PCSK9)
setorder(check2,-PostProb)
check2 = check2[!duplicated(SNP),]
check2[,table(phenotype,pval<0.05)]
check2[PostProb>0.1,]

#' # Get Credible Sets for other locus ####
#' ***
IndepSignals = IndepSignals[candidateGene != "PCSK9",]
IndepSignals[, input_CS := paste0("../data/SumStat_",pheno,"_230728.txt.gz")]

dumTab = foreach(i=1:dim(IndepSignals)[1])%do%{
  #i=1
  myRow = copy(IndepSignals)
  myRow = myRow[i,]
  
  credSet = getCredibleSet(filename = myRow$input_CS,
                           type = "GWAMA",
                           CHR = myRow$Chr,
                           region_start = myRow$region_start,
                           region_end = myRow$region_end,
                           NR = i,
                           SNP = myRow$SNP,
                           pheno = myRow$pheno,
                           candidateGene = myRow$candidateGene)
  
  credSet[,phenotype := myRow$pheno]
  credSet[,locus := i]
  credSet[,leadSNP := myRow$SNP]
  credSet[,candidateGene := myRow$candidateGene]
  credSet[,type := "GWAMA"]
  credSet[,rank := 1:dim(credSet)[1]]
  credSet[,flag := F]
  credSet[1:min(which(credSet[,SumProb]>0.95)),flag:=T]
  credSet
}

credSet_other = rbindlist(dumTab)
credSet_other[,leadSNP := gsub(":.*","",leadSNP)]
leadSNPs = unique(credSet_other$leadSNP)
check3 = credSet_other[,.N,by=candidateGene]
check3

#' # Some checks and plots ####
#' ***
#' ## PCSK9 locus
credSet_PCSK9[,group := paste(phenotype,leadSNP,sep=" ")]
credSet_PCSK9[,group := gsub("PCSK9_","",group)]
credSet_PCSK9[,group := gsub("_"," ",group)]
credSet_PCSK9[,table(group)]

myPlot1 = ggplot(credSet_PCSK9[leadSNP == PCSK9_leadSNPs[1]], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/03_CredSet_PCSK9_signal_rs11583680.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot1
dev.off()

myPlot1 = ggplot(credSet_PCSK9[leadSNP == PCSK9_leadSNPs[2]], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/03_CredSet_PCSK9_signal_rs11591147.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot1
dev.off()

myPlot1 = ggplot(credSet_PCSK9[leadSNP == PCSK9_leadSNPs[3]], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/03_CredSet_PCSK9_signal_rs2495491.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot1
dev.off()

myPlot1 = ggplot(credSet_PCSK9[leadSNP == PCSK9_leadSNPs[4]], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/03_CredSet_PCSK9_signal_rs693668.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot1
dev.off()

#' ## Other loci
credSet_other[,group := paste(phenotype,candidateGene,sep=" ")]
credSet_other[,group := gsub("PCSK9_","",group)]
credSet_other[,group := gsub("_"," ",group)]
credSet_other[,table(group)]

myPlot2 = ggplot(credSet_other[candidateGene %in% c("ALOX5","SLCO1B3","NOS1","PRKAG2","KHDRBS2")], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/03_CredSet_other_signals1.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot2
dev.off()

myPlot2 = ggplot(credSet_other[candidateGene %in% c("APOB","FADS1","HP","JMJD1C","TM6SF2")], aes(x=rank, y=SumProb)) + 
  geom_line() + 
  facet_wrap(facets = "group",scales = "free_x") + 
  ylab("sum of Posterior Probability") + 
  xlab("position in CS") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=14))

tiff(filename = "../results/03_CredSet_other_signals2.tif", 
     width = 1200, height = 900, res=125, compression = 'lzw')
myPlot2
dev.off()

#' # Save ####
#' ***
names(credSet_PCSK9)
setnames(credSet_PCSK9,"SNP","snp")
save(credSet_PCSK9, file="../results/03_CredSet_PCSK9.RData")

setnames(credSet_other,"SNP","snp")
save(credSet_other, file="../results/03_CredSet_otherLoci.RData")

credSet99 = copy(credSet_PCSK9)
credSet99 = credSet99[,c(1,17,8,9)] 

credSet95 = copy(credSet_PCSK9)
credSet95 = credSet95[flag == T,] 
credSet95 = credSet95[,c(1,17,8,9)] 

write.table(credSet99,file="../results/03_CredSet99_PCSK9_SNPs.txt",sep="\t",col.names = T,row.names = F,quote = F)
write.table(credSet95,file="../results/03_CredSet95_PCSK9_SNPs.txt",sep="\t",col.names = T,row.names = F,quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
