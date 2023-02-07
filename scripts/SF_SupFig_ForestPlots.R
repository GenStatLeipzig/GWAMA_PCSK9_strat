#' ---
#' title: "Supplemental Figures: Forest Plots"
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
#' Forest plots for all loci:
#' 
#' * PCSK9: 
#'    * 4 SNPs with all studies in all traits
#' * other SNPs:
#'    * 12 SNPs with all studies in best trait only
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")
source("../helperFunctions/plotForest1.R")
source("../helperFunctions/plotForest2.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/05_GCTA_COJO.RData")

ToDoList = copy(IndepSignals)
ToDoList = ToDoList[,c(1:3,26,22:25,21,15)]
ToDoList = ToDoList[!duplicated(SNP),]
ToDoList

myTab<-fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/step80_meta_singlestudyresults_2023-01-26_PCSK9_sex_strat_v2.txt")
myTab<-myTab[markername %in% ToDoList$SNP,]
setorder(myTab,chr,pos)

#' # Select data ####
#' ***
PCSK9_SNPs = copy(myTab)
PCSK9_SNPs = PCSK9_SNPs[markername %in% ToDoList[1:4,SNP]]
PCSK9_SNPs[,candidateGene := "PCSK9"]

Other_SNPs = copy(myTab)
Other_SNPs[,dumID:= paste(pheno,markername,sep="::")]
ToDoList[,dumID := paste(pheno,SNP,sep="::")]
Other_SNPs = Other_SNPs[dumID %in% ToDoList[5:16,dumID]]
table(Other_SNPs$dumID == ToDoList[5:16,dumID])
Other_SNPs[,candidateGene := ToDoList[5:16,candidateGene]]
Other_SNPs[,dumID:=NULL]

tab1 = rbind(PCSK9_SNPs,Other_SNPs)


setnames(tab1, c( "CochransQ" ,"pCochransQ"), c( "cochQ", "cochQpval"), skip_absent = T)
tab2 = melt(tab1, 
            id.vars = c("markername","pheno","betaFEM","seFEM","pFEM","cochQ","cochQpval","I2","betaREM","seREM","pREM"), 
            measure.vars = patterns("^beta\\.|^se\\.|^p\\."), 
            variable.factor = F )
tab2 = tab2[!is.na(value),]
tab2[grepl("LIFE_Adult",variable),cohort :="LIFE-Adult"]
tab2[grepl("LIFE_Heart",variable),cohort :="LIFE-Heart"]
tab2[grepl("LURIC",variable),cohort :="LURIC"]
tab2[grepl("TWINGENE",variable),cohort :="TwinGene"]
tab2[grepl("CLEANED.males_treated",variable),cohort :="GWAMA 1: males treated"]
tab2[grepl("CLEANED.females_treated",variable),cohort :="GWAMA 1: females treated"]
tab2[grepl("CLEANED.males_free",variable),cohort :="GWAMA 1: males free"]
tab2[grepl("CLEANED.females_free",variable),cohort :="GWAMA 1: females free"]
table(tab2$cohort)
tab2[,stats := sapply(str_split(variable, "^beta\\.|^se\\.|^p\\."), '[', 2)]
tab2[,stats := str_replace_all(variable, paste0(".",stats), "")]
tab2[grepl("male",variable),sex := "males"]
tab2[grepl("femal",variable),sex := "females"]
tab2[grepl("Fem",variable),sex :="females"]
tab2[grepl("Mal",variable),sex :="males"]
tab2[,table(sex,is.na(sex))]
tab2[grepl("treated",variable),statin := "treated"]
tab2[grepl("free",variable),statin := "free"]
tab2[grepl("TwinGene",cohort),statin :="free"]
tab2[,table(statin,is.na(statin))]
tab2[,table(statin,sex)]


#' # Plot PCSK9: rs11591147
#' 
tiff(filename = "../figures/SupplementalFigure_ForestPlots_PCSK9_rs11591147.tiff", 
     width = 2400, height = 1600, res = 200, compression = 'lzw')

par(mfrow=c(4,2))
plotForest1(mysnpstats = tab1[1,],data = tab2)
plotForest1(mysnpstats = tab1[2,],data = tab2)
plotForest1(mysnpstats = tab1[5,],data = tab2)
plotForest1(mysnpstats = tab1[6,],data = tab2)
plotForest1(mysnpstats = tab1[3,],data = tab2)
plotForest1(mysnpstats = tab1[7,],data = tab2)
plotForest1(mysnpstats = tab1[4,],data = tab2)
plotForest1(mysnpstats = tab1[8,],data = tab2)

dev.off()


#' # Plot PCSK9: rs28385704
#' 
tiff(filename = "../figures/SupplementalFigure_ForestPlots_PCSK9_rs28385704.tiff", 
     width = 2400, height = 1600, res = 200, compression = 'lzw')

par(mfrow=c(4,2))
plotForest1(mysnpstats = tab1[9,],data = tab2)
plotForest1(mysnpstats = tab1[10,],data = tab2)
plotForest1(mysnpstats = tab1[13,],data = tab2)
plotForest1(mysnpstats = tab1[14,],data = tab2)
plotForest1(mysnpstats = tab1[11,],data = tab2)
plotForest1(mysnpstats = tab1[15,],data = tab2)
plotForest1(mysnpstats = tab1[12,],data = tab2)
plotForest1(mysnpstats = tab1[16,],data = tab2)

dev.off()

#' # Plot PCSK9: rs553741
#' 
tiff(filename = "../figures/SupplementalFigure_ForestPlots_PCSK9_rs553741.tiff", 
     width = 2400, height = 1600, res = 200, compression = 'lzw')

par(mfrow=c(4,2))
plotForest1(mysnpstats = tab1[17,],data = tab2)
plotForest1(mysnpstats = tab1[18,],data = tab2)
plotForest1(mysnpstats = tab1[21,],data = tab2)
plotForest1(mysnpstats = tab1[22,],data = tab2)
plotForest1(mysnpstats = tab1[19,],data = tab2)
plotForest1(mysnpstats = tab1[23,],data = tab2)
plotForest1(mysnpstats = tab1[20,],data = tab2)
plotForest1(mysnpstats = tab1[24,],data = tab2)

dev.off()

#' # Plot PCSK9: rs693668
#' 
tiff(filename = "../figures/SupplementalFigure_ForestPlots_PCSK9_rs693668.tiff", 
     width = 2400, height = 1600, res = 200, compression = 'lzw')

par(mfrow=c(4,2))
plotForest1(mysnpstats = tab1[25,],data = tab2)
plotForest1(mysnpstats = tab1[26,],data = tab2)
plotForest1(mysnpstats = tab1[29,],data = tab2)
plotForest1(mysnpstats = tab1[30,],data = tab2)
plotForest1(mysnpstats = tab1[27,],data = tab2)
plotForest1(mysnpstats = tab1[31,],data = tab2)
plotForest1(mysnpstats = tab1[28,],data = tab2)
plotForest1(mysnpstats = tab1[32,],data = tab2)

dev.off()


#' # Plot Lipids: 
#' 34,11,3,7,8,9
#' 
tiff(filename = "../figures/SupplementalFigure_ForestPlots_lipids.tiff", 
     width = 2400, height = 1200, res = 200, compression = 'lzw')

par(mfrow=c(3,2))
plotForest2(mysnpstats = tab1[32+2,],data = tab2)
plotForest2(mysnpstats = tab1[32+11,],data = tab2)
plotForest2(mysnpstats = tab1[32+3,],data = tab2)
plotForest2(mysnpstats = tab1[32+7,],data = tab2)
plotForest2(mysnpstats = tab1[32+8],data = tab2)
plotForest2(mysnpstats = tab1[32+9,],data = tab2)

dev.off()

#' # Plot other hits: 
#' 4,5,10,1,6,12
#' 
tiff(filename = "../figures/SupplementalFigure_ForestPlots_notLipids.tiff", 
     width = 2400, height = 1200, res = 200, compression = 'lzw')

par(mfrow=c(3,2))
plotForest2(mysnpstats = tab1[32+4,],data = tab2)
plotForest2(mysnpstats = tab1[32+5,],data = tab2)
plotForest2(mysnpstats = tab1[32+10,],data = tab2)
plotForest2(mysnpstats = tab1[32+1,],data = tab2)
plotForest2(mysnpstats = tab1[32+6,],data = tab2)
plotForest2(mysnpstats = tab1[32+12,],data = tab2)

dev.off()


#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
