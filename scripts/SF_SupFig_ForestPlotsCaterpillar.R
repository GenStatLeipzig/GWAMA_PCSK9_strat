#' ---
#' title: "Supplemental Figures: Forest Plots (Caterpillar)"
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
#'    * 3 SNPs with all subgroups
#' * interaction hit lead SNPs:
#'    * 7 SNPs with four double stratified subgroups / all subgroups
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")
source("../helperFunctions/plotForest3.R")
source("../helperFunctions/plotForest4.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/05_GCTA_COJO.RData")

ToDoList = copy(IndepSignals)
ToDoList = ToDoList[,c(1:3,26,22:25,21,15)]
ToDoList = ToDoList[!duplicated(SNP),]
ToDoList

load("../results/03_InteractionTests_2way.RData")
IATab_2way = IATab_2way[markername %in% ToDoList$SNP,]

load("../results/03_InteractionTests_3way.RData")
IATab_3way = IATab_3way[markername %in% ToDoList$SNP,]

SigIAGenes = IATab_2way[IA_hierarch_fdr5proz == T, unique(candidateGene)]

ToDoList = ToDoList[candidateGene %in% c(SigIAGenes,"PCSK9"),]
ToDoList = ToDoList[SNP != "rs553741:55520408:G:C",]

myTab<-fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/step80_meta_singlestudyresults_2023-01-26_PCSK9_sex_strat_v2.txt")
myTab<-myTab[markername %in% ToDoList$SNP,]
setorder(myTab,chr,pos)
myNames = names(myTab)[c(1:16,53:55,58)]
myNames

colsOut<-setdiff(colnames(myTab),myNames)
myTab[,get("colsOut"):=NULL]
dim(myTab)

myTab[,table(pheno)]
myTab[,phenotype := gsub("PCSK9_","",pheno)]
myTab[,phenotype := gsub("females_","F - ",phenotype)]
myTab[,phenotype := gsub("males_","M - ",phenotype)]
myTab[,phenotype := gsub("females","F - comb",phenotype)]
myTab[,phenotype := gsub("males","M - comb",phenotype)]
myTab[!grepl("-",phenotype),phenotype := paste0("A - ",phenotype)]
myTab[,table(phenotype)]


#' # Plot PCSK9 SNPs ####
#' ***
tiff(filename = "../figures/SupplementalFigure_ForestPlots_Caterpillar_PCSK9_combined.tiff", 
     width = 1400, height = 1400, res = 200, compression = 'lzw')

par(mfrow=c(1,1))
plotForest4(mysnpstats = ToDoList[1,],data = myTab[chr==1])
dev.off()

tiff(filename = "../figures/SupplementalFigure_ForestPlots_Caterpillar_PCSK9_perSNP.tiff", 
     width = 1000, height = 1350, res = 200, compression = 'lzw')

par(mfrow=c(3,1))
plotForest3(mysnpstats = ToDoList[1,],data = myTab)
plotForest3(mysnpstats = ToDoList[2,],data = myTab)
plotForest3(mysnpstats = ToDoList[3,],data = myTab)
dev.off()

#' # Plot interaction Hits ####
#' ***
#' I group the 7 SNPs like mentioned in the paper:
#' 
#' * female-specific & treatment-specific: MACROD2, SLCO1B3, FOSL2, and KHDRBS2
#' * multiple interactions: SASH1, PRKAG2, and gene desert
#' 
#' ## Round 1: all 8 subgroups ####
tiff(filename = "../figures/SupplementalFigure_ForestPlots_Caterpillar_2wayIA.tiff", 
     width = 1000, height = 1800, res = 200, compression = 'lzw')

par(mfrow=c(4,1))

plotForest3(mysnpstats = ToDoList[10,],data = myTab)
plotForest3(mysnpstats = ToDoList[8,],data = myTab)
plotForest3(mysnpstats = ToDoList[4,],data = myTab)
plotForest3(mysnpstats = ToDoList[5,],data = myTab)

dev.off()

tiff(filename = "../figures/SupplementalFigure_ForestPlots_Caterpillar_3wayIA.tiff", 
     width = 1000, height = 1350, res = 200, compression = 'lzw')

par(mfrow=c(3,1))

plotForest3(mysnpstats = ToDoList[6,],data = myTab)
plotForest3(mysnpstats = ToDoList[7,],data = myTab)
plotForest3(mysnpstats = ToDoList[9,],data = myTab)

dev.off()

#' ## Round 2: 4 double-stratified subgroups ####

uniquePhenos = myTab[,unique(phenotype)]
relPhenos = uniquePhenos[!grepl("A - ",uniquePhenos)]
relPhenos = relPhenos[!grepl("- comb",relPhenos)]
relPhenos

tiff(filename = "../figures/SupplementalFigure_ForestPlots_Caterpillar_2wayIA_filt.tiff", 
     width = 1000, height = 1350, res = 200, compression = 'lzw')

par(mfrow=c(4,1))

plotForest3(mysnpstats = ToDoList[10,],data = myTab[phenotype %in% relPhenos])
plotForest3(mysnpstats = ToDoList[8,],data = myTab[phenotype %in% relPhenos])
plotForest3(mysnpstats = ToDoList[4,],data = myTab[phenotype %in% relPhenos])
plotForest3(mysnpstats = ToDoList[5,],data = myTab[phenotype %in% relPhenos])

dev.off()

tiff(filename = "../figures/SupplementalFigure_ForestPlots_Caterpillar_3wayIA_filt.tiff", 
     width = 1000, height = 1000, res = 200, compression = 'lzw')

par(mfrow=c(3,1))

plotForest3(mysnpstats = ToDoList[6,],data = myTab[phenotype %in% relPhenos])
plotForest3(mysnpstats = ToDoList[7,],data = myTab[phenotype %in% relPhenos])
plotForest3(mysnpstats = ToDoList[9,],data = myTab[phenotype %in% relPhenos])

dev.off()



#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
