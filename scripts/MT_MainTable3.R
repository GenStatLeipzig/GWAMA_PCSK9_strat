#' ---
#' title: "Main Table 3: Stratified MR"
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
#' **Table 3: Main results of the stratified Mendelian Randomization.**: 
#' 
#' I want the eight subgroups for setting all (fixed) and without rs11583680 (beta, se, Q, p(Q) each)
#' 
#' => eight rows, nine columns
#' 
#' **Table 4: Main results of the IA Test within MR**
#' 
#' I want the six tests for setting all (fixed) (2x trait, theta, se; dif, dif p)
#' 
#' => six rows, 8 columns 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/06_MR_WaldEstimates_UKBB_230907.RData")
load("../results/06_MR_IVWEstimates_UKBB_230907.RData")
load("../results/06_MR_InteractionTest_UKBB_230907.RData")

#' # Table 3 ####
#' ***
F1 = myMRTab[,min(FStat.PCSK9),by = phenotype]
F2 = myMRTab[SNP != "rs11583680:55505668:C:T",min(FStat.PCSK9),by = phenotype]

tab3 = copy(myMRTab_results)
tab3 = tab3[grepl("all",setting),]
tab3 = tab3[,c(1,4:5,7,8)]
table(tab3$phenotype == F1$phenotype)
tab3[,minFStat_all := F1$V1]
tab3 = tab3[,c(1,6,2:5)]

tab32 = copy(myMRTab_results)
tab32 = tab32[grepl("rs11583680",setting),]
tab32 = tab32[,c(4:5,7,8)]

tab3 = cbind(tab3,tab32)

#' make pretty
tab3[,phenotype := gsub("PCSK9_","",phenotype)]
tab3[,phenotype := gsub("females_","W - ",phenotype)]
tab3[,phenotype := gsub("females","W",phenotype)]
tab3[,phenotype := gsub("males_","M - ",phenotype)]
tab3[,phenotype := gsub("males","M",phenotype)]

names(tab3)[2:6] = gsub("IVW","all",names(tab3)[2:6]) 
names(tab3)[7:10] = gsub("IVW","woRs11583680",names(tab3)[7:10]) 

tab3[,beta_all := round(beta_all,2)]
tab3[,SE_all := round(SE_all,2)]
tab3[,HeteroStat_all := round(HeteroStat_all,1)]
tab3[,HeteroStat_pval_all := signif(HeteroStat_pval_all,3)]

tab3[,beta_woRs11583680 := round(beta_woRs11583680,2)]
tab3[,SE_woRs11583680 := round(SE_woRs11583680,2)]
tab3[,HeteroStat_woRs11583680 := round(HeteroStat_woRs11583680,1)]
tab3[,HeteroStat_pval_woRs11583680 := signif(HeteroStat_pval_woRs11583680,3)]

tab3[,minFStat_all := round(minFStat_all,1)]
tab3

tag = format(Sys.time(), "%Y-%m-%d")
tag2 = gsub("2023-","23-",tag)
tag2 = gsub("-","",tag2)

tosave4 = data.table(data = c("tab3"), 
                     SheetNames = c("MainTable3"))
excel_fn = paste0("../tables/MainTable3_",tag2,".xlsx")

WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Table 4 ####
#' ***
tab4 = copy(IATab_IVW)
tab4 = tab4[setting == "all (fixed)",]
tab4 = tab4[,c(2:4,6:8,10,11,13)]

tab4[,trait1 := gsub("females_","W - ",trait1)]
tab4[,trait1 := gsub("females","W",trait1)]
tab4[,trait1 := gsub("males_","M - ",trait1)]
tab4[,trait1 := gsub("males","M",trait1)]

tab4[,trait2 := gsub("females_","W - ",trait2)]
tab4[,trait2 := gsub("females","W",trait2)]
tab4[,trait2 := gsub("males_","M - ",trait2)]
tab4[,trait2 := gsub("males","M",trait2)]

tab4[,trait1_beta := round(trait1_beta,2)]
tab4[,trait1_SE := round(trait1_SE,2)]
tab4[,trait2_beta := round(trait2_beta,2)]
tab4[,trait2_SE := round(trait2_SE,2)]
tab4[,IA_diff := round(IA_diff,2)]
tab4[,IA_SE := round(IA_SE,2)]
tab4[,IA_pval := signif(IA_pval,3)]

tag = format(Sys.time(), "%Y-%m-%d")
tag2 = gsub("2023-","23-",tag)
tag2 = gsub("-","",tag2)

tosave4 = data.table(data = c("tab4"), 
                     SheetNames = c("MainTable4"))
excel_fn = paste0("../tables/MainTable4_",tag2,".xlsx")

WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


