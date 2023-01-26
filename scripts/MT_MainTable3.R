#' ---
#' title: "Main Table 3: Bidirectional MR"
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
#' **Table 3: Main results of the bidirectional Mendelian Randomization.**: 
#' 
#' * Model (females, males, females_free, males_free, free)
#' * Direction 1: PCSK9 --> LDLC 
#'    * beta_IV
#'    * se_IV
#'    * p_IV
#'    * Q
#' * Direction 2: LDLC --> PCSK9 
#'    * beta_IV
#'    * se_IV
#'    * p_IV
#'    * Q
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/07_MR_LDLC_PCSK9_summary.RData")
load("../results/07_MR_PCSK9_LDLC_summary.RData")

MRTab_PCSK9SNPs
myMRTab_meta

#' ## Filt data for relevant traits
tab3 = copy(MRTab_PCSK9SNPs)
tab3 = tab3[!is.na(Q_IV)]
myMRTab_meta = myMRTab_meta[N==20,]
matched = match(tab3$phenotype,myMRTab_meta$phenotype)
myMRTab_meta = myMRTab_meta[matched]

tab3[,beta_IV_2 := myMRTab_meta$beta_IV]
tab3[,SE_IV_2 := myMRTab_meta$SE_IV]
tab3[,pval_IV_2 := myMRTab_meta$pval_IV]
tab3[,Q_IV_2 := myMRTab_meta$Q_IV]

tab3[,N:=NULL]
tab3[,comment:=NULL]
tab3[,outcome:=NULL]
tab3[,pvalQ_IV:=NULL]

tab3[,phenotype:=c("females (statin-combined)",
                   "females (statin-free)",
                   "statin-free (sex-combined)",
                   "males (statin-combined)",
                   "males (statin-free)")]

tab3 = tab3[c(1,4,2,5,3)]
setnames(tab3,"phenotype","model")
tab3

#' # Round estimates ####
#' ***
tab3[,beta_IV := round(beta_IV,3)]
tab3[,beta_IV_2 := round(beta_IV_2,3)]
tab3[,SE_IV := round(SE_IV,3)]
tab3[,SE_IV_2 := round(SE_IV_2,3)]
tab3[,Q_IV := round(Q_IV,3)]
tab3[,Q_IV_2 := round(Q_IV_2,3)]

tab3[,pval_IV := signif(pval_IV,3)]
tab3[,pval_IV_2 := signif(pval_IV_2,3)]

tab3

#' # Save ####
#' ***
tosave4 = data.table(data = c("tab3"), 
                     SheetNames = c("MainTable3"))
excel_fn = "../tables/MainTable3.xlsx"

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


