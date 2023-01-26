#' ---
#' title: "Main Table 2: Summary of PCSK9 locus"
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
#' **Table 1: Additional SNPs with at least suggestive significance in our GWAMA.**: 
#' 
#' * Cytoband
#' * Candidate Gene
#' * SNP - EA
#' * Associated phenotypes
#' * Sample size
#' * Effect allele frequency
#' * beta
#' * P
#' * Interaction T/F/T
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/05_GCTA_COJO.RData")

#' ## Filt data for PCSK9
tab2 = copy(IndepSignals)
tab2 = tab2[candidateGene !="PCSK9",]

#' ## Filt data for relevant columns
names(tab2)
names(tab2)[c(2,3,21,26,4,15,9,5,6,8)]

tab2 = tab2[,c(2,3,21,26,4,15,9,5,6,8)]
tab2

#' # Add cytoband score ####
#' ***
toplist = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/topliste_2023-01-21_PCSK9_sex_strat_v2.txt")
toplist = toplist[markername %in% tab2$SNP,]
matched = match(tab2$SNP,toplist$markername)
table(tab2$SNP == toplist[matched,markername])
tab2[,cytoband := toplist[matched,cyto]]
tab2[,n := toplist[matched,topn]]
tab2

#' # Add interaction info ####
#' ***
load("../results/03_InteractionTests_3way.RData")
IATab_3way = IATab_3way[markername %in% tab2$SNP]
IATab_3way
load("../results/03_InteractionTests_2way.RData")
IATab_2way = IATab_2way[markername %in% tab2$SNP]
IATab_2way_sex = copy(IATab_2way)
IATab_2way_sex = IATab_2way_sex[type =="sexIA"]
IATab_2way_statin = copy(IATab_2way)
IATab_2way_statin = IATab_2way_statin[type !="sexIA"]

matched = match(tab2$SNP,IATab_3way$markername)
table(tab2$SNP == IATab_3way[matched,markername])
tab2[,IA_3way := IATab_3way[matched,IA_pval]]

matched = match(tab2$SNP,IATab_2way_sex$markername)
table(tab2$SNP == IATab_2way_sex[matched,markername])
tab2[,IA_2way_sex := IATab_2way_sex[matched,IA_pval]]

matched = match(tab2$SNP,IATab_2way_statin$markername)
table(tab2$SNP == IATab_2way_statin[matched,markername])
tab2[,IA_2way_tatin := IATab_2way_statin[matched,IA_pval]]
tab2

#' # Make pretty ####
#' ***
tab2[,SNP := NULL]
tab2[,bp := NULL]
tab2[,pheno := gsub("PCSK9_","",pheno)]
tab2[,pheno := gsub("females_","F - ",pheno)]
tab2[,pheno := gsub("males_","M - ",pheno)]
tab2[,pheno := gsub("females","F - comb",pheno)]
tab2[,pheno := gsub("males","M - comb",pheno)]
tab2[pheno == "free",pheno := "A - free"]
tab2[pheno == "treated",pheno := "A - treated"]

tab2[,freq := round(freq,3)]
tab2[,b := round(b,3)]
tab2[,p := signif(p,3)]
tab2[,IA_3way := signif(IA_3way,3)]
tab2[,IA_2way_sex := signif(IA_2way_sex,3)]
tab2[,IA_2way_tatin := signif(IA_2way_tatin,3)]
tab2[,SNP_EA := paste(rsID, refA, sep=" - ")]
tab2[,region := c(rep("A",5),"B",rep("C",6))]

tab2 = tab2[,c(9,1,13,4:8,10:12)]
tab2

#' # Save ####
#' ***
tosave4 = data.table(data = c("tab2"), 
                     SheetNames = c("MainTable2"))
excel_fn = "../tables/MainTable2.xlsx"

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


