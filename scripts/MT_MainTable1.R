#' ---
#' title: "Main Table 1: Summary of PCSK9 locus"
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
#' **Table 1: Summary of independent SNPs per phenotype at the PCSK9 gene locus.**: 
#' 
#' * SNP - EA
#' * CADD Score
#' * Associated phenotypes
#' * Sample size
#' * Effect allele frequency
#' * beta_joint (or beta in case of oncond)
#' * P_joint (or p in case of uncond)
#' * CS 99 size
#' * Interaction T/F
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
tab1 = copy(IndepSignals)
tab1 = tab1[candidateGene =="PCSK9",]

#' ## Filt data for relevant columns
names(tab1)
names(tab1)[c(2,3,26,4,15,9,5,11,13)]

tab1 = tab1[,c(2,3,26,4,15,9,5,11,13)]
setorder(tab1,bp)
tab1

#' # Add CADD score ####
#' ***
CADD_tab = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/cadd_eigen_func_deleteriousness__2023-01-21_PCSK9_sex_strat_v2.txt")
CADD_tab = CADD_tab[markername %in% tab1$SNP,]
CADD_tab = CADD_tab[,c(1,120)]
CADD_tab = CADD_tab[!duplicated(markername)]
matched = match(tab1$SNP,CADD_tab$markername)
table(tab1$SNP == CADD_tab[matched,markername])
tab1[,CADD := CADD_tab[matched,CADD_scaled]]
tab1

#' # Add CS info ####
#' ***
load("../results/05_CredSet.RData")
credSet = credSet[candidateGene == "PCSK9"]
credSet

x = credSet[,.N,c("group","phenotype")]
x[,rsID := gsub(".::","",group)]
x[,dumID := paste(phenotype,rsID,sep="::")]
x
tab1[,dumID := paste(pheno,rsID,sep="::")]
table(duplicated(tab1$dumID))
matched = match(tab1$dumID,x$dumID)
table(tab1$dumID == x[matched,dumID])
tab1[,CS := x[matched,N]]
tab1

#' # Add interaction info ####
#' ***
load("../results/03_InteractionTests_3way.RData")
IATab_3way = IATab_3way[markername %in% tab1$SNP]
IATab_3way
load("../results/03_InteractionTests_2way.RData")
IATab_2way = IATab_2way[markername %in% tab1$SNP]
IATab_2way_sex = copy(IATab_2way)
IATab_2way_sex = IATab_2way_sex[type =="sexIA"]
IATab_2way_statin = copy(IATab_2way)
IATab_2way_statin = IATab_2way_statin[type !="sexIA"]

matched = match(tab1$SNP,IATab_3way$markername)
table(tab1$SNP == IATab_3way[matched,markername])
tab1[,IA_3way := IATab_3way[matched,IA_pval]]

matched = match(tab1$SNP,IATab_2way_sex$markername)
table(tab1$SNP == IATab_2way_sex[matched,markername])
tab1[,IA_2way_sex := IATab_2way_sex[matched,IA_pval]]

matched = match(tab1$SNP,IATab_2way_statin$markername)
table(tab1$SNP == IATab_2way_statin[matched,markername])
tab1[,IA_2way_tatin := IATab_2way_statin[matched,IA_pval]]
tab1

#' # Update Sample size ####
#' ***
#' When using GCTA-COJO, the effective sample size is calculated. But I want to report the actual sample size. 
#' 
#' I use the coloc-input data to get the sample size per SNP and trait.
#' 
load("../temp/06_GWAS_allLoci.RData")
data_GWAS = data_GWAS[markername %in% tab1$SNP]
data_GWAS[,dumID := paste(phenotype,markername,sep="::")]
tab1[,dumID := paste(pheno,SNP,sep="::")]

matched = match(tab1$dumID,data_GWAS$dumID)
table(tab1$dumID == data_GWAS[matched,dumID])
tab1[,n := data_GWAS[matched,nSamples]]
tab1

#' # Make pretty ####
#' ***
tab1[,SNP := NULL]
tab1[,dumID := NULL]
tab1[,bp := NULL]
tab1[,pheno := gsub("PCSK9_","",pheno)]
tab1[,pheno := gsub("females_","F - ",pheno)]
tab1[,pheno := gsub("males_","M - ",pheno)]
tab1[,pheno := gsub("females","F - comb",pheno)]
tab1[,pheno := gsub("males","M - comb",pheno)]
tab1[pheno == "free",pheno := "A - free"]
tab1[pheno == "treated",pheno := "A - treated"]

tab1[,freq := round(freq,3)]
tab1[,bJ := round(bJ,3)]
tab1[,pJ := signif(pJ,3)]
tab1[,IA_3way := signif(IA_3way,3)]
tab1[,IA_2way_sex := signif(IA_2way_sex,3)]
tab1[,IA_2way_tatin := signif(IA_2way_tatin,3)]
tab1[,SNP_EA := paste(rsID, refA, sep=" - ")]
tab1[,region := c(rep("A",5),"B",rep("C",6))]

tab1 = tab1[,c(14,13,3:9,11,12,10)]
tab1

#' # Save ####
#' ***
tosave4 = data.table(data = c("tab1"), 
                     SheetNames = c("MainTable1"))
excel_fn = "../tables/MainTable1.xlsx"

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

