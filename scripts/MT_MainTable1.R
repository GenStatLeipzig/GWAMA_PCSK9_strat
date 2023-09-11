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
#' **Table 1: Summary of lead SNPs (best association per SNP)**: 
#' 
#' * cytoband
#' * lead SNP
#' * effect allele
#' * EAF
#' * beta
#' * pvalue
#' * explained variance
#' * best associated phenotype
#' * candidate gene
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../temp/04_IATest_input.RData")
result.5 = copy(result.2)
setorder(result.5,pval)
result.5 = result.5[!duplicated(markername),]

gwas_annot = fread(paste(path_GenStatPipeline,"synopsis/topliste_tabdelim/topliste_2023-07-26_PCSK9_strat.txt"))
gwas_annot = gwas_annot[markername %in% result.5$markername]

table(gwas_annot$markername == result.5$markername)
result.5[,cytoband := gwas_annot$cyto]
result.5[,nearbyGene := gwas_annot$nearestgenes]

#' ## Filt coulmns
tab1 = copy(result.5)
names(tab1)
tab1[,r2 := beta^2 / (beta^2 + nSamples*SE^2)]
tab1 = tab1[,c(26,18,4,6,10,12,28,16,27)]
tab1

#' # Make pretty ####
#' ***
names(tab1)
tab1[,phenotype := gsub("PCSK9_","",phenotype)]
tab1[,phenotype := gsub("females_","F - ",phenotype)]
tab1[,phenotype := gsub("males_","M - ",phenotype)]
tab1[,phenotype := gsub("males","M",phenotype)]

tab1[,EAF := round(EAF,3)]
tab1[,beta := round(beta,3)]
tab1[,pval := signif(pval,3)]
tab1[,r2 := 100*r2]
tab1[,r2 := round(r2,2)]

tab1[,nearbyGene := gsub("[(].*","",nearbyGene)]
setnames(tab1,"nearbyGene","nearestGene")
tab1[,candidateGene := result.5$candidateGene]

#' # Save ####
#' ***
tag = format(Sys.time(), "%Y-%m-%d")
tag2 = gsub("2023-","23-",tag)
tag2 = gsub("-","",tag2)

tosave4 = data.table(data = c("tab1"), 
                     SheetNames = c("MainTable1"))
excel_fn = paste0("../tables/MainTable1_",tag2,".xlsx")

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

