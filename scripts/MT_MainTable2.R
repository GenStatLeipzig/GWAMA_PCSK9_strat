#' ---
#' title: "Main Table 2: Interaction Test results"
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
#' **Table 2: Significant interactions of the indenpendent SNPs of the valid loci**: 
#' 
#' * Cytoband
#' * Candidate Gene
#' * SNP
#' * Best associated phenotype
#' * Type (3way, sexIA, statinIA)
#' * Difference 
#' * SE
#' * p-value
#' * q-value (after hierFDR!)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/04_IATest_2way_filtered.RData")

#' ## Filt data for significant IAs
tab2 = copy(IATab_filtered)
tab2 = tab2[IA_hierarch_fdr5proz == T,]

#' ## Filt data for relevant columns
names(tab2)
names(tab2)[c(6,1,7,8,10,11,13,14)]

tab2 = tab2[,c(6,1,7,8,10,11,13,14)]
tab2

#' # Add cytoband score ####
#' ***
toplist = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/topliste_2023-07-26_PCSK9_strat.txt")
toplist = toplist[markername %in% tab2$markername,]
matched = match(tab2$markername,toplist$markername)
table(tab2$markername == toplist[matched,markername])
tab2[,cytoband := toplist[matched,cyto]]
tab2 = tab2[,c(9,1:8)]
tab2

#' # Make pretty ####
#' ***
setnames(tab2,"markername","SNP")
tab2[,SNP := gsub(":.*","",SNP)]

tab2[,bestPheno := gsub("PCSK9_","",bestPheno)]
tab2[,bestPheno := gsub("females_","F - ",bestPheno)]
tab2[,bestPheno := gsub("males_","M - ",bestPheno)]
tab2[,bestPheno := gsub("males","M",bestPheno)]

tab2[,type := gsub("IA","",type)]

tab2[,IA_diff := round(IA_diff,3)]
tab2[,IA_SE := round(IA_SE,3)]
tab2[,IA_pval := signif(IA_pval,3)]
tab2[,IA_pval_adj := signif(IA_pval_adj,3)]
setnames(tab2,"IA_pval_adj","IA_qval")

tab2

setorder(tab2,type,IA_qval)

#' # Save ####
#' ***
tosave4 = data.table(data = c("tab2"), 
                     SheetNames = c("MainTable2"))
excel_fn = "../tables/MainTable2_230907.xlsx"

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


