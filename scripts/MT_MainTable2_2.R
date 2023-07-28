#' ---
#' title: "Main Table 2: Highlights of interaction tests"
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
#' New Table 2?
#' 
#' **Table 2: Summary of significant 2-way interactions.**: 
#' 
#' * cytoband	
#' * lead SNP	
#' * topPheno	
#' * type	
#' * beta	
#' * pval (topPheno)	
#' * beta	(topPheno)
#' * pval	(other trait)
#' * pIA (other trait)
#' * qIA
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
load("../results/03_InteractionTests_2way.RData")
load("../results/03_InteractionTests_3way.RData")

#' ## Filt data for PCSK9
tab2 = copy(IndepSignals)
tab2 = tab2[candidateGene !="PCSK9",]

#' ## Filt data for relevant SNPs
IATab_2way = IATab_2way[markername %in% tab2$SNP,]
IATab_3way = IATab_3way[markername %in% tab2$SNP,]

relGenes = IATab_2way[IA_hierarch_fdr5proz ==T,unique(candidateGene)]

IATab_2way = IATab_2way[candidateGene %in% relGenes,]
IATab_3way = IATab_3way[candidateGene %in% relGenes,]

tab2 = copy(IATab_2way)

names(tab2)
names(tab2)[c(1,7,8,10,14,16,20,22,27,28)]

tab2 = tab2[,c(1,7,8,10,14,16,20,22,27,28)]
tab2

#' # Add cytoband score ####
#' ***
toplist = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/topliste_2023-01-26_PCSK9_sex_strat_v2.txt")
toplist = toplist[markername %in% tab2$markername,]
matched = match(tab2$markername,toplist$markername)
table(tab2$markername == toplist[matched,markername])
tab2[,cytoband := toplist[matched,cyto]]
tab2 = tab2[,c(11,1:10)]

#' # Make pretty ####
#' ***
tab2[,bestPheno := gsub("PCSK9_","",bestPheno)]
tab2[,bestPheno := gsub("females_","F - ",bestPheno)]
tab2[,bestPheno := gsub("males_","M - ",bestPheno)]
tab2[,bestPheno := gsub("females","F - comb",bestPheno)]
tab2[,bestPheno := gsub("males","M - comb",bestPheno)]
tab2[bestPheno == "free",bestPheno := "A - free"]
tab2[bestPheno == "treated",bestPheno := "A - treated"]

tab2[,markername := gsub(":.*","",markername)]

tab2[,trait1_beta := round(trait1_beta,3)]
tab2[,trait2_beta := round(trait2_beta,3)]
tab2[,trait1_pval := signif(trait1_pval,3)]
tab2[,trait2_pval := signif(trait2_pval,3)]
tab2[,IA_pval := signif(IA_pval,3)]
tab2[,IA_pval_adj := signif(IA_pval_adj,3)]

tab2

#' # Save ####
#' ***
tosave4 = data.table(data = c("tab2"), 
                     SheetNames = c("MainTable2"))
excel_fn = "../tables/MainTable2_230324.xlsx"

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


