#' ---
#' title: "Look-up of sex-biased gene expression"
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
#' I just want to check the sex-stratified eQTL data for sex-biased eQTLs
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_forostar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

#' # Check eQTLs ####
#' ***
tab = fread(paste0(path_GTExv8_sexStrat, "/GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt"))
tab = tab[hugo_gene_id  %in% c("PCSK9","APOB","KHDRBS2" ,"PRKAG2","MARCHF8",
                               "ALOX5","JMJD1C","FADS1", "FADS2","SLCO1B3",
                               "KSR2", "NOS1","HPR", "HP","TM6SF2","CNKSR2")]
tab2 = tab[pval_nominal_sb<0.05,]
tab2
tab2 = tab2[pvals.corrected<0.05,]
tab2

#' # Check gene expression ####
#' ***
tab4 = fread(paste0(path_GTExv8_sexStrat, "/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt"))
tab4 = tab4[gene  %in% tab$ensembl_gene_id,]
matched = match(tab4$gene,tab$ensembl_gene_id)
table(tab4$gene == tab[matched,ensembl_gene_id])
tab4[,hugo_gene_id := tab[matched,hugo_gene_id]]
tab4[,.N,by=hugo_gene_id]
tab4

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
