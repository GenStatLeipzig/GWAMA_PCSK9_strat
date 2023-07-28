#' ---
#' title: "Get summary statistics"
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
#' Load the step10 data objects from our pipeline and save as .txt as it would be uploaded to zenodo or LHA. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

#' # Create ToDo ####
#' ***
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../../2307_GWAMA/06_Annotation2/results/",pattern = "step10")]
ToDoList[,statistic_path := paste0("../../2307_GWAMA/06_Annotation2/results/",statistic)]

ToDoList[,pheno := gsub("step10_2_primaryANDsecondary_","",statistic)]
ToDoList[,pheno := gsub(".RData","",pheno)]

ToDoList[,sex := c(rep("female",3),"combined",rep("male",3),"combined")]
ToDoList[,statin := c("free","treated","combined","free","free","treated","combined","treated")]

ToDoList

#' # Loop over ToDo ####
#' ***

dumTab<-foreach(i=1:dim(ToDoList)[1]) %do% {
  #i=1
  myRow<-ToDoList[i,]
  pheno<-myRow$pheno
  message("Working on trait ",pheno)
  
  # load step10 data
  loaded1<-load(myRow$statistic_path)
  stopifnot(loaded1=="erg1")
  erg2<-copy(erg1)
  dim(erg2)
  message("                 ",pheno," - finished loading (dim: ",dim(erg2)[1]," SNPs, ",dim(erg2)[2]," columns)")
  
  # filter columns
  myNames = names(erg2)[c(2:4,17,18,12,15,8,7,20:22,11,32,33)]
  colsOut<-setdiff(colnames(erg2),myNames)
  erg2[,get("colsOut"):=NULL]
  setcolorder(erg2,myNames)
  dim(erg2)
  
  # rename columns
  setnames(erg2,"pos","bp_hg19")
  setnames(erg2,"effect_allele","EA")
  setnames(erg2,"other_allele","OA")
  setnames(erg2,"eaf","EAF")
  setnames(erg2,"beta_score","beta")
  setnames(erg2,"se_score","SE")
  setnames(erg2,"p_score","pval")
  setnames(erg2,"n","nSamples")
  setnames(erg2,"numberStudies","nStudies")
  setnames(erg2,"reason4exclusion_topl","reason4exclusion")
  setnames(erg2,"invalid_assocs_topl","invalidAssocs")
  
  # add phenotype and cytoband
  erg2[, phenotype:=pheno]
  message("                 ",pheno," - finished filtering and renaming columns")
  
  # order by chr and pos
  setorder(erg2,chr,bp_hg19)
  head(erg2)
  
  # save
  tag = format(Sys.time(), "%Y-%m-%d")
  tag2 = gsub("2023-","23-",tag)
  tag2 = gsub("-","",tag2)
  outFile = paste0("../data/SumStat_", pheno,"_",tag2,".txt")
  fwrite(erg2, file = outFile, quote = F, sep = "\t", col.names = T, row.names = F, na = NA, dec = ".")
  gzip(outFile,destname = paste0(outFile, ".gz"))
  message("                 ",pheno," - finished saving")
  
  # return: phenotype, NR SNPs with QC ok, min NR Samples, median NR Samples, max NR Samples
  myRow[,SumStats:=outFile]
  myRow[,NR_SNPs := dim(erg2[invalidAssocs==F,])[1]]
  myRow[,NR_Samples_min := min(erg2[invalidAssocs==F,nSamples])]
  myRow[,NR_Samples_med := median(erg2[invalidAssocs==F,nSamples])]
  myRow[,NR_Samples_max := max(erg2[invalidAssocs==F,nSamples])]
  
  myRow
  
}
todo2<-rbindlist(dumTab)
todo2[,c(1,4,8:11)]

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

