#' ---
#' title: "Co-localization Part 7: Run coloc cond against lipids"
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
#' *Co-localization analysis between lipids and cond PCSK9 statistics*: 
#' 
#' 
#'   
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))
source("../helperFunctions/colocFunction_jp.R")

#' # Load data ####
#' ***
load("../temp/05_PCSK9cond.RData")
data_GWAS[,.N,by=pheno]

load("../temp/05_OtherGWASs.RData")
myOtherGWAS[,.N,by=phenotype]

myPhenos = unique(data_GWAS$pheno)
myPhenos

myOtherGWASPhenos = c("CAD","HDL","LDL","logTG","nonHDL","TC")

ToDoList = data.table(trait1 = rep(myPhenos,each=6),
                      trait2 = rep(myOtherGWASPhenos,8))
ToDoList[grepl("_female",trait1),trait2 := paste0(trait2,"_INV_FEMALE")]
ToDoList[grepl("_male",trait1),trait2 := paste0(trait2,"_INV_MALE")]
ToDoList[!grepl("INV",trait2),trait2 := paste0(trait2,"_INV_ALL")]
ToDoList[grepl("CAD",trait2),trait2 := "CAD"]

table(is.element(ToDoList$trait2,myOtherGWAS$phenotype))
mySNPs = unique(gsub(".*__","",myPhenos))

#' # Prep data ####
#' ***
data_GWAS[,MAF := freq]
data_GWAS[,chrPos_b37 := paste0("chr1:",bp)]

myOtherGWAS = myOtherGWAS[candidateGene == "PCSK9"]
myOtherGWAS[,table(phenotype)]
setnames(myOtherGWAS,"bp_hg19","bp")
setnames(myOtherGWAS,"beta","bC")
setnames(myOtherGWAS,"SE","bC_se")
setnames(myOtherGWAS,"pval","pC")
setnames(myOtherGWAS,"EA","refA")
myOtherGWAS[,refA2 := refA]
setnames(myOtherGWAS,"nSamples","N")


#' # Run Coloc ####
#' ***
data_GWAS[,ChrPos := paste(Chr,bp,sep=":")]
myOtherGWAS[,ChrPos := paste(chr,bp,sep=":")]

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = copy(ToDoList)
  myRow = myRow[i,]
  
  data_GWAS1 = copy(data_GWAS)
  data_GWAS1 = data_GWAS1[pheno == myRow$trait1 ,]
  
  data_GWAS2 = copy(myOtherGWAS)
  data_GWAS2 = data_GWAS2[phenotype == myRow$trait2 ,]
  
  data_GWAS1 = data_GWAS1[ChrPos %in% data_GWAS2$ChrPos,]
  data_GWAS2 = data_GWAS2[ChrPos %in% data_GWAS1$ChrPos,]
  
  setorder(data_GWAS1,bp)
  setorder(data_GWAS2,bp)
  
  res = colocFunction_jp(tab1 = data_GWAS1,tab2 = data_GWAS2,
                         trait1 = myRow$trait1,trait2 = myRow$trait2,
                         locus = "1p32.3",locus_name = "PCSK9",plotting = F,
                         col_SNPID = "ChrPos", col_pos = "bp",
                         col_beta = "bC",col_se = "bC_se",col_P = "pC",
                         col_N = "N",col_MAF = "MAF",
                         col_effect = "refA",col_other="refA2")
  x2<-as.data.table(res)
  x3<-t(x2)
  x4<-as.data.table(x3)
  names(x4)<-names(res)
  x4[,gene:= "PCSK9"]
  x4[,trait1:= myRow$trait1]
  x4[,trait2:= myRow$trait2]
  x4
  
}
ColocTable = rbindlist(dumTab1)
ColocTable[,SNP := gsub(".*__","",trait1)]
ColocTable[,trait1 := gsub("__.*","",trait1)]
ColocTable[,trait3 := gsub("_INV.*","",trait2)]
ColocTable[,table(PP.H4.abf>=0.75,trait3,SNP)]
ColocTable[,table(PP.H3.abf>=0.75)]
ColocTable[,table(PP.H2.abf>=0.75)]
ColocTable[,table(PP.H1.abf>=0.75)]
ColocTable[,table(PP.H0.abf>=0.75)]

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H3.abf>=0.75,]
ColocTable[PP.H1.abf>=0.75,]

#' # Save results ####
#' ***
ColocTable = ColocTable[,c(10,8,9,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Independent variant", 
                                         "Tested PCSK9 subgroup ",
                                         "Tested other GWAS trait",
                                         "Number of SNPs included in co-localization analysis per test",
                                         "Posterior probability for hypothesis 0: neither trait associated",
                                         "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                         "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                         "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                         "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/05_7_coloc_otherGWAScond.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/05_7_coloc_otherGWAScond.xlsx"

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

