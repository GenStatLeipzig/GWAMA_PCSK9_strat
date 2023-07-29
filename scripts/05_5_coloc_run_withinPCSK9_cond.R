#' ---
#' title: "Co-localization Part 5: Run coloc within PCSK9 data cond"
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
#' *Co-localization analysis between strata*: 
#' 
#' Used comparisons:
#' 
#' * combined: males vs. females
#' * free: males vs. females
#' * combined: free vs. treated
#' * females: free vs. treated
#' * males: free vs. treated
#' 
#'   
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))
source("../helperFunctions/colocFunction_jp.R")

#' # Get ToDoList ####
#' ***
load("../temp/05_PCSK9cond.RData")
data_GWAS[,.N,by=pheno]

#' Okay, how many tests will this create?
#' 
#' 1) SNP 1 - males vs females
#' 2) SNP 1 - free vs treated
#' 3) SNP 1 - males_free vs females_free
#' 4) SNP 1 - males_treated vs females_treated
#' 5) SNP 1 - males_free vs males_treated
#' 6) SNP 1 - females_free vs females_treated
#' 
#' 24 tests (4 SNPs)
#' 
data_GWAS_cond = copy(data_GWAS)

#' # Run Coloc ####
#' ***
myPhenos = unique(data_GWAS$pheno)
myPhenos
ToDoList = data.table(fix = c(rep("combined",2),"free","treated","males","females"),
                      trait1 = c("males","free","males_free","males_treated",
                                 "males_free","females_free"),
                      trait2 = c("females","treated","females_free","females_treated","males_treated","females_treated"))
ToDoList
mySNPs = unique(gsub(".*__","",myPhenos))

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=4
  myRow = copy(ToDoList)
  myRow = myRow[i,]
  
  dumTab2 = foreach(j=1:length(mySNPs))%do%{
    #j=1
    mySNP = mySNPs[j]
    trait1 = paste0("PCSK9_",myRow$trait1,"__",mySNP)
    trait2 = paste0("PCSK9_",myRow$trait2,"__",mySNP)
    
    data_GWAS_trait1 = copy(data_GWAS)
    data_GWAS_trait1 = data_GWAS_trait1[pheno == trait1,]
    data_GWAS_trait2 = copy(data_GWAS)
    data_GWAS_trait2 = data_GWAS_trait2[pheno == trait2,]
    
    # filter for same SNPs
    data_GWAS_trait1 = data_GWAS_trait1[SNP %in% data_GWAS_trait2$SNP,]
    data_GWAS_trait2 = data_GWAS_trait2[SNP %in% data_GWAS_trait1$SNP,]
    setorder(data_GWAS_trait1,bp)
    setorder(data_GWAS_trait2,bp)
    stopifnot(data_GWAS_trait1$bp == data_GWAS_trait2$bp)
    stopifnot(data_GWAS_trait1$SNP == data_GWAS_trait2$SNP)
    
    res = colocFunction_jp(tab1 = data_GWAS_trait1,tab2 = data_GWAS_trait2,
                           trait1 = trait1,trait2 = trait2,
                           locus = "1p32.3",locus_name = "PCSK9",plotting = F,
                           col_SNPID = "SNP", col_pos = "bp",
                           col_beta = "bC",col_se = "bC_se",col_P = "pC",
                           col_N = "N",col_MAF = "MAF",
                           col_effect = "refA",col_other="refA2")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,locus:= "1p32.3"]
    x4[,gene:= "PCSK9"]
    x4[,trait1:=myRow$trait1]
    x4[,trait2:=myRow$trait2]
    x4[,fixed:=myRow$fix]
    x4[,SNP := mySNP]
    x4
  }
  dumTab2 = rbindlist(dumTab2)
  dumTab2
}

ColocTable = rbindlist(dumTab1)

ColocTable[,table(PP.H4.abf>=0.75)]
ColocTable[,table(PP.H3.abf>=0.75)]
ColocTable[,table(PP.H2.abf>=0.75)]
ColocTable[,table(PP.H1.abf>=0.75)]
ColocTable[,table(PP.H0.abf>=0.75)]

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H3.abf>=0.75,]
ColocTable[PP.H1.abf>=0.75,]

#' # Save results ####
#' ***
ColocTable = ColocTable[,c(12,7,8,11,9,10,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Indep SNP","cytoband","gene","fixed strata",
                                         "tested conditional strata 1", "tested conditional strata 2", 
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/05_5_coloc_withinPCSK9_cond.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/05_5_coloc_withinPCSK9_cond.xlsx"

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

