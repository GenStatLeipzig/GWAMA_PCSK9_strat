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
load("../results/05_GCTA_COJO.RData")

ToDoList = copy(IndepSignals)
ToDoList = ToDoList[candidateGene == "PCSK9",]
ToDoList

#' Okay, how many tests will this create?
#' 
#' 1) combined: males_1 vs. females_1
#' 2) combined: males_1 vs. females_2
#' 3) combined: males_2 vs. females_1
#' 4) combined: males_2 vs. females_2
#' 7) free: males_1 vs. females_1
#' 8) free: males_1 vs. females_2
#' 9) free: males_2 vs. females_1
#' 10) free: males_2 vs. females_2
#' 11) combined: free_1 vs. treated
#' 12) combined: free_2 vs. treated
#' 13) males: free_1 vs. treated
#' 14) males: free_2 vs. treated
#' 
#' # Load and filter data ####
#' ***
load("../temp/06_GWAS_allLoci.RData")
data_GWAS_uncond = copy(data_GWAS)
data_GWAS_uncond = data_GWAS_uncond[candidateGene == "PCSK9",]
load("../temp/06_PCSK9cond.RData")
data_GWAS_cond = copy(data_GWAS)

names(data_GWAS_uncond)
names(data_GWAS_cond)
data_GWAS_uncond = data_GWAS_uncond[,c(2,1,3,4,6,10:12,8,6,10:12,16)]
names(data_GWAS_uncond)
data_GWAS_uncond[,MAF := EAF]
data_GWAS_uncond[EAF>0.5,MAF := 1-EAF]
data_GWAS_uncond[,N := ceiling(nSamples)]
data_GWAS_uncond[,refA2 := EA]
names(data_GWAS_uncond) = names(data_GWAS_cond)

data_GWAS = rbind(data_GWAS_cond, data_GWAS_uncond)


#' # Run Coloc ####
#' ***
myPhenos = unique(data_GWAS$pheno)
myPhenos
ToDoList2 = data.table(fix = c(rep("combined",4),rep("free",4),rep("combined",2),rep("male",2)),
                       trait1 = c(myPhenos[7],myPhenos[7],myPhenos[8],myPhenos[8],
                                  myPhenos[9],myPhenos[9],myPhenos[10],myPhenos[10],
                                  myPhenos[14],myPhenos[14],
                                  myPhenos[18],myPhenos[18]),
                       trait2 = c(myPhenos[1],myPhenos[2],myPhenos[1],myPhenos[2],
                                  myPhenos[3],myPhenos[4],myPhenos[3],myPhenos[4],
                                  myPhenos[5],myPhenos[6],
                                  myPhenos[9],myPhenos[10]))
ToDoList2

dumTab1 = foreach(i=1:dim(ToDoList2)[1])%do%{
  #i=1
  myRow = copy(ToDoList2)
  myRow = myRow[i,]
  
  data_GWAS_trait1 = copy(data_GWAS)
  data_GWAS_trait1 = data_GWAS_trait1[pheno == myRow$trait1,]
  data_GWAS_trait2 = copy(data_GWAS)
  data_GWAS_trait2 = data_GWAS_trait2[pheno == myRow$trait2,]
  
  # filter for same SNPs
  data_GWAS_trait1 = data_GWAS_trait1[SNP %in% data_GWAS_trait2$SNP,]
  data_GWAS_trait2 = data_GWAS_trait2[SNP %in% data_GWAS_trait1$SNP,]
  setorder(data_GWAS_trait1,bp)
  setorder(data_GWAS_trait2,bp)
  stopifnot(data_GWAS_trait1$bp == data_GWAS_trait2$bp)
  stopifnot(data_GWAS_trait1$SNP == data_GWAS_trait2$SNP)
  
  trait1 = myRow$trait1
  trait2 = myRow$trait2
  
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
  x4[,trait1:=trait1]
  x4[,trait2:=trait2]
  x4[,fixed:=myRow$fix]
  x4[,NR :=i]
  x4
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
                         description = c("NR of test","cytoband","gene","fixed strata",
                                         "tested conditional strata 1", "tested conditional strata 2", 
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/06_5_coloc_withinPCSK9_cond.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/06_5_coloc_withinPCSK9_cond.xlsx"

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

