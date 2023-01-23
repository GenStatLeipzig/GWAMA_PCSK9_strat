#' ---
#' title: "Co-localization Part 3: Run coloc cond against eQTLs"
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
#' *Co-localization analysis between eQTLs and cond PCSK9 statistics*: 
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

#' # Get ToDoList ####
#' ***
load("../results/05_GCTA_COJO.RData")

ToDoList = copy(IndepSignals)
ToDoList = ToDoList[multipleSignals == T,]
ToDoList

load("../results/06_a_usedGenes.RData")
myGenTab = myGenTab[genename == "PCSK9"]
filt = !is.na(myGenTab)
myGenTab = myGenTab[,filt,with=F]
myTissues = names(myGenTab)[8:48]

#' Okay, how many tests will this create? 
#' 
#' --> for each of the 10 signals all tissues with PCSK9 eQTLs available (41 tissues)!
#' 
#' # Load and filter data ####
#' ***
#' I need 10 data sets 
#' 

dumTab = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  myRow
  
  data = fread(myRow$input_CS)
  data[,pheno := paste(myRow$pheno,myRow$rsID,sep="_")]
  data = data[!is.na(bC),]
  data
}
data_GWAS = rbindlist(dumTab)

data_GWAS[,MAF := freq]
data_GWAS[freq>0.5,MAF := 1-freq]
data_GWAS[,N := ceiling(n)]
data_GWAS[,refA2 := refA]
save(data_GWAS,file = "../temp/06_PCSK9cond.RData")
data_GWAS[,chrPos_b37 := paste0("chr1:",bp)]

#' # Run Coloc ####
#' ***
myPhenos = unique(data_GWAS$pheno)
myPhenos

dumTab1 = foreach(i=1:length(myTissues))%do%{
  #i=1
  myTissue = myTissues[i]
  message("Working on tissue ",myTissue, ", ",i, " of 41")
  
  load(paste0("../temp/06_coloc/GTEx_v8_filtered_",myTissue,".RData"))
  data0 = data0[gene == "PCSK9"]
  setnames(data0,"pos_b37","bp")
  setnames(data0,"beta","bC")
  setnames(data0,"se","bC_se")
  setnames(data0,"pval","pC")
  setnames(data0,"effect_allele","refA")
  data0[,refA2 := refA]
  setnames(data0,"maf","MAF")
  setnames(data0,"n_samples","N")
  
  dumTab2 = foreach(j=1:length(myPhenos))%do%{
    #j=1
    myPheno = myPhenos[j]
    message("     Working on phenotype ",myPheno)
    
    data_GWAS1 = copy(data_GWAS)
    data_GWAS1 = data_GWAS1[pheno == myPheno,]
    data_GWAS1 = data_GWAS1[!is.na(bC),]
    
    data_eQTL1 = copy(data0)
    data_GWAS1 = data_GWAS1[chrPos_b37 %in% data_eQTL1$chrPos_b37,]
    data_eQTL1 = data_eQTL1[chrPos_b37 %in% data_GWAS1$chrPos_b37,]
    
    setorder(data_GWAS1,bp)
    setorder(data_eQTL1,bp)
    stopifnot(data_GWAS1$bp == data_eQTL1$bp)
    stopifnot(data_GWAS1$chrPos_b37 == data_eQTL1$chrPos_b37)
    
    trait1 = myPheno
    trait2 = myTissue
    
    res = colocFunction_jp(tab1 = data_GWAS1,tab2 = data_eQTL1,
                           trait1 = trait1,trait2 = trait2,
                           locus = "1p32.3",locus_name = "PCSK9",plotting = F,
                           col_SNPID = "chrPos_b37", col_pos = "bp",
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
ColocTable = ColocTable[,c(7,9,8,10,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("cytoband","phenotype setting and condition","gene","tissue",
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (PCSK9 trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/06_3_coloc_eQTLs_cond.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/06_3_coloc_eQTLs_cond.xlsx"

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

