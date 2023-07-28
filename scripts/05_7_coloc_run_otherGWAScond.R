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

#' # Get ToDoList ####
#' ***
load("../results/05_GCTA_COJO.RData")

ToDoList = copy(IndepSignals)
ToDoList = ToDoList[multipleSignals == T,]
ToDoList[,otherGWAS := "lipids"]
ToDoList2 = copy(ToDoList)
ToDoList2[,otherGWAS := "CAD"]
ToDoList = rbind(ToDoList2,ToDoList)

#' Okay, how many tests will this create? 
#' 
#' --> for each of the 10 signals all 5 lipid traits (sex-matched)
#' --> for each of the 10 signal CAD
#' 
#' # Load PCSK9 data ####
#' ***
load("../temp/06_PCSK9cond.RData")
data_GWAS[,MAF := freq]
data_GWAS[,chrPos_b37 := paste0("chr1:",bp)]

#' # Load lipid data ####
#' ***
load("../temp/06_OtherGWASs.RData")
myOtherGWAS
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
myOtherGWAS[,phenotype := gsub("MALE","MALES",phenotype)]

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=11
  myRow = copy(ToDoList)
  myRow = myRow[i,]
  sex = unlist(strsplit(myRow$pheno,"_"))[2]
  sex = toupper(sex)
  if(sex %in% c("FREE","TREATED"))sex = "ALL"
  myPheno = paste(myRow$pheno,myRow$rsID,sep="_")
  
  message("working on phenotype ",myPheno, " and candidate gene ",myRow$candidateGene)
  
  data_GWAS1 = copy(data_GWAS)
  data_GWAS1 = data_GWAS1[pheno == myPheno ,]
  
  data_GWAS2 = copy(myOtherGWAS)
  
  if(myRow$otherGWAS == "lipids"){
    myLipidTraits = unique(data_GWAS2$phenotype)
    myLipidTraits = myLipidTraits[grepl(paste0("_",sex),myLipidTraits)]
    
    dumTab2 = foreach(j=1:length(myLipidTraits))%do%{
      #j=1
      myLipidTrait = myLipidTraits[j]
      message("    working on lipid ",myLipidTrait)
      
      data_GWAS3 = copy(data_GWAS2)
      data_GWAS3 = data_GWAS3[phenotype == myLipidTrait,]
      
      data_GWAS4 = copy(data_GWAS1)
      data_GWAS4 = data_GWAS4[ChrPos %in% data_GWAS3$ChrPos,]
      data_GWAS3 = data_GWAS3[ChrPos %in% data_GWAS4$ChrPos,]
      
      setorder(data_GWAS3,bp)
      setorder(data_GWAS4,bp)
      
      res = colocFunction_jp(tab1 = data_GWAS4,tab2 = data_GWAS3,
                             trait1 = myPheno,trait2 = myLipidTrait,
                             locus = "1p32.3",locus_name = "PCSK9",plotting = F,
                             col_SNPID = "ChrPos", col_pos = "bp",
                             col_beta = "bC",col_se = "bC_se",col_P = "pC",
                             col_N = "N",col_MAF = "MAF",
                             col_effect = "refA",col_other="refA2")
      x2<-as.data.table(res)
      x3<-t(x2)
      x4<-as.data.table(x3)
      names(x4)<-names(res)
      x4[,gene:= myRow$candidateGene]
      x4[,trait1:= myPheno]
      x4[,trait2:=myLipidTrait]
      x4
    }
    dumTab2 = rbindlist(dumTab2)
  }else{
    message("    working on trait ",myRow$otherGWAS)
    
    data_GWAS3 = copy(data_GWAS2)
    data_GWAS3 = data_GWAS3[phenotype == myRow$otherGWAS]
    trait2 = unique(data_GWAS3$phenotype)
    
    data_GWAS4 = copy(data_GWAS1)
    data_GWAS4 = data_GWAS4[ChrPos %in% data_GWAS3$ChrPos,]
    data_GWAS3 = data_GWAS3[ChrPos %in% data_GWAS4$ChrPos,]
    
    setorder(data_GWAS3,bp)
    setorder(data_GWAS4,bp)
    
    if(sum(!is.na(data_GWAS3$MAF))==0){
      data_GWAS3[,MAF := NULL]
    }
    
    res = colocFunction_jp(tab1 = data_GWAS4,tab2 = data_GWAS3,
                           trait1 = myPheno,trait2 = myLipidTrait,
                           locus = "1p32.3",locus_name = "PCSK9",plotting = F,
                           col_SNPID = "ChrPos", col_pos = "bp",
                           col_beta = "bC",col_se = "bC_se",col_P = "pC",
                           col_N = "N",col_MAF = "MAF",
                           col_effect = "refA",col_other="refA2")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myRow$candidateGene]
    x4[,trait1:= myPheno]
    x4[,trait2:=trait2]
    dumTab2 = x4
  }
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
ColocTable = ColocTable[,c(7:9,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Candidate Gene of tested region", 
                                         "Tested PCSK9 subgroup ",
                                         "Tested lipid ",
                                         "Number of SNPs included in co-localization analysis per test",
                                         "Posterior probability for hypothesis 0: neither trait associated",
                                         "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                         "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                         "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                         "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/06_7_coloc_otherGWAScond.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/06_7_coloc_otherGWAScond.xlsx"

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

