#' ---
#' title: "Co-localization Part 4: Run coloc within PCSK9 data"
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
#' Used comparisons (pairwise non-overlapping samples):
#' 
#' * combined: males vs. females
#' * free: males vs. females
#' * treated: males vs. females
#' * combined: free vs. treated
#' * females: free vs. treated
#' * males: free vs. treated
#' 
#' For locus *PCSK9*, I will test all 6 possible comparisons. 
#'  
#' For the other loci, I will test 2 colocs, one sex-specific, one statin-specific. E.g. best phenotype is PCSK9_male_free, then I test free: males vs. females and males: free vs. treated
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
ToDoList = ToDoList[c(1,13:24),c(1:5,15,19:21,26)]
ToDoList

#' # Load and filter data ####
#' ***
myPhenos = unique(ToDoList$pheno)

dumTab1 = foreach(i = 1:length(myPhenos))%do%{
  # i=1
  myPheno = myPhenos[i]
  message("Working on phenotype ",myPheno)
  
  myFileName = paste0("../data/SumStat_",myPheno,"_230120.txt.gz")
  
  data_GWAS1 = fread(myFileName)
  data_GWAS1 = data_GWAS1[invalidAssocs==F,]
  
  dumTab2 = foreach(j = 1:dim(ToDoList)[1])%do%{
    # j=1
    myRow = copy(ToDoList)
    myRow = myRow[j,]
    message("     Working on gene ",myRow$candidateGene)
    
    data_GWAS2 = copy(data_GWAS1)
    data_GWAS2 = data_GWAS2[chr == myRow$Chr,]
    data_GWAS2 = data_GWAS2[bp_hg19 >= myRow$region_start,]
    data_GWAS2 = data_GWAS2[bp_hg19 <= myRow$region_end,]
    data_GWAS2[,candidateGene := myRow$candidateGene]
    
    data_GWAS2[,setting_sex := gsub("PCSK9_","",phenotype)]
    data_GWAS2[,setting_statin := setting_sex]
    
    data_GWAS2[,setting_sex := gsub("_.*","",setting_sex)]
    data_GWAS2[,setting_statin := gsub(".*_","",setting_statin)]
    data_GWAS2
    
  }
  dumTab2 = rbindlist(dumTab2)
  stopifnot(length(unique(dumTab2$candidateGene))==13)
  dumTab2
}
data_GWAS = rbindlist(dumTab1)
data_GWAS[setting_statin == "females", setting_statin := "combined"]
data_GWAS[setting_statin == "males", setting_statin := "combined"]
data_GWAS[setting_sex == "free", setting_sex := "combined"]
data_GWAS[setting_sex == "treated", setting_sex := "combined"]
data_GWAS[,table(phenotype, setting_sex)]
data_GWAS[,table(phenotype, setting_statin)]

data_GWAS[,MAF := EAF]
data_GWAS[EAF>0.5,MAF := 1-EAF]
save(data_GWAS,file = "../temp/06_GWAS_allLoci.RData")


#' # Run Coloc ####
#' ***
dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=2
  myRow = copy(ToDoList)
  myRow = myRow[i,]
  
  data_GWAS1 = copy(data_GWAS)
  data_GWAS1 = data_GWAS1[candidateGene == myRow$candidateGene,]
  
  if(i!=1){
    myBestPheno = myRow$pheno
    myBestSex = data_GWAS1[phenotype == myBestPheno,unique(setting_sex)]
    myBestStatin = data_GWAS1[phenotype == myBestPheno,unique(setting_statin)]
    
    # Check 1: males vs females in same statin setting
    data_GWAS2 = copy(data_GWAS1)
    data_GWAS2 = data_GWAS2[setting_statin == myBestStatin,]
    data_GWAS2_males = copy(data_GWAS2)
    data_GWAS2_males = data_GWAS2_males[setting_sex == "males"]
    data_GWAS2_females = copy(data_GWAS2)
    data_GWAS2_females = data_GWAS2_females[setting_sex == "females"]
    
    data_GWAS2_males = data_GWAS2_males[markername %in% data_GWAS2_females$markername,]
    data_GWAS2_females = data_GWAS2_females[markername %in% data_GWAS2_males$markername,]
   
    setorder(data_GWAS2_males,bp_hg19)
    setorder(data_GWAS2_females,bp_hg19)
    stopifnot(data_GWAS2_males$bp_hg19 == data_GWAS2_females$bp_hg19)
    stopifnot(data_GWAS2_males$markername == data_GWAS2_females$markername)
   
    trait1 = paste0("PCSK9_males_",myBestStatin)
    trait2 = paste0("PCSK9_females_",myBestStatin)
    res = colocFunction_jp(tab1 = data_GWAS2_males,tab2 = data_GWAS2_females,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myRow$candidateGene]
    x4[,trait1:=trait1]
    x4[,trait2:=trait2]
    
    # Check 2: free vs treated in same sex setting
    data_GWAS3 = copy(data_GWAS1)
    data_GWAS3 = data_GWAS3[setting_sex == myBestSex,]
    data_GWAS3_free = copy(data_GWAS3)
    data_GWAS3_free = data_GWAS3_free[setting_statin == "free"]
    data_GWAS3_treated = copy(data_GWAS3)
    data_GWAS3_treated = data_GWAS3_treated[setting_statin == "treated"]
    
    data_GWAS3_free = data_GWAS3_free[markername %in% data_GWAS3_treated$markername,]
    data_GWAS3_treated = data_GWAS3_treated[markername %in% data_GWAS3_free$markername,]
    
    setorder(data_GWAS3_free,bp_hg19)
    setorder(data_GWAS3_treated,bp_hg19)
    stopifnot(data_GWAS3_free$bp_hg19 == data_GWAS3_treated$bp_hg19)
    stopifnot(data_GWAS3_free$markername == data_GWAS3_treated$markername)
    
    trait1 = paste0("PCSK9_",myBestSex,"_free")
    trait2 = paste0("PCSK9_",myBestSex,"_treated")
    res = colocFunction_jp(tab1 = data_GWAS3_free,tab2 = data_GWAS3_treated,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x5<-as.data.table(x3)
    names(x5)<-names(res)
    x5[,gene:= myRow$candidateGene]
    x5[,trait1:=trait1]
    x5[,trait2:=trait2]
    
    # Combine results and return
    x6 = rbind(x4,x5)
    x6
  }else{
    # Check 1: males vs females in combined setting
    data_GWAS2 = copy(data_GWAS1)
    data_GWAS2 = data_GWAS2[setting_statin == "combined",]
    data_GWAS2_males = copy(data_GWAS2)
    data_GWAS2_males = data_GWAS2_males[setting_sex == "males"]
    data_GWAS2_females = copy(data_GWAS2)
    data_GWAS2_females = data_GWAS2_females[setting_sex == "females"]
    
    data_GWAS2_males = data_GWAS2_males[markername %in% data_GWAS2_females$markername,]
    data_GWAS2_females = data_GWAS2_females[markername %in% data_GWAS2_males$markername,]
    
    setorder(data_GWAS2_males,bp_hg19)
    setorder(data_GWAS2_females,bp_hg19)
    stopifnot(data_GWAS2_males$bp_hg19 == data_GWAS2_females$bp_hg19)
    stopifnot(data_GWAS2_males$markername == data_GWAS2_females$markername)
    
    trait1 = paste0("PCSK9_males_combined")
    trait2 = paste0("PCSK9_females_combined")
    res = colocFunction_jp(tab1 = data_GWAS2_males,tab2 = data_GWAS2_females,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myRow$candidateGene]
    x4[,trait1:=trait1]
    x4[,trait2:=trait2]
    check1 = copy(x4)
    
    # Check 2: males vs females in free setting
    data_GWAS2 = copy(data_GWAS1)
    data_GWAS2 = data_GWAS2[setting_statin == "free",]
    data_GWAS2_males = copy(data_GWAS2)
    data_GWAS2_males = data_GWAS2_males[setting_sex == "males"]
    data_GWAS2_females = copy(data_GWAS2)
    data_GWAS2_females = data_GWAS2_females[setting_sex == "females"]
    
    data_GWAS2_males = data_GWAS2_males[markername %in% data_GWAS2_females$markername,]
    data_GWAS2_females = data_GWAS2_females[markername %in% data_GWAS2_males$markername,]
    
    setorder(data_GWAS2_males,bp_hg19)
    setorder(data_GWAS2_females,bp_hg19)
    stopifnot(data_GWAS2_males$bp_hg19 == data_GWAS2_females$bp_hg19)
    stopifnot(data_GWAS2_males$markername == data_GWAS2_females$markername)
    
    trait1 = paste0("PCSK9_males_free")
    trait2 = paste0("PCSK9_females_free")
    res = colocFunction_jp(tab1 = data_GWAS2_males,tab2 = data_GWAS2_females,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myRow$candidateGene]
    x4[,trait1:=trait1]
    x4[,trait2:=trait2]
    check2 = copy(x4)
    
    # Check 3: males vs females in treated setting
    data_GWAS2 = copy(data_GWAS1)
    data_GWAS2 = data_GWAS2[setting_statin == "treated",]
    data_GWAS2_males = copy(data_GWAS2)
    data_GWAS2_males = data_GWAS2_males[setting_sex == "males"]
    data_GWAS2_females = copy(data_GWAS2)
    data_GWAS2_females = data_GWAS2_females[setting_sex == "females"]
    
    data_GWAS2_males = data_GWAS2_males[markername %in% data_GWAS2_females$markername,]
    data_GWAS2_females = data_GWAS2_females[markername %in% data_GWAS2_males$markername,]
    
    setorder(data_GWAS2_males,bp_hg19)
    setorder(data_GWAS2_females,bp_hg19)
    stopifnot(data_GWAS2_males$bp_hg19 == data_GWAS2_females$bp_hg19)
    stopifnot(data_GWAS2_males$markername == data_GWAS2_females$markername)
    
    trait1 = paste0("PCSK9_males_treated")
    trait2 = paste0("PCSK9_females_treated")
    res = colocFunction_jp(tab1 = data_GWAS2_males,tab2 = data_GWAS2_females,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myRow$candidateGene]
    x4[,trait1:=trait1]
    x4[,trait2:=trait2]
    check3 = copy(x4)
    
    # Check 4: free vs treated in combined setting
    data_GWAS3 = copy(data_GWAS1)
    data_GWAS3 = data_GWAS3[setting_sex == "combined",]
    data_GWAS3_free = copy(data_GWAS3)
    data_GWAS3_free = data_GWAS3_free[setting_statin == "free"]
    data_GWAS3_treated = copy(data_GWAS3)
    data_GWAS3_treated = data_GWAS3_treated[setting_statin == "treated"]
    
    data_GWAS3_free = data_GWAS3_free[markername %in% data_GWAS3_treated$markername,]
    data_GWAS3_treated = data_GWAS3_treated[markername %in% data_GWAS3_free$markername,]
    
    setorder(data_GWAS3_free,bp_hg19)
    setorder(data_GWAS3_treated,bp_hg19)
    stopifnot(data_GWAS3_free$bp_hg19 == data_GWAS3_treated$bp_hg19)
    stopifnot(data_GWAS3_free$markername == data_GWAS3_treated$markername)
    
    trait1 = paste0("PCSK9_combined_free")
    trait2 = paste0("PCSK9_combined_treated")
    res = colocFunction_jp(tab1 = data_GWAS3_free,tab2 = data_GWAS3_treated,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x5<-as.data.table(x3)
    names(x5)<-names(res)
    x5[,gene:= myRow$candidateGene]
    x5[,trait1:=trait1]
    x5[,trait2:=trait2]
    check4 = copy(x5)
    
    # Check 5: free vs treated in males setting
    data_GWAS3 = copy(data_GWAS1)
    data_GWAS3 = data_GWAS3[setting_sex == "males",]
    data_GWAS3_free = copy(data_GWAS3)
    data_GWAS3_free = data_GWAS3_free[setting_statin == "free"]
    data_GWAS3_treated = copy(data_GWAS3)
    data_GWAS3_treated = data_GWAS3_treated[setting_statin == "treated"]
    
    data_GWAS3_free = data_GWAS3_free[markername %in% data_GWAS3_treated$markername,]
    data_GWAS3_treated = data_GWAS3_treated[markername %in% data_GWAS3_free$markername,]
    
    setorder(data_GWAS3_free,bp_hg19)
    setorder(data_GWAS3_treated,bp_hg19)
    stopifnot(data_GWAS3_free$bp_hg19 == data_GWAS3_treated$bp_hg19)
    stopifnot(data_GWAS3_free$markername == data_GWAS3_treated$markername)
    
    trait1 = paste0("PCSK9_males_free")
    trait2 = paste0("PCSK9_males_treated")
    res = colocFunction_jp(tab1 = data_GWAS3_free,tab2 = data_GWAS3_treated,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x5<-as.data.table(x3)
    names(x5)<-names(res)
    x5[,gene:= myRow$candidateGene]
    x5[,trait1:=trait1]
    x5[,trait2:=trait2]
    check5 = copy(x5)
    
    # Check 6: free vs treated in females setting
    data_GWAS3 = copy(data_GWAS1)
    data_GWAS3 = data_GWAS3[setting_sex == "females",]
    data_GWAS3_free = copy(data_GWAS3)
    data_GWAS3_free = data_GWAS3_free[setting_statin == "free"]
    data_GWAS3_treated = copy(data_GWAS3)
    data_GWAS3_treated = data_GWAS3_treated[setting_statin == "treated"]
    
    data_GWAS3_free = data_GWAS3_free[markername %in% data_GWAS3_treated$markername,]
    data_GWAS3_treated = data_GWAS3_treated[markername %in% data_GWAS3_free$markername,]
    
    setorder(data_GWAS3_free,bp_hg19)
    setorder(data_GWAS3_treated,bp_hg19)
    stopifnot(data_GWAS3_free$bp_hg19 == data_GWAS3_treated$bp_hg19)
    stopifnot(data_GWAS3_free$markername == data_GWAS3_treated$markername)
    
    trait1 = paste0("PCSK9_females_free")
    trait2 = paste0("PCSK9_females_treated")
    res = colocFunction_jp(tab1 = data_GWAS3_free,tab2 = data_GWAS3_treated,
                           trait1 = trait1,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "markername", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x5<-as.data.table(x3)
    names(x5)<-names(res)
    x5[,gene:= myRow$candidateGene]
    x5[,trait1:=trait1]
    x5[,trait2:=trait2]
    check6 = copy(x5)
    
    # Combine results and return
    x6 = rbind(check1,check2,check3,check4,check5,check6)
    x6
  }
 
  ColocTable = x6
  ColocTable
}

ColocTable = rbindlist(dumTab1)

ColocTable[,table(PP.H4.abf>=0.75)]
ColocTable[,table(PP.H3.abf>=0.75)]
ColocTable[,table(PP.H2.abf>=0.75)]
ColocTable[,table(PP.H1.abf>=0.75)]
ColocTable[,table(PP.H0.abf>=0.75)]

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H2.abf>=0.75,]
ColocTable[PP.H1.abf>=0.75,]

#' **Summary**: 
#' 
#' * APOB & PCSK9 have the same signals in males and females (yes, weaker in one sex, but shared signal!)
#' * PLB1, CYP51A1P3, & SLCO1B1 have only signal in females (not in males)
#' * MYNC has signal only in males & only in free setting
#' * PRKAG2 & NOS1 have only signal in males (not in females)
#' * TM6SF2 has signal only in free setting
#' 
#' # Save results ####
#' ***
ColocTable = ColocTable[,c(7:9,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Candidate Gene of tested region", 
                                         "Tested PCSK9 subgroup 1 ",
                                        "Tested PCSK9 subgroup 1 ",
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/6_4_coloc_withinPCSK9.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/6_4_coloc_withinPCSK9.xlsx"

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

