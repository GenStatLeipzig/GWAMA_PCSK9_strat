#' ---
#' title: "Interaction Tests"
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
#' Perform series of interaction tests: 
#' 
#'  1. 3-way interaction test: SNP x sex x statin
#'  2. 2-way interaction test: SNP x sex
#'  3. 2-way interaction test: SNP x statin
#' 
#' Use only the 14 independent SNPs as identified by GCTA COJO. I will not correct for spurious correlation between the strata. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/TwoWayInteractionTest_jp.R")
source("../helperFunctions/ThreeWayInteractionTest_jp.R")

#' # Load data ####
#' ***
#' Load data and filter for suggestive significant and valid SNPs. Finally, merge all significant SNPs into one object. 
#' 
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]

load("../results/03_GCTA_COJO_filtered.RData")
IndepSignals = IndepSignals_filtered[!duplicated(SNP),]

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Working on data ",myRow$pheno)
  data = fread(myRow$statistic_path)
  data = data[markername %in% IndepSignals$SNP,]
  data
}

result.0 = rbindlist(dumTab)
matched = match(result.0$markername,IndepSignals$SNP)
result.0[,candidateGene := IndepSignals[matched,candidateGene]]
result.0[,rsID := gsub(":.*","",markername)]
setorder(result.0,chr,bp_hg19)

#' # Get Strata info per SNP ####
#' ***
#' 
result.1 = copy(result.0)
setorder(result.1,pval)
result.1 = result.1[!duplicated(markername)]
matched = match(result.0$markername,result.1$markername)
result.0[,bestPheno := result.1[matched,phenotype]]

result.2 = copy(result.0)
result.2[,table(phenotype,bestPheno)]

result.2[, bestSex := "combined"]
result.2[grepl("_male",bestPheno), bestSex := "males"]
result.2[grepl("_female",bestPheno), bestSex := "females"]
result.2[,table(bestPheno,bestSex)]

result.2[, bestStatin := "combined"]
result.2[grepl("_free",bestPheno), bestStatin := "free"]
result.2[grepl("_treated",bestPheno), bestStatin := "treated"]
result.2[,table(bestPheno,bestStatin)]

result.2[, sex := "combined"]
result.2[grepl("_male",phenotype), sex := "males"]
result.2[grepl("_female",phenotype), sex := "females"]
result.2[,table(bestPheno,sex)]

result.2[, statin := "combined"]
result.2[grepl("_free",phenotype), statin := "free"]
result.2[grepl("_treated",phenotype), statin := "treated"]
result.2[,table(bestPheno,statin)]

result.2[,corSex := 0]
result.2[,corStatin := 0]

#' # Check EAFs ####
#' ***
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_females_free",EAF])
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_females_treated",EAF])
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_males",EAF])
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_males_free",EAF])
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_males_treated",EAF])
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_free",EAF])
plot(result.2[phenotype == "PCSK9_females",EAF],result.2[phenotype == "PCSK9_treated",EAF])

save(result.2,file="../temp/04_IATest_input.RData")

#' # Interaction Tests ####
#' ***
#' ## Run 3-way test ####
IATab_3way = ThreeWayInteractionTest_jp(data = result.2,
                                  pheno1 = "PCSK9_males_treated",
                                  pheno2 = "PCSK9_females_treated",
                                  pheno3 = "PCSK9_males_free",
                                  pheno4 = "PCSK9_females_free")

matched = match(IATab_3way$markername,result.2$markername)
table(is.na(matched))
table(IATab_3way$markername == result.2$markername[matched])
IATab_3way[,candidateGene := result.2$candidateGene[matched]]

#' ## Run 2-way test for sex ####
IATab_2way_sex = TwoWayInteractionTest_jp(data = result.2,
                                          pheno1 = "males",
                                          pheno2 = "females",
                                          type = "sexIA",
                                          corCol = "corSex")

matched = match(IATab_2way_sex$markername,result.2$markername)
table(is.na(matched))
table(IATab_2way_sex$markername == result.2$markername[matched])
IATab_2way_sex[,candidateGene := result.2$candidateGene[matched]]

#' ## Run 2-way test for statin ####
IATab_2way_statin = TwoWayInteractionTest_jp(data = result.2,
                                             pheno1 = "treated",
                                             pheno2 = "free",
                                             type = "statinIA",
                                             corCol = "corStatin")

matched = match(IATab_2way_statin$markername,result.2$markername)
table(is.na(matched))
table(IATab_2way_statin$markername == result.2$markername[matched])
IATab_2way_statin[,candidateGene := result.2$candidateGene[matched]]

#' ## Combine data sets ####
matched = match(IATab_3way$markername,result.1$markername)
IATab_3way[,bestPheno := result.1[matched,phenotype]]
IATab_3way[,type := "3way"]
IATab_3way[,fix := NA]

IATab_2way_sex[,cor := NULL]
IATab_2way_statin[,cor := NULL]

IATab = rbind(IATab_3way,IATab_2way_sex,IATab_2way_statin,fill = T)
IATab = IATab[,c(1:5,34:37,30:33,6:29)]
setorder(IATab,chr,bp_hg19)
IATab[,c(1:13)]

#' ## Correct for multiple testing ####
myFDR = addHierarchFDR(pvalues = IATab$IA_pval, categs = IATab$markername,quiet = F)
IATab[,IA_pval_adj := myFDR$fdr_level1]
IATab[,IA_hierarch_fdr5proz := myFDR$hierarch_fdr5proz]
IATab[,table(IA_hierarch_fdr5proz,candidateGene)]

#' ## Summary ####
IATab[IA_hierarch_fdr5proz==T & type=="3way",]

#' **Summary** 3-way interaction:
#' 
#' - **PCSK9, rs11583680**: significant 3-way interaction, no association in statin-treated females
#' - **ALOX5**: significant 3-way interaction, only significant association in statin-free males
#'        
IATab[IA_hierarch_fdr5proz==T & type=="sexIA",c(1:25,38,39)]

#' **Summary** 2-way sex interaction:
#' 
#' - **PCSK9, rs693668**: sex-related, stronger in males
#' - **PRKAG2**: male-specific association
#' - **ALOX5**: male-specific association
#' - **SLCO1B1**: female-specfic association
#' - **NOS1**: male-specific association
#' 
IATab[IA_hierarch_fdr5proz==T & type=="statinIA",c(1:25,38,39)]

#' **Summary** 2-way statin interaction:
#' 
#' - **PCSK9, rs11591147**: sex-related, stronger in statin-free individuals
#' - **PCSK9, rs11583680**: sex-related, stronger in statin-treated individuals
#' - **APOB**: free-specific effect
#' - **KHDRBS2**: treated-specific effect
#' - **PRKAG2**: free-specific effect
#' - **ALOX5**: free-specific effect
#' 
#' ## Save ####
IATab = IATab[,c(1:13,38,39,14:37)]
save(IATab,file="../results/04_IATest_complete.RData")

IATab_filtered = copy(IATab)
IATab_filtered = IATab_filtered[,c(1:15)]
save(IATab_filtered,file="../results/04_IATest_filtered.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

