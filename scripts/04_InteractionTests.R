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
#' Use only the 4 independent SNPs PCSK9 as identified by GCTA COJO and the 10 lead SNPs from the GWAS. For the chr 12 locus, I choose the SNP with an rsID (easier to reference). I will not correct for spurious correlation between the strata. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/TwoWayInteractionTest_jp.R")

#' # Load data ####
#' ***
#' Load data and filter for suggestive significant and valid SNPs. Finally, merge all significant SNPs into one object. 
#' 
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]

PCSK9_SNPs = fread("../results/03_GCTA_COJO_joint/PCSK9.snplist",header=F)
load("../results/02_LociOverallPhenotypes_filtered.RData")
Other_SNPs = result.5[-1,markername]
Other_SNPs[grep("12:21117156:GAA:GA",Other_SNPs)] = "rs4762806:21067768:T:C"
mySNPs = c(PCSK9_SNPs$V1,Other_SNPs)

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Working on data ",myRow$pheno)
  data = fread(myRow$statistic_path)
  data = data[markername %in% mySNPs,]
  data
}

result.0 = rbindlist(dumTab)
result.5[grep("12:21117156:GAA:GA",markername),markername := "rs4762806:21067768:T:C"]
matched = match(result.0$markername,result.5$markername)
result.0[,candidateGene := result.5[matched,candidateGene]]
result.0[is.na(candidateGene),candidateGene := "PCSK9"]
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
#' ## Run 2-way test for sex ####
IATab_2way_sex = TwoWayInteractionTest_jp(data = result.2,
                                          pheno1 = "males",
                                          pheno2 = "females",
                                          type = "sexIA",
                                          useBestPheno = T,
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
                                             useBestPheno = T,
                                             corCol = "corStatin")

matched = match(IATab_2way_statin$markername,result.2$markername)
table(is.na(matched))
table(IATab_2way_statin$markername == result.2$markername[matched])
IATab_2way_statin[,candidateGene := result.2$candidateGene[matched]]

#' ## Combine data sets ####
IATab_2way_sex[,cor := NULL]
IATab_2way_statin[,cor := NULL]

IATab = rbind(IATab_2way_sex,IATab_2way_statin,fill = T)
IATab = IATab[,c(1:9,22:25,10:21)]
setorder(IATab,chr,bp_hg19)
IATab[,c(1:13)]

#' ## Correct for multiple testing ####
myFDR = addHierarchFDR(pvalues = IATab$IA_pval, categs = IATab$markername,quiet = F)
IATab[,IA_pval_adj := myFDR$fdr_level1]
IATab[,IA_hierarch_fdr5proz := myFDR$hierarch_fdr5proz]
IATab[,table(IA_hierarch_fdr5proz,candidateGene)]

IATab[,diff2 := abs(trait1_beta) - abs(trait2_beta)]
IATab[,diff3 := sign(diff2)]

IATab[,type2 := "unspecific"]

IATab[IA_hierarch_fdr5proz==T & type=="sexIA" & diff3==1 & trait2_pval>1e-6, type2 := "male-specific"]
IATab[IA_hierarch_fdr5proz==T & type=="sexIA" & diff3==-1 & trait1_pval>1e-6, type2 := "female-specific"]
IATab[IA_hierarch_fdr5proz==T & type=="sexIA" & diff3==1 & trait1_pval<=1e-6 & trait2_pval<=1e-6, type2 := "sex-related (stronger in males)"]
IATab[IA_hierarch_fdr5proz==T & type=="sexIA" & diff3==-1 & trait1_pval<=1e-6 & trait2_pval<=1e-6, type2 := "sex-related (stronger in females)"]

IATab[IA_hierarch_fdr5proz==T & type=="statinIA" & diff3==1 & trait2_pval>1e-6, type2 := "treated-specific"]
IATab[IA_hierarch_fdr5proz==T & type=="statinIA" & diff3==-1 & trait1_pval>1e-6, type2 := "free-specific"]
IATab[IA_hierarch_fdr5proz==T & type=="statinIA" & diff3==1 & trait1_pval<=1e-6 & trait2_pval<=1e-6, type2 := "statin-related (stronger in treated)"]
IATab[IA_hierarch_fdr5proz==T & type=="statinIA" & diff3==-1 & trait1_pval<=1e-6 & trait2_pval<=1e-6, type2 := "statin-related (stronger in free)"]

IATab[IA_hierarch_fdr5proz==T,.N,c("type2","diff3")]

#' ## Summary ####
IATab[IA_hierarch_fdr5proz==T & type=="sexIA",c(1:13,26,27,30)]

#' **Summary** 2-way sex interaction:
#' 
#' - **PCSK9, rs693668**: sex-related, stronger in males
#' - **PLPP3**: male-specific association
#' - **PRKAG2**: male-specific association
#' - **ALOX5**: male-specific association
#' - **SLCO1B1**: female-specfic association
#' - **NOS1**: male-specific association
#' 
IATab[IA_hierarch_fdr5proz==T & type=="statinIA",c(1:13,26,27,30)]

#' **Summary** 2-way statin interaction:
#' 
#' - **PCSK9, rs11591147**: statin-related, stronger in statin-free individuals
#' - **PCSK9, rs11583680**: statin-related, stronger in statin-treated individuals
#' - **APOB**: free-specific effect
#' - **KHDRBS2**: treated-specific effect
#' - **PRKAG2**: free-specific effect
#' - **ALOX**: free-specific effect
#' - **JMJD1C**: free-specific effect
#' 
#' ## Save ####
IATab = IATab[,c(1:13,26,27,30,14:25)]
save(IATab,file="../results/04_IATest_2way_complete.RData")

IATab_filtered = copy(IATab)
IATab_filtered = IATab_filtered[,c(1:16)]
save(IATab_filtered,file="../results/04_IATest_2way_filtered.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

