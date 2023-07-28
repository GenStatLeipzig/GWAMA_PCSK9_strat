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
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/TwoWayInteractionTest_jp.R")
source("../helperFunctions/ThreeWayInteractionTest_jp.R")
source("../helperFunctions/getCorrelation.R")

#' # Load data ####
#' ***
#' Load data and filter for suggestive significant and valid SNPs. Finally, merge all significant SNPs into one object. 
#' 
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Working on data ",myRow$pheno)
  data = fread(myRow$statistic_path)
  data
}

result.0 = rbindlist(dumTab)

cor_sex_combined = getCorrelation(data = result.0,pheno1 = "PCSK9_males",pheno2 = "PCSK9_females")
cor_sex_free = getCorrelation(data = result.0,pheno1 = "PCSK9_males_free",pheno2 = "PCSK9_females_free")
cor_sex_treated = getCorrelation(data = result.0,pheno1 = "PCSK9_males_treated",pheno2 = "PCSK9_females_treated")
cor_statin_combined = getCorrelation(data = result.0,pheno1 = "PCSK9_treated",pheno2 = "PCSK9_free")
cor_statin_males = getCorrelation(data = result.0,pheno1 = "PCSK9_males_treated",pheno2 = "PCSK9_males_free")
cor_statin_females = getCorrelation(data = result.0,pheno1 = "PCSK9_females_treated",pheno2 = "PCSK9_females_free")
corTab = rbind(cor_sex_combined,cor_sex_free,cor_sex_treated,
               cor_statin_combined,cor_statin_males,cor_statin_females)
corTab
save(corTab, file = "../temp/03_CorrelationTable.RData")

#' # Filter data ####
#' ***
#' I only want to check the 11 loci that have more than associated 2 SNPs (pval<1e-6). 
#' 
load("../results/02_LociOverallPhenotypes_filtered.RData")

SNPs_sig = result.0[invalidAssocs  == F & pval<=1e-6,unique(markername)]
result.1 = copy(result.0)
result.1 = result.1[markername %in% SNPs_sig,]

result.2 = foreach(i=1:dim(result.5)[1])%do%{
  #i=1
  myRow = result.5[i,]
  dummy = copy(result.1)
  dummy = dummy[chr == myRow$chr]
  dummy = dummy[bp_hg19 >= myRow$region_start]
  dummy = dummy[bp_hg19 <= myRow$region_end]
  dummy[,region := myRow$region]
  dummy[,candidateGene := myRow$candidateGene]
  dummy[,bestPheno := myRow$phenotype]
  dummy
}
result.2 = rbindlist(result.2)
table(result.2$phenotype,result.2$chr)
setorder(result.2,chr,bp_hg19)

#' # Get Strata info per SNP ####
#' ***
#' 
result.2[,table(phenotype,bestPheno)]
dummy = copy(result.2)
setorder(dummy,pval)
dummy = dummy[!is.na(pval),]
dummy = dummy[!duplicated(markername),]
matched = match(result.2$markername,dummy$markername)
result.2[,bestPheno2 := dummy[matched,phenotype]]
result.2[,table(bestPheno, bestPheno2)]

result.2[, bestSex := "combined"]
result.2[grepl("_male",bestPheno2), bestSex := "males"]
result.2[grepl("_female",bestPheno2), bestSex := "females"]
result.2[,table(bestPheno2,bestSex)]

result.2[, bestStatin := "combined"]
result.2[grepl("_free",bestPheno2), bestStatin := "free"]
result.2[grepl("_treated",bestPheno2), bestStatin := "treated"]
result.2[,table(bestPheno2,bestStatin)]

result.2[, sex := "combined"]
result.2[grepl("_male",phenotype), sex := "males"]
result.2[grepl("_female",phenotype), sex := "females"]
result.2[,table(bestPheno2,sex)]

result.2[, statin := "combined"]
result.2[grepl("_free",phenotype), statin := "free"]
result.2[grepl("_treated",phenotype), statin := "treated"]
result.2[,table(bestPheno2,statin)]

# result.2[,corSex := corTab[1,correlation]]
# result.2[bestStatin == "free",corSex := corTab[2,correlation]]
# result.2[bestStatin == "treated",corSex := corTab[3,correlation]]
# result.2[,corStatin := corTab[4,correlation]]
# result.2[bestSex == "males",corStatin := corTab[5,correlation]]
# result.2[bestSex == "females",corStatin := corTab[6,correlation]]
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

save(result.2,file="../temp/03_IATest_input.RData")

#' # 3-way Interaction Tests ####
#' ***
#' ## Run Test ####
IATab_3way = ThreeWayInteractionTest_jp(data = result.2,
                                  pheno1 = "PCSK9_males_treated",
                                  pheno2 = "PCSK9_females_treated",
                                  pheno3 = "PCSK9_males_free",
                                  pheno4 = "PCSK9_females_free")

matched = match(IATab_3way$markername,result.2$markername)
table(is.na(matched))
table(IATab_3way$markername == result.2$markername[matched])
IATab_3way[,candidateGene := result.2$candidateGene[matched]]

save(IATab_3way,file="../temp/03_IATest_output_3way.RData")

#' ## Correct for multiple testing ####
IATab_3way = IATab_3way[!is.na(IA_diff),]

myFDR1 = addHierarchFDR(pvalues = IATab_3way$IA_pval, categs = IATab_3way$candidateGene,quiet = F)
table(myFDR1$hierarch_fdr5proz,myFDR1$category)

IATab_3way[,IA_pval_adj := myFDR1$fdr_level1]
IATab_3way[,IA_hierarch_fdr5proz := myFDR1$hierarch_fdr5proz]
IATab_3way[,table(IA_hierarch_fdr5proz,candidateGene)]
IATab_3way[,table(IA_pval<0.05,candidateGene)]

IATab_3way[,type2 := "unspecific"]
IATab_3way[IA_pval<0.05,type2 := "nom. sig. interaction"]
IATab_3way[IA_hierarch_fdr5proz==T,type2 := "interaction"]

result.6 = IATab_3way[markername %in% result.5$markername,c(1,30:36)]
result.6[,IA_pval_adj2 := p.adjust(p=IA_pval,method = "fdr")]
result.6

#' **Summary** 3-way interaction:
#' 
#' - **ALOX5**: significant 3-way interaction in 38 SNPs, including the lead SNP
#'    - no association in statin-free females
#'    - no association in statin-treated females
#'    - no association in statin-treated males
#'    - significant association in statin-free males
#' - using lead SNPs only does not change the result
#'      

save(IATab_3way,file = "../results/03_InteractionTests_3way.RData")

#' # 2-way Interaction Tests ####
#' ***
#' The two way interaction should be tested with respect to the best-associated phenotype per loci. E.g. if the best-associated phenotype is statin-free females, then I compare the effects in statin-free males and females, and in statin-treated and statin-free females.  
#' 
#' ## Run Test for sex ####
IATab_2way_sex = TwoWayInteractionTest_jp(data = result.2,
                                          pheno1 = "males",
                                          pheno2 = "females",
                                          type = "sexIA",
                                          corCol = "corSex")

matched = match(IATab_2way_sex$markername,result.2$markername)
table(is.na(matched))
table(IATab_2way_sex$markername == result.2$markername[matched])
IATab_2way_sex[,candidateGene := result.2$candidateGene[matched]]

#' ## Correct for multiple testing for sex ####
IATab_2way_sex = IATab_2way_sex[!is.na(IA_diff),]

myFDR2 = addHierarchFDR(pvalues = IATab_2way_sex$IA_pval, categs = IATab_2way_sex$candidateGene,quiet = F)
table(myFDR2$hierarch_fdr5proz,myFDR2$category)

IATab_2way_sex[,IA_pval_adj := myFDR2$fdr_level1]
IATab_2way_sex[,IA_hierarch_fdr5proz := myFDR2$hierarch_fdr5proz]
IATab_2way_sex[,table(IA_hierarch_fdr5proz,candidateGene)]
IATab_2way_sex[,table(IA_pval<0.05,candidateGene)]

IATab_2way_sex[,type2 := "unspecific"]
IATab_2way_sex[IA_hierarch_fdr5proz==T,type2 := "sex-interaction"]

result.7 = IATab_2way_sex[markername %in% result.5$markername,c(1,24:27,7,28:29)]
result.7[,IA_pval_adj2 := p.adjust(p=IA_pval,method = "fdr")]
result.7

IATab_2way_sex[IA_hierarch_fdr5proz==T,.N,candidateGene]

#' **Summary** 2-way sex interaction:
#' 
#' - **PRKAG2**: male-specific association
#' - **ALOX5**: male-specific association
#' - **SLCO1B1**: female-specfic association
#' - **NOS1**: male-specific association
#' - **PCSK9**: 34 mixed interactions (consistent effect direction)
#'    - 9 SNPs with genome-wide significance in both, but weaker effects in females
#'    - 9 SNPs with significant effect in males, but no effect in females
#'    - 16 SNPs with significant effect in males, but only nominal effect in females

IATab_2way_sex[IA_hierarch_fdr5proz==T & candidateGene %in% c("PRKAG2","ALOX5","NOS1"),type2 := "male-specific"]
IATab_2way_sex[IA_hierarch_fdr5proz==T & candidateGene %in% c("SLCO1B3"),type2 := "female-specific"]
table(IATab_2way_sex$type2)
save(IATab_2way_sex,file="../temp/03_IATest_output_2way_sex.RData")

#' ## Run Test for statin ####
IATab_2way_statin = TwoWayInteractionTest_jp(data = result.2,
                                             pheno1 = "treated",
                                             pheno2 = "free",
                                             type = "statinIA",
                                             corCol = "corStatin")

matched = match(IATab_2way_statin$markername,result.2$markername)
table(is.na(matched))
table(IATab_2way_statin$markername == result.2$markername[matched])
IATab_2way_statin[,candidateGene := result.2$candidateGene[matched]]

#' ## Correct for multiple testing for statin ####
IATab_2way_statin = IATab_2way_statin[!is.na(IA_diff),]

myFDR3 = addHierarchFDR(pvalues = IATab_2way_statin$IA_pval, categs = IATab_2way_statin$candidateGene,quiet = F)
table(myFDR3$hierarch_fdr5proz,myFDR3$category)

IATab_2way_statin[,IA_pval_adj := myFDR3$fdr_level1]
IATab_2way_statin[,IA_hierarch_fdr5proz := myFDR3$hierarch_fdr5proz]
IATab_2way_statin[,table(IA_hierarch_fdr5proz,candidateGene)]
IATab_2way_statin[,table(IA_pval<0.05,candidateGene)]

IATab_2way_statin[,type2 := "unspecific"]
IATab_2way_statin[IA_hierarch_fdr5proz==T,type2 := "interaction"]

result.8 = IATab_2way_statin[markername %in% result.5$markername,c(1,24:27,7,28:29)]
result.8[,IA_pval_adj2 := p.adjust(p=IA_pval,method = "fdr")]
result.8

IATab_2way_statin[IA_hierarch_fdr5proz==T,.N,candidateGene]

#' **Summary** 2-way statin interaction:
#' 
#' Here it makes a difference if you use all SNPs or just the lead SNPs!
#' 
#' - **KHDRBS2**: 7 significant interactions, lead SNP has only effect in statin-treated individuals
#' - **PRKAG2**: 4 significant interactions, lead SNP has only effect in statin-free males
#' - **ALOX5**: 60 significant interactions, lead SNP has only effect in statin-free males
#' - **PCSK9**: 28 nominal significant interactions, lead SNP has stronger effect in statin-free individuals
#' - **APOB**: 109 nominal significant interactions, lead SNP has only effect in statin-free individuals
#' - **JMJD1C**: 108 nominal significant interactions, lead SNP has only effect in statin-free individuals
#' - **TM6SF2**: 20 nominal significant interactions, lead SNP has only effect in statin-free individuals
#' 
IATab_2way_statin[IA_hierarch_fdr5proz==T & candidateGene %in% c("ALOX5","PRKAG2"),type2 := "free-specific"]
IATab_2way_statin[IA_hierarch_fdr5proz==T & candidateGene %in% c("KHDRBS2"),type2 := "treated-specific"]

table(IATab_2way_statin$type2)

save(IATab_2way_statin,file="../temp/03_IATest_output_2way_statin.RData")

#' ## Summary ####
#' ***
#' 
#' RERUN THIS WHEN YOU HAVE THE INDEPENDENT VARIANTS FROM GCTA - MAYBE THEN EASIER?
#' 
IATab_2way = rbind(IATab_2way_sex,IATab_2way_statin)
save(IATab_2way,file = "../results/03_InteractionTests_2way.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

