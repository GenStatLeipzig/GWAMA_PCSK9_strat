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
ToDoList[,pheno := gsub("_230120.txt.gz","",pheno)]

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
#' I only want to check the 13 loci that I used in the GCTA and Credible Set Analyses. 
#' 
load("../results/02_LociOverallPhenotypes_filtered.RData")

SNPs_sig = result.0[invalidAssocs  == F & pval<=1e-6,unique(markername)]
result.1 = copy(result.0)
result.1 = result.1[markername %in% SNPs_sig,]
result.1 = result.1[chr %in% result.5$chr,]

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
result.1[markername %nin% result.2$markername,]
table(result.2$phenotype,result.2$chr)
setorder(result.2,chr,bp_hg19)

#' # Get Strata info per SNP ####
#' ***
#' 
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

result.2[,corSex := corTab[1,correlation]]
result.2[bestStatin == "free",corSex := corTab[2,correlation]]
result.2[bestStatin == "treated",corSex := corTab[3,correlation]]
result.2[,corStatin := corTab[4,correlation]]
result.2[bestSex == "males",corStatin := corTab[5,correlation]]
result.2[bestSex == "females",corStatin := corTab[6,correlation]]


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

#' **Summary** 3-way interaction:
#' 
#' - **PCSK9**: 12 SNPs with nominal significant interactions 
#'    - 8 at the beginning of the locus: gw sig for **treated males**, nom sig for free males and females, not sig for treated females
#'    - 4 at the end of the locus: sug sig for **free females**, nom sig for treated males, not sig for treated females and free males
#' - **MYNC**: 1 SNP with nominal significant interaction (sug sig in **free males**, not sig in the three other traits)
#' - **SLCO1B1**: 214 SNPs with nominal significant interaction (sug sig in **free females**, not sig in the other three traits)
#' - **SASH1**: 6 SNPs with significant interaction (sug sig in **treated females**, not sig in the other three traits)
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

IATab_2way_sex[IA_hierarch_fdr5proz==T,.N,candidateGene]

#' **Summary** 2-way sex interaction:
#' 
#' - **MYNC**: male-specific association
#' - **PLB1**: female-specific association
#' - **SASH1**: female-specific association
#' - **PRKAG2**: male-specific association
#' - **SLCO1B1**: female-specfic association
#' - **gene desert**: male-specific association
#' - **MACROD2,SNRPB2**: female-specfic association
#' 

IATab_2way_sex[IA_hierarch_fdr5proz==T & candidateGene %in% c("MYNC","PRKAG2","gene desert"),type2 := "male-specific"]
IATab_2way_sex[IA_hierarch_fdr5proz==T & candidateGene %in% c("PLB1","SASH1","SLCO1B1","MACROD2,SNRPB2"),type2 := "female-specific"]
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

IATab_2way_statin[IA_hierarch_fdr5proz==T,.N,candidateGene]

#' **Summary** 2-way statin interaction:
#' 
#' - **PCSK9**: treatment-related association
#' - **new1**: free-specific association
#' - **SASH1**: treatment-specific association
#' - **PRKAG2**: free-specific association
#' - **gene desert**: treatment-specific association
#' 
IATab_2way_statin[IA_hierarch_fdr5proz==T & candidateGene %in% c("new1","PRKAG2"),type2 := "free-specific"]
IATab_2way_statin[IA_hierarch_fdr5proz==T & candidateGene %in% c("SASH1","gene desert"),type2 := "treated-specific"]
IATab_2way_statin[IA_hierarch_fdr5proz==T & candidateGene %in% c("PCSK9"),type2 := "treated-related"]

table(IATab_2way_statin$type2)

save(IATab_2way_statin,file="../temp/03_IATest_output_2way_statin.RData")

#' ## Summary ####
#' ***
#' 
#' * no interaction whatsoever: APOB, TM6SF2, NOS1
#' * only sex interaction: PLB1, MYNC, MACROD2, SLCO1B1
#' * sex and statin interaction: CYP51A1P3, PRKAG2, gene desert
#' * sometimes interaction: PCSK9
#' 
IATab_2way = rbind(IATab_2way_sex,IATab_2way_statin)
save(IATab_2way,file = "../results/03_InteractionTests_2way.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

