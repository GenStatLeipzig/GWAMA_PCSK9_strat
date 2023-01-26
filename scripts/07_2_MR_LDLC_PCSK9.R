#' ---
#' title: "Get associated loci for LDLC in males and females"
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
#' Define loci for my PCSK9 results
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/getSmallestDist.R")
source("../helperFunctions/MRfunction_jp.R")

#' # Load data ####
#' ***
#' Load data and filter for suggestive significant and valid SNPs. Finally, merge all significant SNPs into one object. 
#' 
LDLC_fem = fread(paste0(path_lipids,"/sex_specific_summary_stats/meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_FEMALE_with_N_1.gz"))
LDLC_fem = LDLC_fem[POOLED_ALT_AF >= 0.01,]
LDLC_fem = LDLC_fem[POOLED_ALT_AF <= 0.99,]
LDLC_fem = LDLC_fem[pvalue_neg_log10 >= 7.3,]
LDLC_fem[,phenotype := "LDLC_female"]
LDLC_fem[,METAL_Pvalue2 := as.numeric(METAL_Pvalue)]
LDLC_fem = LDLC_fem[METAL_Pvalue2 < 5e-8,]

LDLC_mal = fread(paste0(path_lipids,"/sex_specific_summary_stats/meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_MALE_with_N_1.gz"))
LDLC_mal = LDLC_mal[POOLED_ALT_AF >= 0.01,]
LDLC_mal = LDLC_mal[POOLED_ALT_AF <= 0.99,]
LDLC_mal = LDLC_mal[pvalue_neg_log10 >= 7.3,]
LDLC_mal[,phenotype := "LDLC_male"]
LDLC_mal[,METAL_Pvalue2 := as.numeric(METAL_Pvalue)]
LDLC_mal = LDLC_mal[METAL_Pvalue2 < 5e-8,]

LDLC_all = fread(paste0(path_lipids,"/trans_ancestry/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz"))
LDLC_all = LDLC_all[POOLED_ALT_AF >= 0.01,]
LDLC_all = LDLC_all[POOLED_ALT_AF <= 0.99,]
LDLC_all = LDLC_all[pvalue_neg_log10 >= 7.3,]
LDLC_all[,phenotype := "LDLC_all"]
LDLC_all[,METAL_Pvalue2 := as.numeric(METAL_Pvalue)]
LDLC_all = LDLC_all[METAL_Pvalue2 < 5e-8,]
LDLC_all[,lnBF := NULL]

result.1 = rbind(LDLC_fem,LDLC_mal,LDLC_all)
result.1[,chrPosPheno := paste(CHROM, POS_b37,phenotype,sep=":")]
table(duplicated(result.1$chrPosPheno))
table(result.1$phenotype,result.1$CHROM)
save(result.1, file = "../temp/07_IndexSNPs_LDLC.RData")
load("../temp/07_IndexSNPs_LDLC.RData")

#' # Get Top SNP per 500 kb range ###
#' ***
#' I want to reduce the list of associated SNPs by assigning SNPs to the top SNP of the region per phenotype. To do this, I will order the SNPs by their p-value, choose the SNP with lowest p-value as lead SNP, and assign all SNPs within 500 kb around the lead SNP to this region. This will be repeated until no SNPs can be assign, as the minimal pairwise distance is larger than 500 kb.  
#' 
subset1 = unique(result.1$phenotype)

result.2 = foreach(s1 = subset1) %do% {
  # s1 = subset1[1]
  subdata = copy(result.1)
  subdata = subdata[phenotype == s1, ]
  subset2 = unique(subdata$CHROM)
  
  result.22 = foreach(s2 = subset2) %do% {
    # s2 = subset2[1]
    subdata2 = copy(subdata)
    subdata2 = subdata2[CHROM == s2, ]
    
    setkey(subdata2, POS_b37)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, NR_SNPs := 0]
    }else{
      subdata2[, keep := NA]
      subdata2[, NR_SNPs := as.numeric(NA)]
      
      smallestDist = getSmallestDist(subdata2[,POS_b37])
      while(smallestDist < 500000) {
        maxLogP = max(subdata2[is.na(keep), pvalue_neg_log10])
        myPOS_b37 = subdata2[maxLogP == pvalue_neg_log10 & is.na(keep), POS_b37]
        if(length(myPOS_b37)>1){
          myPOS_b37 = myPOS_b37[1]
        }
        subdata2[POS_b37 == myPOS_b37, keep := T]
        
        #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
        myFilt = (subdata2[, POS_b37] < (myPOS_b37 - 500000)) | 
          (subdata2[, POS_b37] > (myPOS_b37 + 500000)) | 
          subdata2[, keep] 
        myFilt[is.na(myFilt)] = FALSE
        subdata2 = subdata2[myFilt == TRUE, ]
        
        subdata2[POS_b37 == myPOS_b37, NR_SNPs := sum(myFilt==F)]
        smallestDist = getSmallestDist(subdata2[, POS_b37])
      }
      
      #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
      subdata2[is.na(keep), NR_SNPs := 0]
      subdata2[is.na(keep), keep := TRUE]
    }
    
    subdata2
  }
  result.22 = rbindlist(result.22)
  table(result.22$phenotype)
  table(result.22$CHROM)
  result.22
}
result.2 = rbindlist(result.2)
table(result.2$phenotype,result.2$CHROM)

#add regions
result.2[, region_start := POS_b37 - 500000]
result.2[, region_end := POS_b37 + 500000]
result.2

#' # Locus collapsing by phenotype ###
#' ***
#' I want to collapse loci if their regions are overlapping, although their base positions differ more then 500 kb. This could be the case in large associated regions. 
#' 

result.3 = foreach(s1 = subset1) %do% {
  # s1 = subset1[1]
  subdata = copy(result.2)
  subdata = subdata[phenotype == s1, ]
  subset2 = unique(subdata$CHROM)
  
  result.32 = foreach(s2 = subset2) %do% {
    # s2 = subset2[1]
    subdata2 = copy(subdata)
    subdata2 = subdata2[CHROM == s2, ]
    
    setkey(subdata2, POS_b37)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, region := 1]
    }else{
      subdata2[, keep := NA]
      subdata2[, region := as.numeric(NA)]
      
      subdata2[1, region := 1]
      foreach(l = c(2:nrow(subdata2))) %do% {
        if (subdata2[l, region_start] <= subdata2[l-1, region_end]) {
          subdata2[l, region := subdata2[l-1, region]]
        } else { subdata2[l, region := subdata2[l-1, region] + 1]}
      }
      
      #keep best SNPs per region
      allRegions = unique(subdata2[, region])
      
      subresult = foreach(r = allRegions) %do% {
        subdata.2 = copy(subdata2)
        subdata.2 = subdata.2[region == r, ]
        maxLogP = max(subdata.2[, pvalue_neg_log10])
        subdata.2[pvalue_neg_log10 == maxLogP, region_start := min(subdata.2[, region_start])]
        subdata.2[pvalue_neg_log10 == maxLogP, region_end := max(subdata.2[, region_end])]
        subdata.2[pvalue_neg_log10 == maxLogP, NR_SNPs := sum(subdata.2[,NR_SNPs]) + nrow(subdata.2) - 1]
        subdata.2[pvalue_neg_log10 == maxLogP, ]
      }
      subresult = rbindlist(subresult)
      subresult[,region := paste(s2,region,sep="::")]
      subresult
    }
    
  }
  result.32 = rbindlist(result.32,fill = T)
  table(result.32$phenotype)
  table(result.32$CHROM)
  result.32
  
}
result.3 = rbindlist(result.3,fill = T)
table(result.3$phenotype,result.3$CHROM)
result.3[, region2 := c(1:nrow(result.3))]
result.3

names(result.3)
result.3

save(result.3,file="../temp/07_IndexSNPs_LDLC2.RData")
load("../temp/07_IndexSNPs_LDLC2.RData")

#' # Check for overlapping ###
#' ***
#' I only want those SNPs that show no genome-wide significance with PCSK9, so I exclude all SNPs in my PCSK9 regions
load("../results/05_GCTA_COJO.RData")
IndepSignals
IndepSignals = IndepSignals[c(1,11:20),]

result.4 = foreach(i=1:dim(IndepSignals)[1])%do%{
  #i=1
  myRow = IndepSignals[i,]
  
  result.41 = copy(result.3)
  result.41 = result.41[CHROM == myRow$Chr,]
  result.41 = result.41[POS_b37 >= myRow$region_start,]
  result.41 = result.41[POS_b37 <= myRow$region_end,]
  result.41

}
result.4 = rbindlist(result.4,fill = T)
result.4

result.5 = copy(result.3)
result.5 = result.5[region2 %nin% result.4$region2,]
result.5[,table(phenotype)]

#' Okay, these 468 SNPs can be used in my MR approach :)
#' 
#' # Load & match PCSK9 data 
fem_adj = fread("../data/SumStat_PCSK9_females_230120.txt.gz") 
fem_adj = fem_adj[invalidAssocs  == F ,]
fem_fre = fread("../data/SumStat_PCSK9_females_free_230120.txt.gz") 
fem_fre = fem_fre[invalidAssocs  == F ,]

mal_adj = fread("../data/SumStat_PCSK9_males_230120.txt.gz") 
mal_adj = mal_adj[invalidAssocs  == F ,]
mal_fre = fread("../data/SumStat_PCSK9_males_free_230120.txt.gz") 
mal_fre = mal_fre[invalidAssocs  == F ,]

fre = fread("../data/SumStat_PCSK9_free_230120.txt.gz") 
fre = fre[invalidAssocs  == F ,]

result.5[,chrPosEAOA := paste0(CHROM,POS_b37,ALT,REF, sep=":")]
result.5[,chrPosOAEA := paste0(CHROM,POS_b37,REF,ALT, sep=":")]
result.5_fem = copy(result.5)
result.5_fem = result.5_fem[grepl("female",phenotype)]
result.5_mal = copy(result.5)
result.5_mal = result.5_mal[grepl("_male",phenotype)]
result.5_all = copy(result.5)
result.5_all = result.5_all[grepl("all",phenotype)]
result.5_fem[,table(duplicated(POS_b37))]
result.5_mal[,table(duplicated(POS_b37))]
result.5_all[,table(duplicated(POS_b37))]

fem_adj[,chrPosEAOA := paste0(chr,bp_hg19, EA, OA, sep=":")]
fem_fre[,chrPosEAOA := paste0(chr,bp_hg19, EA, OA, sep=":")]
mal_adj[,chrPosEAOA := paste0(chr,bp_hg19, EA, OA, sep=":")]
mal_fre[,chrPosEAOA := paste0(chr,bp_hg19, EA, OA, sep=":")]
fre[,chrPosEAOA := paste0(chr,bp_hg19, EA, OA, sep=":")]

fem_adj = fem_adj[chrPosEAOA %in% c(result.5_fem$chrPosEAOA,result.5_fem$chrPosOAEA),]
fem_fre = fem_fre[chrPosEAOA %in% c(result.5_fem$chrPosEAOA,result.5_fem$chrPosOAEA),]
mal_adj = mal_adj[chrPosEAOA %in% c(result.5_mal$chrPosEAOA,result.5_mal$chrPosOAEA),]
mal_fre = mal_fre[chrPosEAOA %in% c(result.5_mal$chrPosEAOA,result.5_mal$chrPosOAEA),]
fre = fre[chrPosEAOA %in% c(result.5_all$chrPosEAOA,result.5_all$chrPosOAEA),]

fem_adj[,table(duplicated(bp_hg19))]
fem_fre[,table(duplicated(bp_hg19))]
mal_adj[,table(duplicated(bp_hg19))]
mal_fre[,table(duplicated(bp_hg19))]
fre[,table(duplicated(bp_hg19))]



#' ## Matching female adjusted ####
#' 
matched1 = match(fem_adj$bp_hg19,result.5_fem$POS_b37)
table(is.na(matched1))
fem_adj = fem_adj[,c(16,1:4,6,10:12,8)]
names(fem_adj)
names(fem_adj) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
fem_adj[,EA.LDLC := result.5_fem[matched1,ALT]]
fem_adj[,EAF.LDLC := result.5_fem[matched1,POOLED_ALT_AF]]
fem_adj[,beta.LDLC := result.5_fem[matched1,METAL_Effect]]
fem_adj[,SE.LDLC := result.5_fem[matched1,METAL_StdErr]]
fem_adj[,logP.LDLC := result.5_fem[matched1,pvalue_neg_log10]]
fem_adj[,N.LDLC := result.5_fem[matched1,N]]
fem_adj[EA.LDLC !=EA, EAF.LDLC := 1-EAF.LDLC]
fem_adj[EA.LDLC !=EA, beta.LDLC := (-1)*beta.LDLC]
fem_adj[,EA:=NULL]
plot(fem_adj$EAF.PCSK9,fem_adj$EAF.LDLC)

#' ## Matching female free ####
#' 
matched1 = match(fem_fre$bp_hg19,result.5_fem$POS_b37)
table(is.na(matched1))
fem_fre = fem_fre[,c(16,1:4,6,10:12,8)]
names(fem_fre)
names(fem_fre) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
fem_fre[,EA.LDLC := result.5_fem[matched1,ALT]]
fem_fre[,EAF.LDLC := result.5_fem[matched1,POOLED_ALT_AF]]
fem_fre[,beta.LDLC := result.5_fem[matched1,METAL_Effect]]
fem_fre[,SE.LDLC := result.5_fem[matched1,METAL_StdErr]]
fem_fre[,logP.LDLC := result.5_fem[matched1,pvalue_neg_log10]]
fem_fre[,N.LDLC := result.5_fem[matched1,N]]
fem_fre[EA.LDLC !=EA, EAF.LDLC := 1-EAF.LDLC]
fem_fre[EA.LDLC !=EA, beta.LDLC := (-1)*beta.LDLC]
fem_fre[,EA:=NULL]
plot(fem_fre$EAF.PCSK9,fem_fre$EAF.LDLC)

#' ## Matching male adjusted ####
#' 
matched1 = match(mal_adj$bp_hg19,result.5_mal$POS_b37)
table(is.na(matched1))
mal_adj = mal_adj[,c(16,1:4,6,10:12,8)]
names(mal_adj)
names(mal_adj) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
mal_adj[,EA.LDLC := result.5_mal[matched1,ALT]]
mal_adj[,EAF.LDLC := result.5_mal[matched1,POOLED_ALT_AF]]
mal_adj[,beta.LDLC := result.5_mal[matched1,METAL_Effect]]
mal_adj[,SE.LDLC := result.5_mal[matched1,METAL_StdErr]]
mal_adj[,logP.LDLC := result.5_mal[matched1,pvalue_neg_log10]]
mal_adj[,N.LDLC := result.5_mal[matched1,N]]
mal_adj[EA.LDLC !=EA, EAF.LDLC := 1-EAF.LDLC]
mal_adj[EA.LDLC !=EA, beta.LDLC := (-1)*beta.LDLC]
mal_adj[,EA:=NULL]
plot(mal_adj$EAF.PCSK9,mal_adj$EAF.LDLC)

#' ## Matching male free ####
#' 
matched1 = match(mal_fre$bp_hg19,result.5_mal$POS_b37)
table(is.na(matched1))
mal_fre = mal_fre[,c(16,1:4,6,10:12,8)]
names(mal_fre)
names(mal_fre) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
mal_fre[,EA.LDLC := result.5_mal[matched1,ALT]]
mal_fre[,EAF.LDLC := result.5_mal[matched1,POOLED_ALT_AF]]
mal_fre[,beta.LDLC := result.5_mal[matched1,METAL_Effect]]
mal_fre[,SE.LDLC := result.5_mal[matched1,METAL_StdErr]]
mal_fre[,logP.LDLC := result.5_mal[matched1,pvalue_neg_log10]]
mal_fre[,N.LDLC := result.5_mal[matched1,N]]
mal_fre[EA.LDLC !=EA, EAF.LDLC := 1-EAF.LDLC]
mal_fre[EA.LDLC !=EA, beta.LDLC := (-1)*beta.LDLC]
mal_fre[,EA:=NULL]
plot(mal_fre$EAF.PCSK9,mal_fre$EAF.LDLC)

#' ## Matching free ####
#' 
matched1 = match(fre$bp_hg19,result.5_all$POS_b37)
table(is.na(matched1))
fre = fre[,c(16,1:4,6,10:12,8)]
names(fre)
names(fre) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
fre[,EA.LDLC := result.5_all[matched1,ALT]]
fre[,EAF.LDLC := result.5_all[matched1,POOLED_ALT_AF]]
fre[,beta.LDLC := result.5_all[matched1,METAL_Effect]]
fre[,SE.LDLC := result.5_all[matched1,METAL_StdErr]]
fre[,logP.LDLC := result.5_all[matched1,pvalue_neg_log10]]
fre[,N.LDLC := result.5_all[matched1,N]]
fre[EA.LDLC !=EA, EAF.LDLC := 1-EAF.LDLC]
fre[EA.LDLC !=EA, beta.LDLC := (-1)*beta.LDLC]
fre[,EA:=NULL]
plot(fre$EAF.PCSK9,fre$EAF.LDLC)

#' ## Save data ####
#' 
myMRTab = rbind(fem_adj,fem_fre,mal_adj,mal_fre,fre)
save(myMRTab,file="../temp/07_DataSet_LDLC_PCSK9.RData")
load("../temp/07_DataSet_LDLC_PCSK9.RData")

#' # Get Ratios ####
#' ***

test_LDL_PCSK9 = myMRTab[,MRfunction_jp(betaX = beta.LDLC,seX = SE.LDLC,betaY = beta.PCSK9,seY = SE.PCSK9)]

myMRTab[,beta.Ratio := test_LDL_PCSK9$beta_IV]
myMRTab[,SEst.Ratio := test_LDL_PCSK9$se_IV1]
myMRTab[,SEnd.Ratio := test_LDL_PCSK9$se_IV2]
myMRTab[,Pst.Ratio := test_LDL_PCSK9$p_IV1]
myMRTab[,Pnd.Ratio := test_LDL_PCSK9$p_IV2]

#' # Get meta per phenotype ####
#' ***
table(myMRTab$phenotype)

mymod_fa<-myMRTab[phenotype == "PCSK9_females",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_ff<-myMRTab[phenotype == "PCSK9_females_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_ma<-myMRTab[phenotype == "PCSK9_males",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_mf<-myMRTab[phenotype == "PCSK9_males_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_f<-myMRTab[phenotype == "PCSK9_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]

dummy1 = myMRTab[,.N,phenotype]
dummy1[,exposure :="LDL-C"]
dummy1[,beta_IV := c(mymod_fa$TE.fixed,mymod_ff$TE.fixed,mymod_ma$TE.fixed,mymod_mf$TE.fixed,mymod_f$TE.fixed)]
dummy1[,SE_IV := c(mymod_fa$seTE.fixed,mymod_ff$seTE.fixed,mymod_ma$seTE.fixed,mymod_mf$seTE.fixed,mymod_f$seTE.fixed)]
dummy1[,pval_IV := c(mymod_fa$pval.fixed,mymod_ff$pval.fixed,mymod_ma$pval.fixed,mymod_mf$pval.fixed,mymod_f$pval.fixed)]
dummy1[,Q_IV := c(mymod_fa$Q,mymod_ff$Q,mymod_ma$Q,mymod_mf$Q,mymod_f$Q)]
dummy1[,pvalQ_IV := c(mymod_fa$pval.Q,mymod_ff$pval.Q,mymod_ma$pval.Q,mymod_mf$pval.Q,mymod_f$pval.Q)]
dummy1

myMRTab1 = copy(myMRTab)
setorder(myMRTab1,-logP.LDLC)
myMRTab1 = myMRTab1[, .SD[1:20], by=phenotype]

mymod_fa<-myMRTab1[phenotype == "PCSK9_females",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_ff<-myMRTab1[phenotype == "PCSK9_females_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_ma<-myMRTab1[phenotype == "PCSK9_males",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_mf<-myMRTab1[phenotype == "PCSK9_males_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_f<-myMRTab1[phenotype == "PCSK9_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]

dummy2 = myMRTab1[,.N,phenotype]
dummy2[,exposure :="LDL-C"]
dummy2[,beta_IV := c(mymod_fa$TE.fixed,mymod_ff$TE.fixed,mymod_ma$TE.fixed,mymod_mf$TE.fixed,mymod_f$TE.fixed)]
dummy2[,SE_IV := c(mymod_fa$seTE.fixed,mymod_ff$seTE.fixed,mymod_ma$seTE.fixed,mymod_mf$seTE.fixed,mymod_f$seTE.fixed)]
dummy2[,pval_IV := c(mymod_fa$pval.fixed,mymod_ff$pval.fixed,mymod_ma$pval.fixed,mymod_mf$pval.fixed,mymod_f$pval.fixed)]
dummy2[,Q_IV := c(mymod_fa$Q,mymod_ff$Q,mymod_ma$Q,mymod_mf$Q,mymod_f$Q)]
dummy2[,pvalQ_IV := c(mymod_fa$pval.Q,mymod_ff$pval.Q,mymod_ma$pval.Q,mymod_mf$pval.Q,mymod_f$pval.Q)]
dummy2

dummy1[,comment:="unfiltered"]
dummy2[,comment:="Top 20 SNPs"]

myMRTab_meta = rbind(dummy1,dummy2)

#' # Compare sexes ####
#' ***
interactionTest(mean1 = dummy1$beta_IV[3],se1 = dummy1$SE_IV[3],
                mean2 = dummy1$beta_IV[1],se2 = dummy1$SE_IV[1])
interactionTest(mean1 = dummy1$beta_IV[4],se1 = dummy1$SE_IV[4],
                mean2 = dummy1$beta_IV[2],se2 = dummy1$SE_IV[2])

interactionTest(mean1 = dummy2$beta_IV[4],se1 = dummy2$SE_IV[4],
                mean2 = dummy2$beta_IV[2],se2 = dummy2$SE_IV[2])
interactionTest(mean1 = dummy2$beta_IV[5],se1 = dummy2$SE_IV[5],
                mean2 = dummy2$beta_IV[3],se2 = dummy2$SE_IV[3])

#' # Make Plots ####
#' ***
#' male/female
#' 
#' adjusted/free
#' 
#' all/ top20
#' 
plotdata1 = copy(myMRTab)
plotdata1[beta.LDLC<0,beta.PCSK9 := (-1)*beta.PCSK9]
plotdata1[beta.LDLC<0,beta.LDLC := (-1)*beta.LDLC]

linedata1 = copy(dummy1)
linedata1[,int := 0]
linedata1[,upper := beta_IV+1.96*SE_IV]
linedata1[,lower := beta_IV-1.96*SE_IV]

myPlot1 = ggplot(plotdata1, aes(x=beta.LDLC, y=beta.PCSK9)) +
  facet_wrap(~phenotype, scales = "free") +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1)+
  #geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", size=1)+
  geom_point(size=2.5,color="steelblue")+
  # geom_errorbar(data = plotdata1,
  #                aes(ymin = beta.PCSK9- 1.96*SE.PCSK9, ymax = beta.PCSK9 + 1.96*SE.PCSK9)) +
  # geom_errorbarh(data = plotdata1,
  #               aes(xmin = beta.LDLC- 1.96*SE.LDLC, xmax = beta.LDLC+ 1.96*SE.LDLC)) +
  geom_abline(data = linedata1,
              aes(intercept = int,slope=beta_IV),
              color="firebrick4", linetype="dashed", size=1.15)+
  geom_abline(data = linedata1,
              aes(intercept = int,slope=lower),
              color="coral", linetype="dotted", size=1)+
  geom_abline(data = linedata1,
              aes(intercept = int,slope=upper),
              color="coral", linetype="dotted", size=1)+
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold")) + 
  labs(x="Effect on LDLC",
       y = "Effect on PCSK9")
myPlot1

##

plotdata2 = copy(myMRTab1)
plotdata2[beta.LDLC<0,beta.PCSK9 := (-1)*beta.PCSK9]
plotdata2[beta.LDLC<0,beta.LDLC := (-1)*beta.LDLC]

linedata2 = copy(dummy2)
linedata2[,int := 0]
linedata2[,upper := beta_IV+1.96*SE_IV]
linedata2[,lower := beta_IV-1.96*SE_IV]

myPlot2 = ggplot(plotdata2, aes(x=beta.LDLC, y=beta.PCSK9)) +
  facet_wrap(~phenotype, scales = "free",ncol = 2) +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1)+
  #geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", size=1)+
  geom_errorbar(data = plotdata2,
                 aes(ymin = beta.PCSK9- 1.96*SE.PCSK9, ymax = beta.PCSK9 + 1.96*SE.PCSK9)) +
  geom_errorbarh(data = plotdata2,
                aes(xmin = beta.LDLC- 1.96*SE.LDLC, xmax = beta.LDLC+ 1.96*SE.LDLC)) +
  geom_point(size=2.5,color="steelblue")+
  geom_abline(data = linedata2,
              aes(intercept = int,slope=beta_IV),
              color="firebrick4", linetype="dashed", size=1.15)+
  geom_abline(data = linedata2,
              aes(intercept = int,slope=lower),
              color="coral", linetype="dotted", size=1)+
  geom_abline(data = linedata2,
              aes(intercept = int,slope=upper),
              color="coral", linetype="dotted", size=1)+
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold")) + 
  labs(x="Effect on LDLC",
       y = "Effect on PCSK9")
myPlot2

tiff(filename = "../figures/SupplementalFigure_MRScatterPlot_allSNPs_reverse.tiff", 
     width = 3000, height = 2700, res=300, compression = 'lzw')
myPlot1
dev.off()

tiff(filename = "../figures/SupplementalFigure_MRScatterPlot_PCSK9SNPs_reverse.tiff", 
     width = 2000, height = 1800, res=300, compression = 'lzw')
myPlot2
dev.off()

#' # Save results ####
#' ***

save(myMRTab,file="../results/07_MR_LDLC_PCSK9_SNPwise.RData")
save(myMRTab_meta,file="../results/07_MR_LDLC_PCSK9_summary.RData")


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

