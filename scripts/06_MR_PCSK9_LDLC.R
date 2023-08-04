#' ---
#' title: "MR 1: PCSK9 --> LDL-C"
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
#' I want to estimate the causal effect of PCSK9 on LDL-C.
#' 
#' I test this four times in all subgroups using the following SNP sets: 
#' 
#' - all 14 SNPs
#' - all SNPs that are at least suggestive significant (different SNPs per group)
#' - only the 4 SNPs at the *PCSK9* gene locus
#' - only the 3 SNPs at the *PCSK9* gene locus without rs11583680 (sig. 3-way interaction)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_forostar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))
source("../helperFunctions/ThreeWayInteractionTest_jp.R")
source("../helperFunctions/TwoWayInteractionTest_jp.R")

#' # Load ####
#' ***
load("../temp/04_IATest_input.RData")

load("../temp/05_OtherGWASs.RData")
myOtherGWAS = myOtherGWAS[grepl("LDL",phenotype)]
myOtherGWAS = myOtherGWAS[SNP %in% result.2$rsID]

data_GWAS = result.2[rsID %in% myOtherGWAS$SNP,]
myPhenos = result.2[,unique(phenotype)]
data_GWAS = data_GWAS[phenotype %in% c(myPhenos)]

#' # Match ####
#' ***
myLipids = copy(myOtherGWAS)
myLipids[,dumID := paste0(bp_hg19,"_male")]
myLipids[grepl("FEMALE",phenotype),dumID := paste0(bp_hg19,"_female")]
myLipids[grepl("ALL",phenotype),dumID := paste0(bp_hg19,"_combined")]

data_GWAS[,dumID := paste0(bp_hg19,"_combined")]
data_GWAS[grepl("_female",phenotype),dumID := paste0(bp_hg19,"_female")]
data_GWAS[grepl("_male",phenotype),dumID := paste0(bp_hg19,"_male")]

matched = match(data_GWAS$dumID,myLipids$dumID)
table(is.na(matched))
myLipids = myLipids[matched]
table(myLipids$dumID == data_GWAS$dumID)

#' # Check allele ####
#' ***
plot(data_GWAS$EAF,myLipids$EAF)
abline(0,1)

table(myLipids$EA == data_GWAS$EA)
table(myLipids$OA == data_GWAS$EA)
filt = myLipids$EA != data_GWAS$EA

EA1 = myLipids$EA
OA1 = myLipids$OA
myLipids[filt, EA := OA1[filt]]
myLipids[filt, OA := EA1[filt]]
myLipids[filt, EAF     := 1-EAF    ]
myLipids[filt, beta:= beta     *(-1)]

plot(data_GWAS$EAF,myLipids$EAF)
abline(0,1)

table(myLipids$EA == data_GWAS$EA)

#' # Merge ####
#' ***
myMRTab = copy(data_GWAS)
myMRTab = myMRTab[,c(17,16,1,2,3,4,6,10,11,12,8)]
names(myMRTab)
names(myMRTab) = c("gene","phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")

myMRTab[,EAF.LDLC := myLipids$EAF]
myMRTab[,beta.LDLC := myLipids$beta]
myMRTab[,SE.LDLC := myLipids$SE]
myMRTab[,P.LDLC := myLipids$pval]
myMRTab[,N.LDLC := myLipids$nSamples]

myMRTab
save(myMRTab,file = "../temp/06_MR_input.RData")

#' # MR analysis ####
#' ***
myPhenos
mySettings = c("all SNPs","sug sig SNPs","PCSK9 SNPs","PCSK9 SNPs filtered")

myMRTab[,setting1 := T]
myMRTab[,setting2 := F]
myMRTab[P.PCSK9<1e-6, setting2 := T]
myMRTab[,setting3 := F]
myMRTab[gene == "PCSK9", setting3 := T]
myMRTab[,setting4 := F]
myMRTab[gene == "PCSK9" & SNP != "rs11583680:55505668:C:T", setting4 := T]

dumTab1 = foreach(i=1:length(myPhenos))%do%{
  #i=2
  message("Working on phenotype ",myPhenos[i])
  myMRTab1 = copy(myMRTab)
  myMRTab1 = myMRTab1[phenotype == myPhenos[i]]
  
  myPheno = myPhenos[i]
  myPheno = gsub("PCSK9_","",myPheno)
  myPheno = gsub("_"," ",myPheno)
  
  dumTab2 = foreach(j=1:length(mySettings))%do%{
    #j=2
    message("   in setting ",j,": ",mySettings[j])
    myMRTab2 = copy(myMRTab1)
    myFilt = paste0("setting",j)
    myMRTab2 = myMRTab2[, setting := get(myFilt)]
    myMRTab2 = myMRTab2[setting ==T,]
    if(dim(myMRTab2)[1]<=1){
      res = data.table(phenotype = myPhenos[i],
                       outcome = "LDLC",
                       setting = mySettings[j],
                       nSNPs = dim(myMRTab2)[1])
      
    }else{
      mrob = mr_input(bx = myMRTab2$beta.PCSK9,
                      bxse = myMRTab2$SE.PCSK9,
                      by = myMRTab2$beta.LDLC,
                      byse = myMRTab2$SE.LDLC,
                      snps = myMRTab2$SNP,
                      exposure = paste0("PCSK9 ",myPheno),
                      outcome = "LDLC")
      res1 = mr_ivw(mrob)
      res2 = mr_egger(mrob)
      
      plot1 = mr_plot(mrob,line = "ivw",orientate = T,interactive = F)
      plot2 = mr_plot(mrob,line = "egger",orientate = T, interactive = F)

      tag = format(Sys.time(), "%Y-%m-%d")
      tag2 = gsub("2023-","23-",tag)
      tag2 = gsub("-","",tag2)

      filename1 = paste0("../figures/07_MRPlots_",myPhenos[i],"_setting",j,"_IVW_",tag2,".tiff")
      filename2 = paste0("../figures/07_MRPlots_",myPhenos[i],"_setting",j,"_egger_",tag2,".tiff")

      tiff(filename = filename1,
           width = 1650, height = 1350, res=250, compression = 'lzw')
      print(plot1)
      dev.off()
      tiff(filename = filename2,
           width = 1650, height = 1350, res=250, compression = 'lzw')
      print(plot2)
      dev.off()
      
      res = data.table(phenotype = myPhenos[i],
                       outcome = c(res1@Outcome),
                       setting = mySettings[j],
                       nSNPs = dim(myMRTab2)[1],
                       beta_IVW = c(res1@Estimate),
                       SE_IVW = c(res1@StdError),
                       pval_IVW = c(res1@Pvalue),
                       HeteroStat_IVW = c(res1@Heter.Stat[1]),
                       HeteroStat_pval_IVW = c(res1@Heter.Stat[2]),
                       beta_egger = c(res2@Estimate),
                       SE_egger = c(res2@StdError.Est),
                       pval_egger = c(res2@Pvalue.Est),
                       beta_egger_int = c(res2@Intercept),
                       SE_egger_int = c(res2@StdError.Int),
                       pval_egger_int = c(res2@Pvalue.Int),
                       HeteroStat_egger = c(res2@Heter.Stat[1]),
                       HeteroStat_pval_egger = c(res2@Heter.Stat[2]))
      
    }
        
    res
  }
  dumTab2 = rbindlist(dumTab2, fill = TRUE)
  dumTab2
}

myMRTab_results = rbindlist(dumTab1)
save(myMRTab_results,file = "../results/06_MRresults.RData")

#' # Interaction Test ####
#' ***
#' I want to perform all 3 tests: 
#' 
#' - 3-way interaction
#' - 2-way interaction for sex
#' - 2-way interaction for statin treatment
#' 
#' ## Run tests #### 
IATab_3way = ThreeWayInteractionTest_jp(data = myMRTab_results,
                                        pheno1 = "PCSK9_males_treated",
                                        pheno2 = "PCSK9_females_treated",
                                        pheno3 = "PCSK9_males_free",
                                        pheno4 = "PCSK9_females_free",
                                        data_type = "MR",
                                        suffix = "IVW")

IATab_2way_sex = TwoWayInteractionTest_jp(data = myMRTab_results,
                                          pheno1 = "males",
                                          pheno2 = "females",
                                          type = "sexIA",
                                          useBestPheno = F,
                                          corCol = "corSex",
                                          data_type = "MR", suffix = "IVW")

IATab_2way_statin = TwoWayInteractionTest_jp(data = myMRTab_results,
                                          pheno1 = "treated",
                                          pheno2 = "free",
                                          type = "statinIA",
                                          useBestPheno = F,
                                          corCol = "corStatin",
                                          data_type = "MR", suffix = "IVW")

#' ## Combine data sets ####
IATab = rbind(IATab_3way,IATab_2way_sex,IATab_2way_statin,fill = T)
IATab = IATab[,c(1,18:21,2:17)]
IATab[,c(1:5)]

#' ## Correct for multiple testing ####
myFDR = addHierarchFDR(pvalues = IATab$IA_pval, categs = IATab$setting,quiet = F)
IATab[,IA_pval_adj := myFDR$fdr_level1]
IATab[,IA_hierarch_fdr5proz := myFDR$hierarch_fdr5proz]
IATab[,table(IA_hierarch_fdr5proz,setting)]

#' ## Summary ####
IATab[IA_hierarch_fdr5proz==T ,]
IATab[IA_pval <0.05 ,]

#' ## Save ####
IATab[,type := c(rep("3way",4),rep("sexIA",4),rep("statinIA",4))]
IATab = IATab[,c(24,1:5,22,23,6:21)]
save(IATab, file = "../results/06_MR_interaction.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
