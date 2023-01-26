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
#' I want to estimate the causal effect of PCSK9 on LDL-C in males and females. 
#' 
#' I test this with all 6 settings, using the independent signals per setting. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))
source("../helperFunctions/MRfunction_jp.R")

#' # Load ####
#' ***
load("../results/05_GCTA_COJO.RData")
myTab = copy(IndepSignals)
myTab = myTab[candidateGene == "PCSK9"]

load("../temp/06_OtherGWASs.RData")
myOtherGWAS = myOtherGWAS[grepl("LDL",phenotype)]
myOtherGWAS = myOtherGWAS[candidateGene == "PCSK9"]
myOtherGWAS = myOtherGWAS[bp_hg19 %in% myTab$bp]

#' # Match ####
#' ***
myLipids = copy(myOtherGWAS)
myLipids[,dumID := paste0(bp_hg19,"_male")]
myLipids[grepl("FEMALE",phenotype),dumID := paste0(bp_hg19,"_female")]
myLipids[grepl("ALL",phenotype),dumID := paste0(bp_hg19,"_combined")]

myTab[,dumID := paste0(bp,"_")]
myTab[num %in% c(1,2),dumID := paste0(bp,"_female")]
myTab[num %in% c(4:6),dumID := paste0(bp,"_male")]
myTab[num %in% c(3,7),dumID := paste0(bp,"_combined")]

matched = match(myTab$dumID,myLipids$dumID)
table(is.na(matched))
myLipids = myLipids[matched]
table(myLipids$dumID == myTab$dumID)

#' # Check allele ####
#' ***
plot(myTab$freq,myLipids$EAF)
abline(0,1)

table(myLipids$EA == myTab$refA)
table(myLipids$OA == myTab$refA)
filt = myLipids$EA != myTab$refA

EA1 = myLipids$EA
OA1 = myLipids$OA
myLipids[filt, EA := OA1[filt]]
myLipids[filt, OA := EA1[filt]]
myLipids[filt, EAF     := 1-EAF    ]
myLipids[filt, beta:= beta     *(-1)]

plot(myTab$freq,myLipids$EAF)
abline(0,1)

table(myLipids$EA == myTab$refA)

#' # Merge ####
#' ***
myMRTab = copy(myTab)
myMRTab = myMRTab[,c(15,2,1,3:9)]
names(myMRTab)
names(myMRTab) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
myMRTab[,N.PCSK9 := ceiling(N.PCSK9)]

myMRTab[,EAF.LDLC := myLipids$EAF]
myMRTab[,beta.LDLC := myLipids$beta]
myMRTab[,SE.LDLC := myLipids$SE]
myMRTab[,P.LDLC := myLipids$pval]
myMRTab[,N.LDLC := myLipids$nSamples]

myMRTab

#' # Get Ratio ####
#' ***

test_PCSK9_LDL = myMRTab[,MRfunction_jp(betaX = beta.PCSK9,seX = SE.PCSK9,betaY = beta.LDLC,seY = SE.LDLC)]
myMRTab[,beta.Ratio := test_PCSK9_LDL$beta_IV]
myMRTab[,SEst.Ratio := test_PCSK9_LDL$se_IV1]
myMRTab[,SEnd.Ratio := test_PCSK9_LDL$se_IV2]
myMRTab[,Pst.Ratio := test_PCSK9_LDL$p_IV1]
myMRTab[,Pnd.Ratio := test_PCSK9_LDL$p_IV2]

#' # Get meta per phenotype ####
#' ***
table(myMRTab$phenotype)

mymod_fa<-myMRTab[phenotype == "PCSK9_females",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_ff<-myMRTab[phenotype == "PCSK9_females_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_ma<-myMRTab[phenotype == "PCSK9_males",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_mf<-myMRTab[phenotype == "PCSK9_males_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod_af<-myMRTab[phenotype == "PCSK9_free",metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]

dummy = myMRTab[,.N,phenotype]
dummy[,outcome :="LDL-C"]
dummy[,beta_IV := c(mymod_fa$TE.fixed,mymod_ff$TE.fixed,mymod_ma$TE.fixed,mymod_mf$TE.fixed, mymod_af$TE.fixed, myMRTab$beta.Ratio[c(11,12)])]
dummy[,SE_IV := c(mymod_fa$seTE.fixed,mymod_ff$seTE.fixed,mymod_ma$seTE.fixed,mymod_mf$seTE.fixed, mymod_af$seTE.fixed, myMRTab$SEnd.Ratio[c(11,12)])]
dummy[,pval_IV := c(mymod_fa$pval.fixed,mymod_ff$pval.fixed,mymod_ma$pval.fixed,mymod_mf$pval.fixed, mymod_af$pval.fixed, myMRTab$Pnd.Ratio[c(11,12)])]
dummy[,Q_IV := c(mymod_fa$Q,mymod_ff$Q,mymod_ma$Q,mymod_mf$Q,mymod_af$Q, NA,NA)]
dummy[,pvalQ_IV := c(mymod_fa$pval.Q,mymod_ff$pval.Q,mymod_ma$pval.Q,mymod_mf$pval.Q,mymod_af$pval.Q, NA,NA)]

MRTab_PCSK9SNPs = copy(dummy)

MRTab_PCSK9SNPs2 = copy(MRTab_PCSK9SNPs)
MRTab_PCSK9SNPs2 = MRTab_PCSK9SNPs2[1:5,]

#' # Compare sexes ####
#' ***
interactionTest(mean1 = MRTab_PCSK9SNPs2$beta_IV[4],se1 = MRTab_PCSK9SNPs2$SE_IV[4],
                mean2 = MRTab_PCSK9SNPs2$beta_IV[1],se2 = MRTab_PCSK9SNPs2$SE_IV[1])
interactionTest(mean1 = MRTab_PCSK9SNPs2$beta_IV[5],se1 = MRTab_PCSK9SNPs2$SE_IV[5],
                mean2 = MRTab_PCSK9SNPs2$beta_IV[2],se2 = MRTab_PCSK9SNPs2$SE_IV[2])


#' # Make Plots ####
#' ***
#' male/female
#' 
#' adjusted/free
#' 
#' all/PCSK9 only

plotdata1 = copy(myMRTab)
plotdata1 = plotdata1[!grepl("treated",phenotype),]
plotdata1[beta.PCSK9<0,beta.LDLC := (-1)*beta.LDLC]
plotdata1[beta.PCSK9<0,beta.PCSK9 := (-1)*beta.PCSK9]

linedata1 = copy(MRTab_PCSK9SNPs)
linedata1 = linedata1[!grepl("treated",phenotype),]
linedata1[,int := 0]
linedata1[,upper := beta_IV+1.96*SE_IV]
linedata1[,lower := beta_IV-1.96*SE_IV]

myPlot1 = ggplot(plotdata1, aes(x=beta.PCSK9, y=beta.LDLC)) +
  facet_wrap(~phenotype, scales = "free",ncol = 2) +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1)+
  #geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", size=1)+
  geom_errorbarh(data = plotdata1,
                 aes(xmin = beta.PCSK9- 1.96*SE.PCSK9, xmax = beta.PCSK9 + 1.96*SE.PCSK9)) +
  geom_errorbar(data = plotdata1,
                aes(ymin = beta.LDLC- 1.96*SE.LDLC, ymax = beta.LDLC+ 1.96*SE.LDLC)) +
  geom_point(size=2.5,color="steelblue")+
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
  labs(x="Effect on PCSK9",
       y = "Effect on LDLC")
myPlot1


tiff(filename = "../figures/SupplementalFigure_MRScatterPlot1_PCSK9SNPs.tiff", 
     width = 3000, height = 2700, res=300, compression = 'lzw')
myPlot1
dev.off()


#' # Save results ####
#' ***
MRTab_PCSK9SNPs[,comment:="PCSK9 SNPs"]

save(myMRTab,file="../results/07_MR_PCSK9_LDLC_SNPwise.RData")
save(MRTab_PCSK9SNPs,file="../results/07_MR_PCSK9_LDLC_summary.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
