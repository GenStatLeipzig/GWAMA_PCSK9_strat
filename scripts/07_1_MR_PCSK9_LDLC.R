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
load("../results/03_GCTA_COJO_filtered.RData")
IndepSignals = IndepSignals_filtered[!duplicated(SNP),]
myTab = copy(IndepSignals)
myTab = myTab[candidateGene == "PCSK9"]

load("../temp/05_OtherGWASs.RData")
myOtherGWAS = myOtherGWAS[grepl("LDL",phenotype)]
myOtherGWAS = myOtherGWAS[candidateGene == "PCSK9"]
myOtherGWAS = myOtherGWAS[bp_hg19 %in% myTab$bp]

load("../temp/04_IATest_input.RData")
data_GWAS = result.2[markername %in% myTab$SNP,]

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
myMRTab = myMRTab[,c(16,1,2,3,4,6,10,11,12,8)]
names(myMRTab)
names(myMRTab) = c("phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")

myMRTab[,EAF.LDLC := myLipids$EAF]
myMRTab[,beta.LDLC := myLipids$beta]
myMRTab[,SE.LDLC := myLipids$SE]
myMRTab[,P.LDLC := myLipids$pval]
myMRTab[,N.LDLC := myLipids$nSamples]

myMRTab

#' # Get Ratio ####
#' ***

test_PCSK9_LDL = myMRTab[,MRfunction_jp(betaX = beta.PCSK9,
                                        seX = SE.PCSK9,
                                        betaY = beta.LDLC,
                                        seY = SE.LDLC)]
myMRTab[,beta.Ratio := test_PCSK9_LDL$beta_IV]
myMRTab[,SEst.Ratio := test_PCSK9_LDL$se_IV1]
myMRTab[,SEnd.Ratio := test_PCSK9_LDL$se_IV2]
myMRTab[,Pst.Ratio := test_PCSK9_LDL$p_IV1]
myMRTab[,Pnd.Ratio := test_PCSK9_LDL$p_IV2]

#' # Get meta per phenotype ####
#' ***
phen = unique(myMRTab$phenotype)

mymod1<-myMRTab[phenotype == phen[1],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod2<-myMRTab[phenotype == phen[2],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod3<-myMRTab[phenotype == phen[3],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod4<-myMRTab[phenotype == phen[4],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod5<-myMRTab[phenotype == phen[5],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod6<-myMRTab[phenotype == phen[6],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod7<-myMRTab[phenotype == phen[7],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]
mymod8<-myMRTab[phenotype == phen[8],
                metagen(TE = beta.Ratio,seTE = SEnd.Ratio,studlab =  SNP)]

dummy = myMRTab[,.N,phenotype]
dummy[,outcome :="LDL-C"]
dummy[,beta_IV := c(mymod1$TE.fixed,mymod2$TE.fixed,mymod3$TE.fixed,mymod4$TE.fixed,
                    mymod5$TE.fixed,mymod6$TE.fixed,mymod7$TE.fixed,mymod8$TE.fixed)]
dummy[,SE_IV := c(mymod1$seTE.fixed,mymod2$seTE.fixed,mymod3$seTE.fixed,mymod4$seTE.fixed,
                  mymod5$seTE.fixed,mymod6$seTE.fixed,mymod7$seTE.fixed,mymod8$seTE.fixed)]
dummy[,pval_IV := c(mymod1$pval.fixed,mymod2$pval.fixed,mymod3$pval.fixed,mymod4$pval.fixed,
                    mymod5$pval.fixed,mymod6$pval.fixed,mymod7$pval.fixed,mymod8$pval.fixed)]
dummy[,Q_IV := c(mymod1$Q,mymod2$Q,mymod3$Q,mymod4$Q,
                 mymod5$Q,mymod6$Q,mymod7$Q,mymod8$Q)]
dummy[,pvalQ_IV := c(mymod1$pval.Q,mymod2$pval.Q,mymod3$pval.Q,mymod4$pval.Q,
                     mymod5$pval.Q,mymod6$pval.Q,mymod7$pval.Q,mymod8$pval.Q)]

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
