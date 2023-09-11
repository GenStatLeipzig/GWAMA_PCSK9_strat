#' ---
#' title: "Cis-MR: PCSK9 --> LDL-C, leave-one-out, more SNPs"
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
#' I test this five times in all subgroups with valid SNPs using the following SNP sets: all SNPs with suggestive significance 
#' 
#' **Using UKBB for outcome data**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_forostar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))
source("../helperFunctions/TwoWayInteractionTest_jp.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag2 = gsub("2023-","23-",tag)
tag2 = gsub("-","",tag2)

#' # Get data ####
#' ***
#' ## Load ####
load("../temp/05_GWAS_allLoci.RData")
data_GWAS = data_GWAS[candidateGene == "PCSK9"]
data_GWAS[,rsID := gsub(":.*","",markername)]
data_GWAS = data_GWAS[markername != "rs149093957:55522741:C:CCTGT",]
loaded = load("../temp/06_LDLC_CAD_assoc_stratified_moreSNPs.RData")
dat2 = get(loaded)
dat2 = dat2[grepl("LDLC",pheno),]

data_GWAS = data_GWAS[rsID %in% dat2$ID,]
myPhenos = data_GWAS[,unique(phenotype)]
data_GWAS = data_GWAS[phenotype %in% c(myPhenos)]

#' # Match ####
#' ***
myLipids = copy(dat2)
myLipids[,dumID := gsub("LDLC_","",pheno)]
myLipids[,dumID := paste(ID,dumID,sep="_")]

data_GWAS[,dumID := gsub("PCSK9_","",phenotype)]
data_GWAS[,dumID := paste(rsID,dumID,sep="_")]

matched = match(data_GWAS$dumID,myLipids$dumID)
table(is.na(matched))
myLipids = myLipids[matched]
table(myLipids$dumID == data_GWAS$dumID)

#' # Check allele ####
#' ***
plot(data_GWAS$EAF,myLipids$A1_FREQ)
abline(0,1)

table(myLipids$A1 == data_GWAS$EA)
filt = myLipids$A1 != data_GWAS$EA

EA1 = myLipids$A1
myLipids[filt, A1_FREQ     := 1-A1_FREQ    ]
myLipids[filt, BETA:= BETA     *(-1)]

plot(data_GWAS$EAF,myLipids$A1_FREQ)
abline(0,1)

#' # Merge ####
#' ***
myMRTab = copy(data_GWAS)
myMRTab = myMRTab[,c(17,16,1,2,3,4,6,10,11,12,8)]
names(myMRTab)
names(myMRTab) = c("gene","phenotype","SNP","chr","pos","EA","EAF.PCSK9","beta.PCSK9","SE.PCSK9","P.PCSK9","N.PCSK9")
myMRTab[, FStat.PCSK9 := (beta.PCSK9/SE.PCSK9)^2]
myMRTab[,.N,by=phenotype]
myMRTab[FStat.PCSK9<10,.N,by=phenotype]

myMRTab[,EAF.LDLC := myLipids$A1_FREQ]
myMRTab[,beta.LDLC := myLipids$BETA]
myMRTab[,SE.LDLC := myLipids$SE]
myMRTab[,P.LDLC := myLipids$P]
myMRTab[,N.LDLC := myLipids$OBS_CT]

myMRTab[, rsID := gsub(":.*","",SNP)]

save(myMRTab,file = paste0("../temp/06_MR_input_UKBB_moreSNPs",tag2,".RData"))
load(paste0("../temp/06_MR_input_UKBB_moreSNPs",tag2,".RData"))

#' # Wald estimates ####
#' ***
myMRTab[,beta.Wald := beta.LDLC/beta.PCSK9]
myMRTab[,SE.Wald := SE.LDLC/sqrt(beta.PCSK9^2)]
myMRTab[,P.Wald := 2*pnorm(-abs(beta.Wald/SE.Wald))]

myMRTab[,table(P.Wald<0.05)]
myMRTab[P.Wald>0.05, c(2,18,8,10,12,19:21)]
myMRTab[FStat.PCSK9<10, c(2,18,8,10,12,19:21)]

save(myMRTab,file = paste0("../results/06_MR_WaldEstimates_UKBB_moreSNPs_",tag2,".RData"))

#' # IVW estimates ####
#' ***
dumTab1 = foreach(i = 1:length(myPhenos))%do%{
  #i=1
  myPheno = myPhenos[i]
  message("Working on phenotype ",myPheno)
  myMRTab1 = copy(myMRTab)
  myMRTab1 = myMRTab1[phenotype == myPheno]
  # myPheno = gsub("PCSK9_","",myPheno)
  # myPheno = gsub("_"," ",myPheno)
  
  dumTab2 = foreach(j = 1:1)%do%{
    #j=1
    if(j==1){
      setting = "all (fixed)"
      mrob = mr_input(bx = myMRTab1$beta.PCSK9, bxse = myMRTab1$SE.PCSK9,
                      by = myMRTab1$beta.LDLC, byse = myMRTab1$SE.LDLC,
                      snps = myMRTab1$rsID, exposure = myPheno, outcome = "LDLC")
    }else{
      k=j-1
      setting = paste0("w/o ",myMRTab1$rsID[k])
      mrob = mr_input(bx = myMRTab1$beta.PCSK9[-k], bxse = myMRTab1$SE.PCSK9[-k],
                       by = myMRTab1$beta.LDLC[-k], byse = myMRTab1$SE.LDLC[-k],
                       snps = myMRTab1$rsID[-k], exposure = myPheno, outcome = "LDLC")
    }
    mod1 = mr_ivw(mrob,model="fixed")
    mod2 = mr_egger(mrob)
    
    res = data.table(phenotype = mod1@Exposure,
                     setting = setting,
                     outcome = c(mod1@Outcome),
                     beta_IVW = c(mod1@Estimate),
                     SE_IVW = c(mod1@StdError),
                     pval_IVW = c(mod1@Pvalue),
                     HeteroStat_IVW = c(mod1@Heter.Stat[1]),
                     HeteroStat_pval_IVW = c(mod1@Heter.Stat[2]),
                     beta_egger = c(mod2@Estimate),
                     SE_egger = c(mod2@StdError.Est),
                     pval_egger = c(mod2@Pvalue.Est),
                     beta_egger_int = c(mod2@Intercept),
                     SE_egger_int = c(mod2@StdError.Int),
                     pval_egger_int = c(mod2@Pvalue.Int),
                     HeteroStat_egger = c(mod2@Heter.Stat[1]),
                     HeteroStat_pval_egger = c(mod2@Heter.Stat[2]))
    res
  } 
  dumTab2 = rbindlist(dumTab2)
  dumTab2

}
myMRTab_results = rbindlist(dumTab1)

#' Check Egger intercept
myMRTab_results[pval_egger_int<0.05,]

#' When removing the lead SNP some pleitropic effects occur!
myMRTab_results[,plot(beta_IVW,beta_egger)]
abline(0,1)
myMRTab_results[pval_egger_int>0.05,plot(beta_IVW,beta_egger)]
abline(0,1)
myMRTab_results[setting == "all (fixed)",plot(beta_IVW,beta_egger)]
abline(0,1)
# myMRTab_results[setting == "w/o rs11583680",plot(beta_IVW,beta_egger)]
# abline(0,1)

save(myMRTab_results,file = paste0("../results/06_MR_IVWEstimates_UKBB_moreSNPs_",tag2,".RData"))

#' # Forest Plots ####
#' ***
#' ## For setting all
myTab = copy(myMRTab_results)
myTab = myTab[grepl("all",setting)]
myTab = rbind(myTab,myTab)
myTab[,setting := rep(c("IVW","Egger"),each=8)]
myTab[setting == "Egger", beta_IVW := beta_egger]
myTab[setting == "Egger", SE_IVW := SE_egger]
myTab[setting == "Egger", pval_IVW := pval_egger]
myTab[setting == "Egger", HeteroStat_IVW := HeteroStat_egger]
myTab[setting == "Egger", HeteroStat_pval_IVW := HeteroStat_pval_egger]
myTab = myTab[,c(1:8)]

myTab[,lowerCI95 := beta_IVW-1.96*SE_IVW]
myTab[,upperCI95 := beta_IVW+1.96*SE_IVW]
myTab[,dumID := paste0(phenotype," (",setting,")")]
myTab[,dumID := gsub("PCSK9_","",dumID)]
myTab[,dumID := gsub("_"," ",dumID)]

myTab[,rank := c(6,4,3,7,5,8,1,2,6,4,3,7,5,8,1,2)]
setorder(myTab, rank,-setting)

filename1 = paste0("../figures/MRPlots/06_MR_ForerstPlots_acrossPhenos_inclEgger_UKBB_moreSNPs_",tag2,".tiff")
xmin = min(c(myTab$lowerCI95, myTab$upperCI95), na.rm = T)
xmax = max(c(myTab$lowerCI95, myTab$upperCI95), na.rm = T)
xrange = xmax - xmin
xmin_margin = xmin - 0.65*xrange
xmax_margin = xmax + 0.5*xrange

tiff(filename = filename1,
     width = 2400, height = 1200, res=250, compression = 'lzw')
par(mar=c(5,6,0,4))
par(font=1)
dets = forest(x=myTab[,beta_IVW],refline = NA, 
              sei = myTab[,SE_IVW],
              xlab="Causal estimate for effect of PCSK9 on LDLC",
              showweights=F, 
              slab = myTab[,dumID],
              ylim = c(-0, nrow(myTab)+3) , 
              cex =1, 
              xlim = c(xmin_margin, xmax_margin), 
              alim = c(xmin, xmax))
par(font=4)
text(min(dets$xlim), max(dets$ylim-1.5), "Subgroup (method)", pos=4)
text(max(dets$xlim), max(dets$ylim-1.5), "Estimate (95% CI)", pos=2)
par(font=2)
dev.off()

filename2 = paste0("../figures/MRPlots/06_MR_ForerstPlots_acrossPhenos_UKBB_moreSNPs_",tag2,".tiff")
myTab = myTab[grepl("IVW",setting),]
myTab[,dumID := gsub(" [(].*","",dumID)]
xmin = min(c(myTab$lowerCI95, myTab$upperCI95), na.rm = T)
xmax = max(c(myTab$lowerCI95, myTab$upperCI95), na.rm = T)
xrange = xmax - xmin
xmin_margin = xmin - 0.65*xrange
xmax_margin = xmax + 0.5*xrange

tiff(filename = filename2,
     width = 2400, height = 1200, res=250, compression = 'lzw')
par(mar=c(5,6,0,4))
par(font=1)
dets = forest(x=myTab[,beta_IVW],refline = NA, 
              sei = myTab[,SE_IVW],
              xlab="Causal estimate for effect of PCSK9 on LDLC",
              showweights=F, 
              slab = myTab[,dumID],
              ylim = c(-0, nrow(myTab)+3) , 
              cex =1, 
              xlim = c(xmin_margin, xmax_margin), 
              alim = c(xmin, xmax))
par(font=4)
text(min(dets$xlim), max(dets$ylim-1.5), "Subgroup", pos=4)
text(max(dets$xlim), max(dets$ylim-1.5), "IVW estimate (95% CI)", pos=2)
par(font=2)
dev.off()

#' ## Per phenotype
myTab = copy(myMRTab)
myTab[,lowerCI95 := beta.Wald-1.96*SE.Wald]
myTab[,upperCI95 := beta.Wald+1.96*SE.Wald]
myTab[,dumID := paste(phenotype,gsub(":.*","",SNP),sep=" - ")]
myTab[,dumID := gsub("PCSK9_","",dumID)]
setorder(myTab, phenotype,beta.Wald)
myTab2 = copy(myMRTab_results)
myTab2 = myTab2[!grepl("random",setting)]

for(i in 1:8){
  #i=1
  message("Working on phenotype ",myPhenos[i])
  filename1 = paste0("../figures/MRPlots/06_MR_ForerstPlots",myPhenos[i],"_UKBB_moreSNPs_",tag2,".tiff")
  myTab_filt = copy(myTab)
  myTab_filt = myTab_filt[phenotype == myPhenos[i]]
  myTab2_filt = copy(myTab2)
  myTab2_filt = myTab2_filt[phenotype == myPhenos[i]]
  myTab_filt[,rsID := gsub(":.*","",SNP)]
  myTab_filt[,rsID := paste("w/o",rsID,sep=" ")]
  myTab2_filt[,rsID := gsub("PCSK9 w/o ","",setting)]
  myTab2_filt[1,rsID := NA]
  # matched = match(myTab2_filt$rsID,myTab_filt$rsID)
  # matched[c(1,2)] = c(-1,0) 
  # matched = matched +2
  # matched = match(myTab_filt$rsID,myTab2_filt$rsID)
  # matched = c(1,matched)
  # myTab2_filt = myTab2_filt[matched,]
  
  xmin = min(c(0,myTab_filt$lowerCI95, myTab_filt$upperCI95), na.rm = T)
  xmax = max(c(0,myTab_filt$lowerCI95, myTab_filt$upperCI95), na.rm = T)
  xrange = xmax - xmin
  xmin_margin = xmin - 0.65*xrange
  xmax_margin = xmax + 0.5*xrange
  
  tiff(filename = filename1,
       width = 2400, height = 1200, res=250, compression = 'lzw')
  par(mar=c(5,6,0,4))
  par(font=1)
  dets = forest(x=myTab_filt[,beta.Wald], 
                sei = myTab_filt[,SE.Wald],
                xlab="Causal estimate for effect of PCSK9 on LDLC",
                showweights=F, 
                slab = myTab_filt[,dumID],
                ylim = c(-2, nrow(myTab_filt)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin, xmax))
  par(font=4)
  text(min(dets$xlim), max(dets$ylim-1.5), "Strata - SNP", pos=4)
  text(max(dets$xlim), max(dets$ylim-1.5), "             Wald estimate (95% CI)", pos=2)
  par(font=2)
  addpoly(x = c(myTab2_filt[,beta_IVW],myTab2_filt[,beta_egger]),
          sei =  c(myTab2_filt[,SE_IVW],myTab2_filt[,SE_egger]),
          row=c(-0.5,-1.5),
          cex=1,
          mlab=c("IVW","Egger"))
  dev.off()
  
  
  
}


#' # Interaction Test ####
#' ***
#' I want to perform all 6 tests although the tests with females treated might have a power problem 
#' 
#' - 2-way interaction for sex (combined and statin free)
#' - 2-way interaction for statin treatment (combined and males)
#' 
ToDoList = data.table(trait1 = c("males", "males_free","males_treated","treated","males_treated","females_treated"),
                      trait2 = c("females", "females_free","females_treated","free","males_free","females_free"),
                      type = c("sexIA","sexIA_free","sexIA_treated","statinIA","statinIA_males","statinIA_females"))

dumTab3 = foreach(i=1:6)%do%{
  #i=1
  IATab = TwoWayInteractionTest_jp(data = myMRTab_results[!grepl("random",setting),],
                                   pheno1 = ToDoList[i,trait1],
                                   pheno2 = ToDoList[i,trait2],
                                   type = ToDoList[i,type],
                                   useBestPheno = F,
                                   corCol = "corSex",
                                   data_type = "MR", suffix = "IVW")
  IATab[,setting := myMRTab_results$setting[1]]
  IATab[,type := ToDoList[i,type]]
  IATab
}
IATab_IVW = rbindlist(dumTab3)
IATab_IVW[setting == "all (fixed)",table(IA_pval<0.05,type)]

dumTab3 = foreach(i=1:6)%do%{
  #i=1
  IATab = TwoWayInteractionTest_jp(data = myMRTab_results[!grepl("random",setting),],
                                   pheno1 = ToDoList[i,trait1],
                                   pheno2 = ToDoList[i,trait2],
                                   type = ToDoList[i,type],
                                   useBestPheno = F,
                                   corCol = "corSex",
                                   data_type = "MR", suffix = "egger")
  IATab[,setting := myMRTab_results$setting[1]]
  IATab[,type := ToDoList[i,type]]
  IATab
}
IATab_egger = rbindlist(dumTab3)
IATab_egger[setting == "all (fixed)",table(IA_pval<0.05,type)]

save(IATab_IVW, file = paste0("../results/06_MR_InteractionTest_UKBB_moreSNPs_",tag2,".RData"))

#' # Interaction Plots ####
#' ***
#' I want a facet plot with the 2-way interaction and setting all and w/o rs11583680 (color = type, shape = setting)
#' 
plotData = copy(IATab_IVW)
#plotData = plotData[setting %in% c("all (fixed)","w/o rs11583680"),]
plotData[,type2 := gsub("_.*","",type)]
plotData[type2 == "sexIA",type2 := "A) sex-interaction"]
plotData[type2 == "statinIA",type2 := "B) statin-interaction"]

plotData[,myX := trait1_beta]
plotData[,myY := trait2_beta]
plotData[,myX_SE := trait1_SE]
plotData[,myY_SE := trait2_SE]

myPlot1 = ggplot(plotData, aes(x=myX, y=myY, color=type)) +
  facet_wrap(~type2, scales = "free", nrow = 1, 
             strip.position = "left", 
             labeller = as_labeller(c("A) sex-interaction" = "Causal effect in females", 
                                      "B) statin-interaction" = "Causal effect in statin-free") ) )+
  #facet_wrap(~type2,scales = "free") +
  #geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  #geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3) + 
  #geom_point(data = IATab, aes(size = abs(IA_diff))) + 
  geom_errorbar(aes(ymin = myY- 1.96*myY_SE, ymax = myY+ 1.96*myY_SE)) +
  geom_errorbarh(aes(xmin = myX- 1.96*myX_SE, xmax = myX + 1.96*myX_SE)) +
  theme_bw(base_size = 10) + 
  scale_colour_manual(values=c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#F0027F","#6A3D9A"),
                      labels=c("statin-combined","statin-free","statin-treated","sex-combined","females","males"))+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  ylab(NULL) +
  labs(x="Causal effect in males                                                      Causal effect in statin-treated")

myPlot1

tiff(filename = paste0("../figures/MRPlots/06_InteractionScatterPlot_UKBB_moreSNPs_",tag2,".tiff"),
     width = 2250, height = 1125, res=200, compression = 'lzw')
myPlot1
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
