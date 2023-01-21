#' ---
#' title: "LD Score Regression Analyses"
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
#' I want to estimate the heritability of PCSK9 in all subsets, and then perform a correlation test within my PCSK9 traits and with lipids. 
#' 
#' In script 08_a, I extracted the GWAS summary statistics and saved them in the format necessary for LDSC. Then, I generated all necessary LDSC commands, which were executed in a putty shell, not in R (reticulate connection with python sometimes unstable). 
#' 
#' In thus R script (08_b), I read in the LDSC log files for heritability and genetic correlation and generate some plots. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_forostar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

load("../results/04_PearsonsCorrelation.RData")

#' # Step 4: Get heritability results ####
#' ***

heritab_log = list.files(path = "../temp/04_LDSC/",pattern = "heritab.log")

myToDoTable = data.table(pheno = gsub("_heritab.log","",heritab_log),
                         file = heritab_log)

heritab = foreach(i = 1:dim(myToDoTable)[1])%do%{
  #i=1
  outfn = myToDoTable[i,file]
  tab = readLines(paste0("../temp/04_LDSC/",outfn))
  filt1 <- grep(pattern = "Total Observed scale h2", tab)
  
  myH2 = tab[filt1]
  myH2 = unlist(strsplit(x=myH2,split = ": "))[2]
  myH2_se = unlist(strsplit(x=myH2,split = "[(]"))[2]
  myH2_se = gsub("[)]","",myH2_se)
  myH2 = unlist(strsplit(x=myH2,split = " [(]"))[1]
  
  myGC = tab[filt1 + 1]
  myGC = unlist(strsplit(x=myGC,split = ": "))[2]
  
  myChi2 = tab[filt1 + 2]
  myChi2 = unlist(strsplit(x=myChi2,split = ": "))[2]
  
  myIntercept = tab[filt1 + 3]
  myIntercept = unlist(strsplit(x=myIntercept,split = ": "))[2]
  myIntercept_se = unlist(strsplit(x=myIntercept,split = "[(]"))[2]
  myIntercept_se = gsub("[)]","",myIntercept_se)
  myIntercept = unlist(strsplit(x=myIntercept,split = " [(]"))[1]
  
  myRatio = tab[filt1 + 4]
  myRatio = unlist(strsplit(x=myRatio,split = ": "))[2]
  myRatio_se = unlist(strsplit(x=myRatio,split = "[(]"))[2]
  myRatio_se = gsub("[)]","",myRatio_se)
  myRatio = unlist(strsplit(x=myRatio,split = " [(]"))[1]
  
  res = data.table::data.table(trait = myToDoTable[i,pheno],
                               h2 = as.numeric(myH2),
                               h2_se = as.numeric(myH2_se),
                               LambdaGC = as.numeric(myGC),
                               meanChi2 = as.numeric(myChi2),
                               Intercept = as.numeric(myIntercept),
                               Intercept_se = as.numeric(myIntercept_se),
                               Ratio = as.numeric(myRatio),
                               Ratio_se = as.numeric(myRatio_se))
  res[,h2_p :=  round(2*pnorm(-abs(h2/h2_se)),digits=5)]
  res[,Intercept_p :=  round(2*pnorm(-abs((Intercept-1)/Intercept_se)),digits=5)]
  res[,Ratio_p :=  round(2*pnorm(-abs(Ratio/Ratio_se)),digits=5)]
  res
}

heritab = rbindlist(heritab)
heritab
save(heritab, file = "../results/08_LDSC_Heritability.RData")

heritab1 = copy(heritab)
heritab1 = heritab1[grepl("PCSK9",trait)]
heritab1[,sex := gsub("PCSK9_","",trait)]
heritab1[,sex := gsub("_.*","",sex)]
heritab1[sex %nin% c("females","males"),sex := "combined"]

heritab1[,statin := gsub("PCSK9_","",trait)]
heritab1[,statin := gsub(".*_","",statin)]
heritab1[statin %nin% c("free","treated"),statin := "combined"]

heritab1[,lowerbound := h2-1.96*h2_se]
heritab1[,upperbound := h2+1.96*h2_se]
heritab1[lowerbound<0,lowerbound:=0]
heritab1[upperbound>1,upperbound:=1]

setorder(heritab1,statin)
p1<- ggplot(heritab1, aes(x=statin, y=h2, fill=sex)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lowerbound, ymax=upperbound), width=.2,
                position=position_dodge(.9))
p1

#' # Step 5: Get genetic correlation results ####
#' ***
#' 
#' 
gencor_log = list.files(path = "../temp/04_LDSC/",pattern = "_gencor.log")

myToDoTable = data.table(pheno = gsub("_gencor.log","",gencor_log),
                         file = gencor_log)

gencortab = foreach(i=1:dim(myToDoTable)[1])%do%{
  #i=1
  outfn = myToDoTable[i,file]
  tab = readLines(paste0("../temp/04_LDSC/",outfn))
  
  filt2 <- grep(pattern = "p1 ", tab)
  dummy2 = tab[filt2+1]
  dummy2 = unlist(strsplit(x=dummy2,split = " "))
  dummy2 = dummy2[dummy2 !=""]
  
  res = data.table::data.table(trait1 = dummy2[1],
                               trait2 = dummy2[2],
                               rg = as.numeric(dummy2[3]),
                               rg_se = as.numeric(dummy2[4]),
                               rg_zScore = as.numeric(dummy2[5]),
                               rg_pval = as.numeric(dummy2[6]))
  res[,trait1 := gsub(projectpath_main,"",trait1)]
  res[,trait1 := gsub("temp/04_LDSC/LDSC_munge_","",trait1)]
  res[,trait1 := gsub(".sumstats.gz","",trait1)]
  res[,trait2 := gsub(projectpath_main,"",trait2)]
  res[,trait2 := gsub("temp/04_LDSC/LDSC_munge_","",trait2)]
  res[,trait2 := gsub(".sumstats.gz","",trait2)]
  
  res                         
}
gencortab = rbindlist(gencortab)
gencortab

gencortab1 = copy(gencortab)
gencortab1 = gencortab1[!grepl("PCSK9",trait2)]

gencortab2 = copy(gencortab)
gencortab2 = gencortab2[grepl("PCSK9",trait2)]

gencortab1[,lipid := gsub("_FEMALE","",trait2)]
gencortab1[,lipid := gsub("_MALE","",lipid)]
gencortab1[,sex := gsub(".*_","",trait2)]
gencortab1[,category := gsub(".*_","",trait1)]
gencortab1[category %in% c("females","males"),category := "combined"]

gencortab1[,lowerbound := rg-1.96*rg_se]
gencortab1[,upperbound := rg+1.96*rg_se]
gencortab1[lowerbound<0,lowerbound:=0]

setorder(gencortab1,trait2)
p2<- ggplot(gencortab1, aes(x=lipid, y=rg, fill=sex)) + 
  facet_wrap(~category) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lowerbound, ymax=upperbound), width=.2,
                position=position_dodge(.9))
p2

save(gencortab, file = "../results/08_LDSC_GeneticCorrelation.RData")

#' # Step 6: Combine plots ####
#' ***
#' 
plotdata = data.table(myX = c(heritab1$statin,gencortab1$lipid),
                      myY = c(heritab1$h2,gencortab1$rg),
                      myFill = c(heritab1$sex,gencortab1$sex),
                      myYmin = c(heritab1$lowerbound,gencortab1$lowerbound),
                      myYmax = c(heritab1$upperbound,gencortab1$upperbound),
                      myCategory = c(rep("heritab",8),gencortab1$category))
plotdata[myCategory == "heritab", myCategory := "A) Heritability (only PCSK9)"]
plotdata[myCategory == "combined", myCategory := "B) Genetic Correlation (combined)"]
plotdata[myCategory == "free", myCategory := "C) Genetic Correlation (free)"]
plotdata[myCategory == "treated", myCategory := "D) Genetic Correlation (treated)"]

plotdata[myFill == "FEMALE",myFill :="females"]
plotdata[myFill == "MALE",myFill :="males"]
#plotdata[myFill == "females",myFill :="a -females"]

p3<- ggplot(plotdata, aes(x=myX, y=myY, fill=myFill)) + 
  facet_wrap(~myCategory,scales = "free_x") + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=myYmin, ymax=myYmax), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(values=c("#82B446","#B2182B","#2166AC"))+
  labs(x="",
       y = "Value",
       fill = "Sex") +
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12,face="bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.1, 0.85),
        legend.background = element_rect(linewidth = 0.5, 
                                         linetype="solid",
                                         colour ="black"))
p3

tiff(filename = "../figures/SupplementalFigure_LDSCresults.tiff", 
     width = 3000, height = 2250, res=300, compression = 'lzw')
p3
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
