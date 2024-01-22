#' ---
#' title: "Main Figure 3: Forest plot of MR estimates"
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
#' Please note: this script ran on my BSU laptop, as I failed to install the forestploter package on the IMISE servers and WANDERTAUBE. 
#' 
#' **Main Figure 3: Forest Plot of MR estimates**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

library(grid)
library(forestploter)
library(data.table)

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2023-","23-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
load("../results/06_MR_IVWEstimates_UKBB_230907.RData")
head(myMRTab_results)

#' # Prep data ####
#' ***
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
myTab[,dumID := gsub("females","Women",dumID)]
myTab[,dumID := gsub("males","Men",dumID)]
myTab[,dumID := gsub(" free"," - free",dumID)]
myTab[,dumID := gsub(" treated"," - treated",dumID)]

dummy = data.table(dumID = unique(myTab$dumID))
dummy[,dumID := gsub(" [(]IVW[)]","",dumID)]
myTab = rbind(myTab,dummy, fill=T)
myTab[,subgroup := setting]
myTab[is.na(setting),subgroup := dumID]

myTab = myTab[c(23,7,15,24,8,16,19,3,11,18,2,10,
                17,1,9,21,5,13,20,4,12,22,6,14),]
myTab = myTab[,c(12,4,9,10,5)]
setDF(myTab)
myTab$subgroup <- ifelse(is.na(myTab$beta_IVW), 
                         myTab$subgroup,
                         paste0("   ", myTab$subgroup))
myTab$` ` <- paste(rep(" ", 20), collapse = " ")
myTab$`Estimate [95% CI]` <- ifelse(is.na(myTab$SE_IVW), "",
                                    sprintf("%.2f [%.2f, %.2f]",
                                            myTab$beta_IVW, myTab$lowerCI95, myTab$upperCI95))
head(myTab)

setDT(myTab)
setnames(myTab,"subgroup","Subgroup")
myTab[4,Subgroup := "Treated"] 
myTab[19,Subgroup := "Free"] 

#' # Plotting ####
#' ***
#' ## Plot 1 ####
#' Plot with subgroups as extra lines
#' 

tm1<- forest_theme(core=list(bg_params=list(fill = c("lightgrey","white","white",
                                                     "lightgrey","white","white",
                                                     "lightgrey","white","white",
                                                     "lightgrey","white","white",
                                                     "lightgrey","white","white",
                                                     "lightgrey","white","white",
                                                     "lightgrey","white","white",
                                                     "lightgrey","white","white"))))
p1<- forest(myTab[,c(1,6,7)],
            est = myTab$beta_IVW,
            lower = myTab$lowerCI95, 
            upper = myTab$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            xlim = c(0, 1.6),
            xlab = "       Causal estimate for effect of PCSK9 on LDL-C",
            theme = tm1)

plot(p1)

filename1 = paste0("../figures/MainFigure3_MRForerstPlots_v1_",tag,".tiff")
tiff(filename = filename1,
     width = 1400, height = 1700, res=250, compression = 'lzw')
plot(p1)
dev.off()

#' ## Plot 2 ####
#' Alternative: no free rows for the subgroup, but different coloring per subgroup (or alternating background per subgroup)

dummy = myTab[c(1,4,7,10,13,16,19,22), Subgroup]
dummy = rep(dummy,each=2)
myTab2 = copy(myTab)
myTab2 = myTab2[-c(1,4,7,10,13,16,19,22),]
myTab2[,Subgroup := gsub(" ","",Subgroup)]
myTab2[,Subgroup := paste0(dummy," (",Subgroup,")")]

tm2<- forest_theme(core=list(bg_params=list(fill = c("lightgrey","lightgrey",
                                                     "white","white",
                                                     "lightgrey","lightgrey",
                                                     "white","white",
                                                     "lightgrey","lightgrey",
                                                     "white","white",
                                                     "lightgrey","lightgrey",
                                                     "white","white"))))
p2<- forest(myTab2[,c(1,6,7)],
            est = myTab2$beta_IVW,
            lower = myTab2$lowerCI95, 
            upper = myTab2$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            xlim = c(0, 1.6),
            xlab = "       Causal estimate for effect of PCSK9 on LDL-C",
            theme = tm2)

plot(p2)

filename2 = paste0("../figures/MainFigure3_MRForerstPlots_v2_",tag,".tiff")
tiff(filename = filename2,
     width = 1400, height = 1200, res=250, compression = 'lzw')
plot(p2)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
