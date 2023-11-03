#' ---
#' title: "Main Figure 4: Scatter plot of MR estimates"
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
#' 
#' **Main Figure 4: Interaction Plot of MR estimates**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_forostar.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2023-","23-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
load("../results/06_MR_InteractionTest_UKBB_230907.RData")
head(IATab_IVW)

#' # Prep data ####
#' ***
plotData = copy(IATab_IVW)
plotData = plotData[setting %in% c("all (fixed)","w/o rs11583680"),]
plotData[,type2 := gsub("_.*","",type)]
plotData[type2 == "sexIA",type2 := "A) sex-interaction"]
plotData[type2 == "statinIA",type2 := "B) statin-interaction"]

plotData[,myX := trait1_beta]
plotData[,myY := trait2_beta]
plotData[,myX_SE := trait1_SE]
plotData[,myY_SE := trait2_SE]

#' # Plotting ####
#' ***
myPlot1 = ggplot(plotData, aes(x=myX, y=myY, color=type,shape=setting)) +
  facet_wrap(~type2, nrow = 1, #scales = "free", 
             strip.position = "left", 
             labeller = as_labeller(c("A) sex-interaction" = "Causal effect in women", 
                                      "B) statin-interaction" = "Causal effect in statin-free individuals") ) )+
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
                      labels=c("statin-combined","statin-free","statin-treated","sex-combined","women","men"))+
  scale_shape_manual(values = c(19,17),
                     labels = c("main analysis \nwith 4 SNPs \n","leave-one-out \nwithout rs11583680"))+
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
  labs(x="Causal effect in men / statin-treated individuals")

myPlot1

filename1 = paste0("../figures/MainFigure4_MRInteraction_",tag,".tiff")
tiff(filename = filename1,width = 2250, height = 1125, res=200, compression = 'lzw')
myPlot1
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
