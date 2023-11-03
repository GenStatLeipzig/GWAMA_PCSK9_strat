#' ---
#' title: "Main Figure 2: 2-way interaction plots"
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
#' **Figure 2: Plot of differences**
#' 
#' I want some sort of panel plot (facets)
#' 
#' - A) 3-way interaction (x=b_(F,treated)-b_(M,treated), y=b_(F,free)-b_(M,free))
#' - B) 2-way sex interaction (x=b_M, y =b_F)
#' - C) 2-way statin interaction (x=b_treated, y=b_free)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load and prep data ####
#' ***
load("../results/04_IATest_2way_complete.RData")

IATab[,myX := trait1_beta]
IATab[,myY := trait2_beta]
IATab[,myX_SE := trait1_SE]
IATab[,myY_SE := trait2_SE]
IATab[,rsID := gsub(":.*","",markername)]
IATab[,myLabel := "x - not sig"]
IATab[IA_hierarch_fdr5proz == T,myLabel := paste0(candidateGene," \n",rsID," \n")]


myPlot1 = ggplot(IATab, aes(x=myX, y=myY, color=myLabel, alpha = IA_hierarch_fdr5proz)) +
  facet_wrap(~type, nrow = 1, 
             strip.position = "left", 
             labeller = as_labeller(c("sexIA" = "SNP effect in women", 
                                      "statinIA" = "SNP effect in statin-free individuals") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = myY- 1.96*myY_SE, ymax = myY+ 1.96*myY_SE),linewidth=1.1) +
  geom_errorbarh(aes(xmin = myX- 1.96*myX_SE, xmax = myX + 1.96*myX_SE),linewidth=1.1) +
  theme_bw(base_size = 10) + 
  scale_colour_manual(values=c(brewer.pal(10, "Paired"),"#000000"))+
  scale_alpha_manual(values = c(0.3,1),labels = c("not significant", "significant"))+
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
  guides(alpha = NULL, size = NULL)+
  labs(x="SNP effect in men / statin-treated individuals", 
       color="", alpha="") 

myPlot1

#' # Save plots ####
#' ***
tiff(filename = "../figures/MainFigure2_InteractionScatterPlot_231101.tiff",
     width = 2250, height = 1125, res=240, compression = 'lzw')
myPlot1
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))

