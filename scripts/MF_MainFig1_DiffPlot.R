#' ---
#' title: "Main Figure 1: 3-way interaction plot"
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
#' **Figure 1: Plot of differences**
#' 
#' Test, I want to check if this results in a meaningful plot!
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load overview data ####
#' ***
load("../results/03_InteractionTests_3way.RData")
table(duplicated(IATab_3way$markername))

#' IA_diff := (trait2_beta - trait1_beta) - (trait4_beta - trait3_beta)
#' 
IATab_3way[,diff1 := trait2_beta - trait1_beta]
IATab_3way[,diff2 := trait4_beta - trait3_beta]
IATab_3way[,gene := ""]
IATab_3way[IA_pval<0.05,gene := candidateGene]

myPlot = ggplot(IATab_3way, aes(x=diff1, y=diff2, color=gene,size=gene)) +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", size=1.15)+
  geom_point()+ 
  theme_bw(base_size = 10) + 
  scale_colour_manual(values=c("#000000","#B2182B","#2166AC","#82B446","#DDA0DD","#7846B4"),
                      labels=c("no \ninteraction","KHDRBS2","PCSK9","PLB1",
                               "SASH1","SLCO1B1"))+
  scale_size_manual(values=c(2,rep(3.5,5)))+
  guides(size="none")+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  labs(x="(females - males) in statin-treated individuals", 
       y = "(females - males) in statin-free individuals",
       color="Candidate genes \nwith 3-way \ninteraction")


myPlot


#' # Save plot ####
#' ***
tiff(filename = "../figures/MainFigure1_3wayIA.tiff",
     width = 2400, height = 1800, res=250, compression = 'lzw')
myPlot
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))

