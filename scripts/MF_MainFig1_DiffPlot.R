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
load("../results/03_InteractionTests_2way.RData")
table(duplicated(IATab_2way$markername))
IATab_2way_sex = copy(IATab_2way)
IATab_2way_statin = copy(IATab_2way)
IATab_2way_statin = IATab_2way_statin[type != "sexIA"]
IATab_2way_sex = IATab_2way_sex[type == "sexIA"]
table(duplicated(IATab_2way_sex$markername))
table(duplicated(IATab_2way_statin$markername))

#' # 3-way Interaction Plot ####
#' ***
#' Here I will plot the two differences:
#' 
#' * diff 1: (females - males) in statin-treated individuals = (trait2_beta - trait1_beta) 
#' * diff 2: (females - males) in statin-free individuals = (trait4_beta - trait3_beta)
#' 
#' I force the effect of females treated to be positive. 
#' 
filt = IATab_3way$trait2_beta < 0
table(filt)
IATab_3way[filt,trait1_beta := trait1_beta *(-1)]
IATab_3way[filt,trait2_beta := trait2_beta *(-1)]
IATab_3way[filt,trait3_beta := trait3_beta *(-1)]
IATab_3way[filt,trait4_beta := trait4_beta *(-1)]
IATab_3way[,diff1 := trait2_beta - trait1_beta]
IATab_3way[,diff2 := trait4_beta - trait3_beta]

matched1 = match(IATab_3way$markername,IATab_2way_sex$markername)
IATab_3way[,sexIA := IATab_2way_sex[matched1,IA_hierarch_fdr5proz]]
matched1 = match(IATab_3way$markername,IATab_2way_sex$markername)
IATab_3way[,statinIA := IATab_2way_statin[matched1,IA_hierarch_fdr5proz]]

IATab_3way[,IA_type:="unspecific"]
IATab_3way[sexIA==T,IA_type:="sex-specific"]
IATab_3way[statinIA==T,IA_type:="statin-specific"]
IATab_3way[statinIA==T & sexIA==T,IA_type:="sex-statin-specific"]
table(IATab_3way$IA_type)

IATab_3way[,gene := ""]
IATab_3way[IA_type!="unspecific",gene := candidateGene]

brewer.pal(7,"BrBG")

myPlot1 = ggplot(IATab_3way, aes(x=diff1, y=diff2, color=gene,shape=IA_type,size = gene)) +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", size=1.15)+
  geom_point()+ 
  theme_bw(base_size = 10) + 
  scale_size_manual(values=c(2,rep(3.5,7)))+
  scale_colour_manual(values=c("#000000","#1B9E77","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"),
                      labels=c("no \ninteraction","FOSL2","gene desert","KHDRBS2","MACROD2","PRKAG2",
                               "SASH1","SLCO1B3"))+
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
       color="Candidate genes \nwith interaction",
       shape="Type of interaction ")


myPlot1


#' # 2-way Interaction Plot ####
#' ***
#' Here I will plot the two 2-way interaction scatter plot using facets:
#' 
#' * facet 1: males vs. females (best statin setting)
#' * facet 2: treated vs free (best sex setting)
#' 
#' I force the effect of trait 2 (females or free) to be positive, and label the genes 
#' 
load("../results/05_GCTA_COJO.RData")

plotdata = copy(IATab_2way)
plotdata = plotdata[markername %in% IATab_3way$markername]
plotdata[,gene := ""]
plotdata[IA_hierarch_fdr5proz==T,gene := candidateGene]
filt = plotdata$trait2_beta < 0
table(filt)
plotdata[filt,trait1_beta := trait1_beta *(-1)]
plotdata[filt,trait2_beta := trait2_beta *(-1)]

plotdata[,type3 := "males vs females"]
plotdata[type =="statinIA",type3 := "treated vs free"]
plotdata[,gene2 := ""]
plotdata[IA_hierarch_fdr5proz==T & markername %in% IndepSignals$SNP,gene2 := candidateGene]

myPlot2 = ggplot(plotdata, aes(x=trait1_beta, y=trait2_beta, color=gene,size = gene,label=gene2)) +
  facet_wrap(~type3)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", size=1.15)+
  geom_point()+ 
  theme_bw(base_size = 10) + 
  scale_size_manual(values=c(2,rep(3.5,7)))+
  scale_colour_manual(values=c("#000000","#1B9E77","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"),
                      labels=c("no \ninteraction","FOSL2","gene desert","KHDRBS2","MACROD2","PRKAG2",
                               "SASH1","SLCO1B3"))+
  guides(size="none",color="none")+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  labs(x="effect in trait 1", 
       y = "effect in trait 2")+
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf)

myPlot2

#' # Save plots ####
#' ***
tiff(filename = "../figures/MainFigure1_3wayIA.tiff",
     width = 2400, height = 1800, res=250, compression = 'lzw')
myPlot1
dev.off()

tiff(filename = "../figures/MainFigure1_2wayIA.tiff",
     width = 2400, height = 1200, res=300, compression = 'lzw')
myPlot2
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))

