#' ---
#' title: "Supplemental Figures: Forest Plots (Caterpillar)"
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
#' Forest plots for PCSK9 locus (all four SNPs)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../temp/04_IATest_input.RData")
myTab = copy(result.2)
#myTab = myTab[candidateGene == "PCSK9",]

myTab[grepl("rs693",markername),beta:= beta*(-1)]

myTab[,lowerCI95 := beta-1.96*SE]
myTab[,upperCI95 := beta+1.96*SE]
myTab[,dumID := gsub("PCSK9_","",phenotype)]
myTab[,dumID := gsub("_"," ",dumID)]
setorder(myTab, beta)

tag = format(Sys.time(), "%Y-%m-%d")
tag2 = gsub("2023-","23-",tag)
tag2 = gsub("-","",tag2)

mySNPs = unique(myTab$rsID)

for(i in 1:length(mySNPs)){
  #i=1
  message("Working on phenotype ",mySNPs[i])
  myTab_filt = copy(myTab)
  myTab_filt = myTab_filt[rsID == mySNPs[i]]
  myGene = unique(myTab_filt[,candidateGene])
  if(myGene == "NOS1") myGene = "KSR2"
  
  filename1 = paste0("../figures/SupplementalFigure5_ForestPlots_",myGene,"_",mySNPs[i],"_",tag2,".tiff")
  xmin = min(c(0,myTab_filt$lowerCI95, myTab_filt$upperCI95), na.rm = T)
  xmax = max(c(0,myTab_filt$lowerCI95, myTab_filt$upperCI95), na.rm = T)
  xrange = xmax - xmin
  xmin_margin = xmin - 0.5*xrange
  xmax_margin = xmax + 0.5*xrange
  
  tiff(filename = filename1,
       width = 2400, height = 1200, res=250, compression = 'lzw')
  par(mar=c(5,6,0,4))
  par(font=1)
  dets = forest(x=myTab_filt[,beta], 
                sei = myTab_filt[,SE],
                xlab=paste0(myGene," - ",mySNPs[i]),
                showweights=F, 
                slab = myTab_filt[,dumID],
                ylim = c(-0, nrow(myTab_filt)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin, xmax))
  par(font=4)
  text(min(dets$xlim), max(dets$ylim-1.5), "Subgroups", pos=4)
  text(max(dets$xlim), max(dets$ylim-1.5), "             beta estimate (95% CI)", pos=2)
  
  dev.off()
  
  
  
}

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
