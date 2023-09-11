#' ---
#' title: "Supplemental Figures: Miami Plots"
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
#' Manhattan Plot colored by stratum association
#' 
#' # Initialize ####
#' ***
time0 = Sys.time()

source("../SourceFile_angmar.R")
source("../helperFunctions/manhattanPlot.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_230.*","",pheno)]

dumTab1 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Working on data ",myRow$pheno)
  data = fread(myRow$statistic_path)
  
  message("    SNPs before filter: ",dim(data)[1])
  data = data[invalidAssocs  == F & pval<=0.01,]
  message("    SNPs after filter: ",dim(data)[1])
  
  data[,chrPosPheno := paste(chr,bp_hg19,phenotype,sep=":")]
  data
}

erg1 = rbindlist(dumTab1)

#' # Add gene names ####
#' ***
loaded1 = load("../results/02_LociOverallPhenotypes_filtered.RData")
loaded1
result.5[candidateGene == "NOS1",candidateGene := "KSR2"]
result.5[candidateGene == "HP",candidateGene := "HP/HPR"]
result.5[candidateGene == "SLCO1B3",candidateGene := "SLCO1B1"]
loaded3 = load("../results/02_LociPerPhenotype.RData")
loaded3
result.3

dumTab2 = foreach(i = 1:dim(result.5)[1])%do%{
  #i=1
  myRow = copy(result.5)
  myRow = myRow[i,]
  
  data = copy(result.3)
  data = data[chr==myRow$chr & bp_hg19<myRow$region_end & bp_hg19>myRow$region_start,]
  data[,candidateGene := myRow$candidateGene]
  data
}
result.4 = rbindlist(dumTab2)
result.4 = result.4[chr %in% result.5$chr]

result.4[,dumID := paste(phenotype,markername,sep="::")]
erg1[,dumID := paste(phenotype,markername,sep="::")]

matched1 = match(erg1$dumID,result.4$dumID)
table(is.na(matched1))
erg1[,candidateGene :=result.4[matched1,candidateGene]]

erg1[, phenotype := gsub("PCSK9_","",phenotype)]
erg1[, phenotype := gsub("_"," ",phenotype)]

erg1[,table(phenotype,candidateGene)]

#' # Generate PlotData ####
#' ***
plotData1 = copy(erg1)
setorder(plotData1,"pval")
plotData1 = plotData1[!duplicated(markername),]
plotData1[,table(phenotype)]
plotData1[!is.na(candidateGene),table(candidateGene,phenotype)]
plotData1[,logP := -log10(pval)]

myNames<-c("markername","chr","bp_hg19","pval","logP","phenotype","candidateGene")
colsOut<-setdiff(colnames(plotData1),myNames)
plotData1[,get("colsOut"):=NULL]

setnames(plotData1,"markername","SNP")
setnames(plotData1,"chr","CHR")
setnames(plotData1,"bp_hg19","BP")
setnames(plotData1,"pval","P")
setorder(plotData1,CHR,BP)

head(plotData1)

save(plotData1, file = "../temp/SupFig_ManhattanPlot_PlotData.RData")

#' # Plot ####
#' ***
ymaxpar11 = ifelse(plotData1[,max(logP,na.rm=T)] <7, 8,plotData1[,max(logP,na.rm=T)]+1)

myPlot1 = manhattanPlot(x=plotData1,
                 ymax = ymaxpar11,
                 ymin = 0,
                 title = "",
                 xlabel = "",
                 ylabel=expression(-log[10](p)),
                 hline1=-log10(5e-8),
                 sugline1=-log10(1e-6),
                 highlight=T, diffsize = T,num_breaks_y=10,
                 returnObject = T, plotGenes=T,
                 overall_max = 20)

myPlot1

#' Save plots as PDF
pdf_from_png(code2parseOrPlot = myPlot1, 
             pdf_filename = "../figures/SupplementalFigure2_ManhattanPlot.pdf",
             weite = 12,
             laenge = 8,
             einheiten = "in",
             resolution = 150)

#' Save plots as tiff
tiff(filename = "../figures/SupplementalFigure2_ManhattanPlot.tiff", 
     width = 4800, height = 2440, res=300, compression = 'lzw')
myPlot1
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


