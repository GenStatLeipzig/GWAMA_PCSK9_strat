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
#' RA plots for all other loci
#' 
#' * page 1: lipid related hits
#' * page 2: other hits
#' 
#' # Initialize ####
#' ***
time0 = Sys.time()

source("../SourceFile_angmar.R")
source("../helperFunctions/miamiPlot.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_230120.txt.gz","",pheno)]

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
result.5
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

erg1[phenotype == "PCSK9_females", phenotype := "PCSK9_females_all"]
erg1[phenotype == "PCSK9_males", phenotype := "PCSK9_males_all"]
erg1[phenotype == "PCSK9_free", phenotype := "PCSK9_all_free"]
erg1[phenotype == "PCSK9_treated", phenotype := "PCSK9_all_treated"]

erg1[,table(phenotype,candidateGene)]

#' # Generate PlotData ####
#' ***
#' Data set 1: females (statin combined, statin free, statin treated) vs males (statin combined, statin free, statin treated)
#' 
#' Data set 2: treated (sex combined, females, males) vs free (sex combined, females, males)
#'  
data11 = copy(erg1)
data11 = data11[grepl("female",phenotype)]
setorder(data11,"pval")
data11 = data11[!duplicated(markername),]
data11[,table(phenotype)]
data11[!is.na(candidateGene),table(candidateGene,phenotype)]

data12 = copy(erg1)
data12 = data12[grepl("_male",phenotype)]
setorder(data12,"pval")
data12 = data12[!duplicated(markername),]
data12[,table(phenotype)]
data12[!is.na(candidateGene),table(candidateGene,phenotype)]

data21 = copy(erg1)
data21 = data21[grepl("treated",phenotype)]
setorder(data21,"pval")
data21 = data21[!duplicated(markername),]
data21[,table(phenotype)]
data21[!is.na(candidateGene),table(candidateGene,phenotype)]

data22 = copy(erg1)
data22 = data22[grepl("free",phenotype)]
setorder(data22,"pval")
data22 = data22[!duplicated(markername),]
data22[,table(phenotype)]
data22[!is.na(candidateGene),table(candidateGene,phenotype)]

#' Pairwise merge 
plotData1 = rbind(data11,data12)
plotData2 = rbind(data21,data22)
plotData1[,logP := -log10(pval)]
plotData2[,logP := -log10(pval)]

myNames<-c("markername","chr","bp_hg19","pval","logP","phenotype","candidateGene")
colsOut<-setdiff(colnames(plotData1),myNames)
plotData1[,get("colsOut"):=NULL]
plotData2[,get("colsOut"):=NULL]

plotData1[grepl("_female",phenotype),flag:="top"]
plotData1[grepl("_male",phenotype),flag:="bottom"]
plotData2[grepl("_treated",phenotype),flag:="top"]
plotData2[grepl("_free",phenotype),flag:="bottom"]

setnames(plotData1,"markername","SNP")
setnames(plotData1,"chr","CHR")
setnames(plotData1,"bp_hg19","BP")
setnames(plotData1,"pval","P")
setorder(plotData1,CHR,BP)
setnames(plotData2,"markername","SNP")
setnames(plotData2,"chr","CHR")
setnames(plotData2,"bp_hg19","BP")
setnames(plotData2,"pval","P")
setorder(plotData2,CHR,BP)

head(plotData1)
head(plotData2)

# check for duplicate candidate genes
plotData1[,dumID := paste(candidateGene,flag, sep="::")]
plotData1[duplicated(dumID) & !is.na(candidateGene)]
plotData1[dumID=="PCSK9::bottom"]
plotData1[dumID=="PCSK9::bottom" & BP == 55506103, candidateGene:=NA]
plotData1[dumID=="SLCO1B3::top"]
plotData1[dumID=="SLCO1B3::top" & BP == 21117156, candidateGene:=NA]

plotData2[,dumID := paste(candidateGene,flag, sep="::")]
plotData2[duplicated(dumID) & !is.na(candidateGene)]
plotData2[dumID=="PCSK9::top"]
plotData2[dumID=="PCSK9::top" & BP == 55521109 , candidateGene:=NA]

save(plotData1, plotData2,file = "../temp/SupFig_MiamiPlot_PlotData.RData")

#' # Plot ####
#' ***
ymaxpar11 = ifelse(plotData1[flag=="top",max(logP,na.rm=T)] <7, 8,plotData1[flag=="top",max(logP,na.rm=T)]+1)
ymaxpar12 = ifelse(plotData1[flag=="bottom",max(logP,na.rm=T)] <7, 8,plotData1[flag=="bottom",max(logP,na.rm=T)]+1)
ymaxpar21 = ifelse(plotData2[flag=="top",max(logP,na.rm=T)] <7, 8,plotData2[flag=="top",max(logP,na.rm=T)]+1)
ymaxpar22 = ifelse(plotData2[flag=="bottom",max(logP,na.rm=T)] <7, 8,plotData2[flag=="bottom",max(logP,na.rm=T)]+1)

myPlot1 = miamiPlot(x=plotData1,
                 ymax = ymaxpar11,
                 ymin = -ymaxpar12,
                 title = "",
                 xlabel = "",
                 ylabel=expression(paste("males: ",log[10](p),"; females: ",-log[10](p))),
                 hline1=-log10(5e-8),hline2=log10(5e-8),
                 sugline1=-log10(1e-6),sugline2=log10(1e-6),
                 highlight=T, diffsize = T,num_breaks_y=10,
                 returnObject = T,
                 overall_max = 15, overall_min = -15)

myPlot1

myPlot2 = miamiPlot(x=plotData2,
                    ymax = ymaxpar21,
                    ymin = -ymaxpar22,
                    title = "",
                    xlabel = "",
                    ylabel=expression(paste("free: ",log[10](p),"; treated: ",-log[10](p))),
                    hline1=-log10(5e-8),hline2=log10(5e-8),
                    sugline1=-log10(1e-6),sugline2=log10(1e-6),
                    highlight=T, diffsize = T,num_breaks_y=10,
                    returnObject = T,
                    overall_max = 15, overall_min = -15)

myPlot2

#' Save plots as PDF
pdf_from_png(code2parseOrPlot = myPlot1, 
             pdf_filename = "../figures/SupplementalFigure_MiamiPlot_females_males.pdf",
             weite = 12,
             laenge = 8,
             einheiten = "in",
             resolution = 150)

pdf_from_png(code2parseOrPlot = myPlot2, 
             pdf_filename = "../figures/SupplementalFigure_MiamiPlot_treated_free.pdf",
             weite = 12,
             laenge = 8,
             einheiten = "in",
             resolution = 150)

#' Save plots as tiff
tiff(filename = "../figures/SupplementalFigure_MiamiPlot_females_males.tiff", 
     width = 4800, height = 2440, res=300, compression = 'lzw')
myPlot1
dev.off()

tiff(filename = "../figures/SupplementalFigure_MiamiPlot_treated_free.tiff", 
     width = 4800, height = 2440, res=300, compression = 'lzw')
myPlot2
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


