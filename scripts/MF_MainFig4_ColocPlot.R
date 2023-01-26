#' ---
#' title: "Co-localization: Coloc Plot"
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
#' I want to plot the co-localization results for the PCSK9 locus (all PCSK9 traits).
#' 
#' **Plot**: lipid loci
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/05_GCTA_COJO.RData")

load("../results/06_6_coloc_otherGWAS.RData")
Coloc_other = copy(ColocTable)

#' # Combine data for PCSK9 plot ####
#' ***
#' ## Filter data 
Coloc_other_1 = copy(Coloc_other)
Coloc_other_1 = Coloc_other_1[grepl("INV",trait2),]
Coloc_other_1 = Coloc_other_1[gene != "PCSK9"]
Coloc_other_1[,trait2 :=gsub("_.*","",trait2)]
relotherGWAS = Coloc_other_1[PP.H4.abf>0.75 | PP.H3.abf>0.75,unique(gene)]
Coloc_other_1 = Coloc_other_1[gene %in% relotherGWAS,]

#' ## Get matrix format 
lipi_H4<-dcast(Coloc_other_1,
               formula = gene ~ trait2,
               value.var = c("PP.H4.abf"),
               sep = "_")
names(lipi_H4)
lipi_H3<-dcast(Coloc_other_1,
               formula = gene ~ trait2,
               value.var = c("PP.H3.abf"),
               sep = "_")
names(lipi_H3)

H4 = cbind(lipi_H4)
H3 = cbind(lipi_H3)
M4<-as.matrix(H4[,-1])
M3<-as.matrix(H3[,-1])

x1 = dim(M4)[1]
x2 = dim(M4)[2]

M<-matrix(0,x1,x2)
for (i in 1:x1){
  for (j in 1:x2){
    m4<-M4[i,j]
    m3<-M3[i,j]
    
    if(is.na(m3)==T){
      M[i,j] = NA
    }else if(m4>m3){
      M[i,j]<-m4
    }else{
      M[i,j]<- -m3
    }
  }
}

#' ## Get plot ####
lipi_H3$gene
rownames(M) = lipi_H3$gene
colnames(M)<-names(H3)[-1]

col_fun = colorRamp2(c(-1,0,1), c("coral","white","steelblue"))
plot1 = Heatmap(M, name = "PP",col=col_fun, 
              column_km = 2,row_km = 2,
              row_title = c("", ""), column_title = c("","")
              )
plot1

tiff(filename = "../figures/MainFigure4_ColocPlot_otherLoci.tiff", 
     width = 2540, height = 1920, res=300, compression = 'lzw')
plot1
dev.off()

        
#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
