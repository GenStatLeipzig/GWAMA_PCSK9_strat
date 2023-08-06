#' ---
#' title: "Main Figure 1: Summary of valid locus"
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
#' **Figure 2: Heatmap of the log10-transformed p-values at all valid loci**
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Load overview data ####
#' ***
load("../temp/04_IATest_input.RData")
result.2[, dumID1 := paste(candidateGene,rsID,sep = " - ")]
result.2[, dumID2 := paste(phenotype,rsID,sep="_")]

result.5 = copy(result.2)
setorder(result.5,pval)
result.5 = result.5[!duplicated(markername),]

negEffectSNPs = result.5[beta<0, markername]

result.2[, beta2 := beta]
result.2[markername %in% negEffectSNPs, beta2 := beta*(-1)]
result.2[, Zscore := beta/SE]
result.2[, Zscore2 := beta2/SE]
result.2[, logP := -log10(pval)]
result.2[, logP2 := 0]
result.2[logP> -log10(0.05), logP2 := 1.3]
result.2[logP> -log10(1e-6), logP2 := 6]
result.2[logP> -log10(5e-8), logP2 := 7.3]
result.2[dumID2 %in% result.5$dumID2 & logP> -log10(5e-8), logP2 := 8]
result.2[dumID2 %in% result.5$dumID2 & logP< -log10(5e-8), logP2 := 7]

#' # Create Heatmap ####
#' ***
#' ## Normal beta ####
plotData1<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("beta"),
                 sep = "_")

SNPID = as.character(plotData1$dumID1)
phenotypes = colnames(plotData1)[-1]
phenotypes = gsub("PCSK9_","",phenotypes)
phenotypes = gsub("females","F",phenotypes)
phenotypes = gsub("males","M",phenotypes)
phenotypes = gsub("_"," - ",phenotypes)

plot.matrix1 = as.matrix(plotData1[, -1])
rownames(plot.matrix1) = SNPID
colnames(plot.matrix1) = phenotypes
#plot.matrix1 = t(plot.matrix1)

col_fun1 = colorRamp2(c(-0.4,0,0.06), c("coral","white","steelblue"))
Heatmap(plot.matrix1, name = "beta",col=col_fun1)

#' ## Positive beta ####
plotData2<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("beta2"),
                 sep = "_")

plot.matrix2 = as.matrix(plotData2[, -1])
rownames(plot.matrix2) = SNPID
colnames(plot.matrix2) = phenotypes

col_fun2 = colorRamp2(c(0,0.4), c("white","steelblue"))
Heatmap(plot.matrix2, name = "abs(beta)",col=col_fun2)

#' ## Z score ####
plotData3<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("Zscore"),
                 sep = "_")

plot.matrix3 = as.matrix(plotData3[, -1])
rownames(plot.matrix3) = SNPID
colnames(plot.matrix3) = phenotypes
#plot.matrix3 = t(plot.matrix3)

col_fun3 = colorRamp2(c(-30,0,15), c("coral","white","steelblue"))
Heatmap(plot.matrix3, name = "Zscore",col=col_fun3)

#' ## Positive Z score ####
plotData4<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("Zscore2"),
                 sep = "_")

plot.matrix4 = as.matrix(plotData4[, -1])
rownames(plot.matrix4) = SNPID
colnames(plot.matrix4) = phenotypes
#plot.matrix4 = t(plot.matrix4)

col_fun4 = colorRamp2(c(0,30), c("white","steelblue"))
Heatmap(plot.matrix4, name = "abs(Zscore)",col=col_fun4)

#' ## Pvalue ####
plotData5<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("logP"),
                 sep = "_")

plot.matrix5 = as.matrix(plotData5[, -1])
rownames(plot.matrix5) = SNPID
colnames(plot.matrix5) = phenotypes
#plot.matrix5 = t(plot.matrix5)

col_fun5 = colorRamp2(c(0,1.3,6,7.3), c("white","yellow","orange","red"))
Heatmap(t(plot.matrix5), name = "-logP",col=col_fun5)

#' ## Pvalue 2 ####
plotData6<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("logP2"),
                 sep = "_")

plot.matrix6 = as.matrix(plotData6[, -1])
rownames(plot.matrix6) = SNPID
colnames(plot.matrix6) = phenotypes
#plot.matrix6 = t(plot.matrix6)

col_fun6 = colorRamp2(c(0,1.3,6,7,7.3,8), c("#FFFFB2", "#FED976", "#FEB24C" ,"#FD8D3C", "#F03B20", "#BD0026"))
Heatmap(plot.matrix6, name = "-logP",col=col_fun6)

Heatmap(plot.matrix6, name = "-logP",col=col_fun6,show_column_names = FALSE, 
        bottom_annotation = HeatmapAnnotation(text = anno_text(phenotypes, 
                                                               rot = 45, 
                                                               location = unit(1, "npc"), 
                                                               just = "right"),
                                              annotation_height = max_text_width(phenotypes)))

plot1 = Heatmap(plot.matrix6, name = "-logP",col=col_fun6, 
                clustering_distance_rows="euclidean",
                clustering_method_columns = "complete",
                show_column_dend = FALSE, show_row_dend = FALSE,
                show_heatmap_legend = F,show_column_names = FALSE, 
                bottom_annotation = HeatmapAnnotation(text = anno_text(phenotypes, 
                                                                       rot = 45, 
                                                                       location = unit(1, "npc"), 
                                                                       just = "right"),
                                                      annotation_height = max_text_width(phenotypes)))
plot1

#' # Save plot ####
#' ***
tiff(filename = "../figures/MainFigure1_Heatmap_logp.tiff",
     width = 1200, height = 1400, res=250, compression = 'lzw')
plot1
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
