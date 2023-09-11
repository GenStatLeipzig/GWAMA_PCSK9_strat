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
result.2[, logP3 := 0]
result.2[logP> -log10(0.05), logP3 := 1.3]
result.2[logP> -log10(1e-6), logP3 := 6]
result.2[logP> -log10(5e-8), logP3 := 7.3]
# result.2[dumID2 %in% result.5$dumID2 & logP> -log10(5e-8), comment := "bestPheno"]
# result.2[dumID2 %in% result.5$dumID2 & logP< -log10(5e-8), comment := "bestPheno"]

#' # Create Heatmap ####
#' ***
plotData6<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("logP2"),
                 sep = "_")

SNPID = as.character(plotData6$dumID1)
phenotypes = colnames(plotData6)[-1]
phenotypes = gsub("PCSK9_","",phenotypes)
phenotypes = gsub("females","F",phenotypes)
phenotypes = gsub("males","M",phenotypes)
phenotypes = gsub("_"," - ",phenotypes)

plot.matrix6 = as.matrix(plotData6[, -1])
rownames(plot.matrix6) = SNPID
colnames(plot.matrix6) = phenotypes
plotData7<-dcast(result.2,
                 formula = dumID1 ~ phenotype,
                 value.var = c("logP3"),
                 sep = "_")
plot.matrix7 = as.matrix(plotData7[, -1])
rownames(plot.matrix7) = SNPID
colnames(plot.matrix7) = phenotypes

#resort phenotypes
plot.matrix7 = plot.matrix7[,c(4,5,6,1,2,8,7,3)]
plot.matrix6 = plot.matrix6[,c(4,5,6,1,2,8,7,3)]

#resort SNPs (as in Table 1)
plot.matrix7 = plot.matrix7[c(9,11,8,10,2,5,3,4,14,6,12,1,13,7),]
plot.matrix6 = plot.matrix6[c(9,11,8,10,2,5,3,4,14,6,12,1,13,7),]

#col_fun6 = colorRamp2(c(0,1.3,6,7,7.3,8), c("#FFFFB2", "#FED976", "#FEB24C" ,"#FD8D3C", "#F03B20", "#BD0026"))
col_fun6 = colorRamp2(c(0,1.3,6,7.3), c("#FFFFB2", "#FED976", "#FD8D3C",  "#BD0026"))
Heatmap(plot.matrix7, name = "-logP",col=col_fun6)

small_mat = plot.matrix6
small_mat[small_mat==8] = "*"
small_mat[small_mat==7] = "*"
small_mat[small_mat!="*"] = ""

Heatmap(plot.matrix7, name = "-logP",col=col_fun6,show_column_names = FALSE, 
        bottom_annotation = HeatmapAnnotation(text = anno_text(phenotypes, 
                                                               rot = 45, 
                                                               location = unit(1, "npc"), 
                                                               just = "right"),
                                              annotation_height = max_text_width(phenotypes)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(small_mat[i, j], x, y, gp = gpar(fontsize = 15))})


plot1 = Heatmap(plot.matrix7, name = "-logP",col=col_fun6, 
                cluster_rows = F, cluster_columns = F,
                # clustering_distance_row = "euclidean",
                # clustering_method_columns = "complete",
                show_column_dend = FALSE, show_row_dend = FALSE,
                show_heatmap_legend = F,show_column_names = FALSE, 
                bottom_annotation = HeatmapAnnotation(text = anno_text(colnames(plot.matrix7), 
                                                                       rot = 45, 
                                                                       location = unit(1, "npc"), 
                                                                       just = "right"),
                                                      annotation_height = max_text_width(colnames(plot.matrix7))),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(small_mat[i, j], x, y, gp = gpar(fontsize = 15))})
plot1

#' # Save plot ####
#' ***
tiff(filename = "../figures/MainFigure1_Heatmap_logp_230904.tiff",
     width = 1200, height = 1400, res=250, compression = 'lzw')
plot1
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
