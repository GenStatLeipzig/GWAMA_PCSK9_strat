#' ---
#' title: "Main Figure 2: Summary of PCSK9 locus"
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
#' **Figure 2: Heatmap of the log10-transformed p-values at the PCSK9 locus**
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
load("../temp/06_GWAS_allLoci.RData")
load("../results/05_GCTA_COJO.RData")
IndepSignals = IndepSignals[candidateGene == "PCSK9",]
table(IndepSignals$pheno)
data_GWAS = data_GWAS[markername %in% IndepSignals$SNP,]
table(data_GWAS$phenotype,data_GWAS$markername)

#' # Load GWAMA data ####
#' ***
#' 
erg_fem_tre = fread("../data/SumStat_PCSK9_females_treated_230120.txt.gz")
erg_fem_tre = erg_fem_tre[markername %in% IndepSignals$SNP,]
erg_fem_tre = erg_fem_tre[invalidAssocs==T,]
erg_fem_tre[,candidateGene := "PCSK9"]
erg_fem_tre[,setting_sex := "females"]
erg_fem_tre[,setting_statin := "treated"]
erg_fem_tre[,MAF := EAF]
erg_fem_tre[EAF>0.5,MAF := 1-EAF]
data_GWAS = rbind(data_GWAS,erg_fem_tre)
table(data_GWAS$phenotype,data_GWAS$markername)

erg_mal_tre = fread("../data/SumStat_PCSK9_males_treated_230120.txt.gz")
erg_mal_tre = erg_mal_tre[markername %in% IndepSignals$SNP,]
erg_mal_tre = erg_mal_tre[invalidAssocs==T,]
erg_mal_tre[,candidateGene := "PCSK9"]
erg_mal_tre[,setting_sex := "males"]
erg_mal_tre[,setting_statin := "treated"]
erg_mal_tre[,MAF := EAF]
erg_mal_tre[EAF>0.5,MAF := 1-EAF]
data_GWAS = rbind(data_GWAS,erg_mal_tre)
table(data_GWAS$phenotype,data_GWAS$markername)

erg_mal = fread("../data/SumStat_PCSK9_males_230120.txt.gz")
erg_mal = erg_mal[markername %in% IndepSignals$SNP,]
erg_mal = erg_mal[invalidAssocs==T,]
erg_mal[,candidateGene := "PCSK9"]
erg_mal[,setting_sex := "males"]
erg_mal[,setting_statin := "combined"]
erg_mal[,MAF := EAF]
erg_mal[EAF>0.5,MAF := 1-EAF]
data_GWAS = rbind(data_GWAS,erg_mal)
table(data_GWAS$phenotype,data_GWAS$markername)

erg_tre = fread("../data/SumStat_PCSK9_treated_230120.txt.gz")
erg_tre = erg_tre[markername %in% IndepSignals$SNP,]
erg_tre = erg_tre[invalidAssocs==T,]
erg_tre[,candidateGene := "PCSK9"]
erg_tre[,setting_sex := "combined"]
erg_tre[,setting_statin := "treated"]
erg_tre[,MAF := EAF]
erg_tre[EAF>0.5,MAF := 1-EAF]
data_GWAS = rbind(data_GWAS,erg_tre)
table(data_GWAS$phenotype,data_GWAS$markername)

#' # Create Heatmap ####
#' ***
#' of unconditional betas
plotData1<-dcast(data_GWAS,
                 formula = markername ~ phenotype,
                 value.var = c("beta"),
                 sep = "_")
SNPID = as.character(plotData1$markername)
SNPID = gsub(":.*","",SNPID)
phenotypes = names(plotData1)[-1]
phenotypes = gsub("PCSK9_","",phenotypes)
phenotypes = gsub("_"," ",phenotypes)
plot.matrix = as.matrix(plotData1[, -1])
rownames(plot.matrix) = SNPID
colnames(plot.matrix) = phenotypes
plot.matrix1 = t(plot.matrix)

col_fun1 = colorRamp2(c(-0.34,0,0.07), c("coral","white","steelblue"))
Heatmap(plot.matrix1, name = "beta",col=col_fun1)
plot1 = Heatmap(plot.matrix1, name = "beta",col=col_fun1, 
                clustering_distance_rows="euclidean", clustering_method_columns = "complete", 
                column_km = 3,row_km = 2,
                row_title = c("", ""), column_title = c("","","")
)
plot1

#' of unconditional log10 transformed pvalues
data_GWAS[,logP := -log10(pval)]
plotData2<-dcast(data_GWAS,
                 formula = markername ~ phenotype,
                 value.var = c("logP"),
                 sep = "_")
SNPID = as.character(plotData2$markername)
SNPID = gsub(":.*","",SNPID)
phenotypes = names(plotData2)[-1]
phenotypes = gsub("PCSK9_","",phenotypes)
phenotypes = gsub("_"," ",phenotypes)
plot.matrix = as.matrix(plotData2[, -1])
rownames(plot.matrix) = SNPID
colnames(plot.matrix) = phenotypes
plot.matrix2 = t(plot.matrix)

col_fun2 = colorRamp2(c(0,6,7.3, 15, 40), c("white","yellow", "orange", "red","darkred"))
plot2 = Heatmap(plot.matrix2, name = "log(p)",col=col_fun2, 
                column_km = 3,row_km = 2,
                row_title = c("", ""), column_title = c("","","")
)
plot2

plot.matr2.2 = copy(plot.matrix2)
plot.matr2.2[plot.matr2.2<1.3]=0 
plot.matr2.2[plot.matr2.2>1.3 & plot.matr2.2<6]=1 
plot.matr2.2[plot.matr2.2>6 & plot.matr2.2<7.3]=2 
plot.matr2.2[plot.matr2.2>7.3]=3 
plot.matr2.2[4,1]=4
plot.matr2.2[7,2]=4
plot.matr2.2[4,3]=4
plot.matr2.2[4,4]=4

col_fun2.2 = colorRamp2(c(0,1,2, 3,4), c("white","yellow", "orange", "red","darkred"))
Heatmap(plot.matr2.2, name = "log(p)",col=col_fun2.2) 
plot2.2 = Heatmap(plot.matr2.2, col=col_fun2.2, 
                column_km = 3,row_km = 2,
                row_title = c("", ""), column_title = c("","",""),
                show_heatmap_legend = F#,
                # show_row_dend = F,
                # show_column_dend = F
)
plot2.2

#' of unconditional Zscores
data_GWAS[,Zscore := beta/SE]
plotData3<-dcast(data_GWAS,
                 formula = markername ~ phenotype,
                 value.var = c("Zscore"),
                 sep = "_")
SNPID = as.character(plotData3$markername)
SNPID = gsub(":.*","",SNPID)
phenotypes = names(plotData3)[-1]
phenotypes = gsub("PCSK9_","",phenotypes)
phenotypes = gsub("_"," ",phenotypes)
plot.matrix = as.matrix(plotData3[, -1])
rownames(plot.matrix) = SNPID
colnames(plot.matrix) = phenotypes
plot.matrix3 = t(plot.matrix)

col_fun3 = colorRamp2(c(-15,0,15), c("coral","white","steelblue"))
Heatmap(plot.matrix3, name = "Zscore",col=col_fun3)
plot3 = Heatmap(plot.matrix3, name = "Zscore",col=col_fun3, 
                column_km = 3,row_km = 2,
                row_title = c("", ""), column_title = c("","","")
)
plot3

#' 
#' # Save plot ####
#' ***
tiff(filename = "../figures/MainFigure2_Heatmap_logp.tiff",
     width = 1650, height = 1350, res=250, compression = 'lzw')
plot2
dev.off()

tiff(filename = "../figures/MainFigure2_Heatmap_betas_230326.tiff",
     width = 1650, height = 1350, res=250, compression = 'lzw')
plot1
dev.off()

tiff(filename = "../figures/MainFigure2_Heatmap_logp_230326.tiff",
     width = 1650, height = 1350, res=250, compression = 'lzw')
plot2.2
dev.off()

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))
