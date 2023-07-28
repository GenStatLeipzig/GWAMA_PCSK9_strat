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
#' **Plot**: PCSK9 (cond & uncond results)
#' 
#' * only the 12 rows for the 12 indep PCSK9 signals
#' * columns eQTLs and lipids
#'    * 9 tissues with H4>0.75
#'    * 3 lipid traits & CAD
#' 
#' => 12 x 13 matrix with PP H3 or H4 
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

load("../results/06_2_coloc_eQTLs.RData")
Coloc_eQTLs = copy(ColocTable)
load("../results/06_3_coloc_eQTLs_cond.RData")
Coloc_eQTLs_cond = copy(ColocTable)
load("../results/06_6_coloc_otherGWAS.RData")
Coloc_other = copy(ColocTable)
load("../results/06_7_coloc_otherGWAScond.RData")
Coloc_other_cond = copy(ColocTable)

#' # Combine data for PCSK9 plot ####
#' ***
#' ## Filter data 
Coloc_eQTLs_cond[,locus := "01p32.3"]
setnames(Coloc_eQTLs_cond,"locus","cytoband")
Coloc_eQTLs_1 = copy(Coloc_eQTLs)
Coloc_eQTLs_1 = Coloc_eQTLs_1[gene == "PCSK9"]
Coloc_eQTLs_1 = Coloc_eQTLs_1[trait1 %in% IndepSignals$pheno[c(11,12)]]
Coloc_eQTLs_1[,trait2 := gsub("GE in ","",trait2)]
Coloc_eQTLs_PCSK9 = rbind(Coloc_eQTLs_1,Coloc_eQTLs_cond,fill=T)
relTissues = Coloc_eQTLs_PCSK9[PP.H4.abf>0.75,unique(trait2)]
Coloc_eQTLs_PCSK9 = Coloc_eQTLs_PCSK9[trait2 %in% relTissues,]

Coloc_other_1 = copy(Coloc_other)
Coloc_other_1 = Coloc_other_1[gene == "PCSK9"]
Coloc_other_1 = Coloc_other_1[trait1 %in% IndepSignals$pheno[c(11,12)]]
Coloc_other_PCSK9 = rbind(Coloc_other_1,Coloc_other_cond,fill=T)
Coloc_other_PCSK9[,trait2 :=gsub("_.*","",trait2)]
relotherGWAS = Coloc_other_PCSK9[PP.H4.abf>0.75,unique(trait2)]
Coloc_other_PCSK9 = Coloc_other_PCSK9[trait2 %in% relotherGWAS,]

#' ## Get matrix format 
eQTL_H4<-dcast(Coloc_eQTLs_PCSK9,
                 formula = trait1 ~ trait2,
                 value.var = c("PP.H4.abf"),
                 sep = "_")
names(eQTL_H4)
eQTL_H3<-dcast(Coloc_eQTLs_PCSK9,
               formula = trait1 ~ trait2,
               value.var = c("PP.H3.abf"),
               sep = "_")
names(eQTL_H3)

lipi_H4<-dcast(Coloc_other_PCSK9,
               formula = trait1 ~ trait2,
               value.var = c("PP.H4.abf"),
               sep = "_")
names(lipi_H4)
lipi_H3<-dcast(Coloc_other_PCSK9,
               formula = trait1 ~ trait2,
               value.var = c("PP.H3.abf"),
               sep = "_")
names(lipi_H3)

table(eQTL_H4$trait1 == lipi_H4$trait1)
table(eQTL_H3$trait1 == lipi_H3$trait1)
H4 = cbind(eQTL_H4,lipi_H4[,-1])
H3 = cbind(eQTL_H3,lipi_H3[,-1])
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
eQTL_H3$trait1
#rownames(M) = eQTL_H3$trait1
rownames(M)<-c("rs11591147: F - free","rs693668: F - free","rs11591147: F - comb",
               "rs693668: F - comb","rs11591147: A - free","rs693668 : A - free",
               "rs11591147: M - free","rs693668 : M - free","rs11591147: M - comb",
               "rs693668: M - comb","rs28385704: M - treated", "rs693668: A - treated")
colnames(M)<-names(H3)[-1]
names(H3)[-1]
colnames(M)[1] = "Brain (CH)"
colnames(M)[2] = "Brain (C)"
colnames(M)[5] = "Skin (not exposed)"
colnames(M)[6] = "Skin (exposed)"
colnames(M)[9] = "Whole Blood"

col_fun = colorRamp2(c(-1,0,1), c("coral","white","steelblue"))
plot1 = Heatmap(M, name = "PP",col=col_fun, 
              column_km = 2,row_km = 2,
              row_title = c("", ""), column_title = c("","")
              )
plot1

tiff(filename = "../figures/MainFigure3_ColocPlot_PCSK9.tiff", 
     width = 2540, height = 1920, res=300, compression = 'lzw')
plot1
dev.off()

        
#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
