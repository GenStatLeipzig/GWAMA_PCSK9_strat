#' ---
#' title: "Supplemental Figures: Regional Association Plots (1)"
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
#' RA plots for PCSK9 locus for all 8 traits in their unconditional settings 
#'  
#' # Initialize ####
#' ***
time0 = Sys.time()

source("../SourceFile_angmar.R")
source(path_RAPlot_function)

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/03_GCTA_COJO.RData")
IndepSignals = IndepSignals[candidateGene=="PCSK9",]
IndepSignals

load("../temp/05_GWAS_allLoci.RData")
erg = copy(data_GWAS)
erg = erg[candidateGene == "PCSK9"]
erg = erg[!is.na(pval)]
names(erg)[c(1:4,6,10:12,16,17)]
erg = erg[,c(1:4,6,10:12,16,17)]
erg[,table(phenotype,candidateGene)]
setorder(erg,pval)

todofile = erg[, .SD[1], by=phenotype]
todofile[,gene := candidateGene]
todofile[,cyto := c("1p32.3")]
todofile[,rsID := gsub(":.*","",markername)]
todofile[,region := 250]
todofile[,resultdir := paste0("../temp/input_RA_PCSK9/",phenotype)]  
todofile[,maf := EAF]  
todofile[EAF>0.5,maf := 1-EAF]  
todofile = todofile[c(1,6,2,4,3,5,7,8),]

#' # RA Plotting ####
#' ***

pdf("../figures/SupplementalFigure_RAPlots_PCSK9.pdf",8,12)
par(mfrow=c(3,2))

dumTab<-foreach(i = 1:dim(todofile)[1])%do%{
  #i=1
  message("Working on phenotype ",todofile$phenotype[i])
  myRow<-todofile[i,]
  
  # prep data
  myDat = copy(erg)
  myDat = myDat[phenotype == myRow$phenotype]
  setnames(myDat,"markername","snp")
  setnames(myDat,"bp_hg19","position")
  setnames(myDat,"pval","pvalue")
  myDat<-myDat[!is.na(myDat$pvalue),]
  
  # create input
  input = createInputfilesRegAssocPlot1kg19Ensembl(myDat, 
                                                   doplink=T, 
                                                   r_on_server = T,
                                                   gene_validation_level = c("KNOWN"),
                                                   path_ldreference = path_1000Genomes, 
                                                   path_ldreference_snps = path_1000Genomes_snps_fn,
                                                   leadsnp = myRow$markername,
                                                   resultdir = myRow$resultdir)
  
  # plotting
  subText=paste0(myRow$rsID," (MAF: ",round(myRow$maf,3),")")
  mainText=gsub("PCSK9_","",myRow$phenotype )
  mainText=gsub(myRow$rsID,"",mainText )
  mainText=gsub("_"," ",mainText )
  #par(mfrow=c(1,1))
  
  dummy = input$myLocus
  dummy = dummy[dummy$POS>=55443166,]
  dummy = dummy[dummy$POS<=55599526,]
  input$myLocus = dummy
  
  dummy2 = input$myGenes
  dummy2 = dummy2[dummy2$biotype=="protein_coding",]
  input$myGenes = dummy2
  
  customASplot_woBROAD(locus = input$myLocus,
                       lead_snp = input$leadsnp, 
                       map = input$myMap,
                       genes = input$myGenes,
                       center_lead_snp = F,
                       dolegend = F,
                       legendpos = "topright",
                       shownregion_kb = myRow$region,
                       #shownregion_kb = 250,
                       maintitle = mainText,
                       subtitle = subText,
                       weakR2 = 0.1,
                       cex_genname=1, title_size=1.5, subtitle_size = 0.8)
  
  Indep2 = copy(IndepSignals)
  Indep2 = Indep2[pheno == myRow$phenotype,]
  Indep2 = Indep2[SNP != myRow$markername,]
  points(Indep2$bp,-log10(Indep2$p), 
         pch=21,cex=2.5,lwd=3, col="red")
  
}
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
