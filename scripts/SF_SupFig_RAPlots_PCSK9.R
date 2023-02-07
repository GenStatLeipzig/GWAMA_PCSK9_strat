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
#' RA plots for PCSK9 locus for all 7 traits in their conditional settings 
#' 
#' * page 1: females, females free, free
#' * page 2: males, males free, males treated, treated
#' 
#' # Initialize ####
#' ***
time0 = Sys.time()

source("../SourceFile_angmar.R")
source(paste0(basicpath,"/07_programme/rtools/1404_regional_assoc_plot_WO_broad/RegAssocPlot_hg19_EnsemblGenes_noBROAD_9.10_variablegengroesse.R"))

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/05_GCTA_COJO.RData")
IndepSignals = IndepSignals[candidateGene=="PCSK9",]
IndepSignals

dumTab = foreach(i = 1:dim(IndepSignals)[1])%do%{
  #i=1
  myRow = copy(IndepSignals)
  myRow = myRow[i,]
  message("Working on ",myRow$pheno," with signal ",myRow$SNP)
  
  erg1 = fread(myRow$input_CS)
  
  if(myRow$multipleSignals==T){
    rsID = gsub(":.*","",myRow$SNP)
    erg1[,phenotype := paste(myRow$pheno,rsID,sep="_")]
    setnames(erg1,"SNP","markername")
    setnames(erg1,"Chr","chr")
    setnames(erg1,"bp","bp_hg19")
    setnames(erg1,"refA","EA")
    setnames(erg1,"freq","EAF")
    setnames(erg1,"bC","beta")
    setnames(erg1,"bC_se","SE")
    setnames(erg1,"pC","pval")
    erg1 = erg1[!is.na(beta),]
  }else{
    erg1 = erg1[invalidAssocs==F,]
    erg1 = erg1[chr==1 & bp_hg19>= 55505647 - 500000 & bp_hg19<= 55505647 + 500000,] 
  }
  
  erg1
}

erg = rbindlist(dumTab,use.names=T,fill=T)

erg = erg[!is.na(pval)]
names(erg)[c(1:5,11:14)]
erg = erg[,c(1:5,11:14)]
setorder(erg,pval)
todofile = erg[, .SD[1], by=phenotype]
todofile[,gene := "PCSK9"]
todofile[,cyto := "1p32.3"]
todofile[,rsID := gsub(":.*","",markername)]
todofile[,region := 250]
todofile[,resultdir := paste0("../temp/input_RA_PCSK9/",phenotype)]  
todofile[,resultdir := gsub("PCSK9_","",resultdir)]  
todofile[,maf := EAF]  
todofile[EAF>0.5,maf := 1-EAF]  
setorder(todofile,phenotype)
todofile

#' # RA Plotting ####
#' ***
#' * PDF 1: females free, females, free (3x2 in typical PDF format)
#' * PDF 2: males free, males, males treated, treated (3x2 in typical PDF format)
#' 
uniqueSNPs = unique(todofile$markername)
todofile[,otherSignals := uniqueSNPs[c(4,4,4,4,4,4,4,4,4,4,3,4)]]

pdf("../figures/SupplementalFigure_RAPlots_PCSK9.pdf",8,12)
par(mfrow=c(3,2))

dumTab<-foreach(i = 1:dim(todofile)[1])%do%{
  #i=11
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
  input$myLocus = dummy
  customASplot_woBROAD(locus = input$myLocus,
                       lead_snp = input$leadsnp, 
                       map = input$myMap,
                       genes = input$myGenes,
                       center_lead_snp = F,
                       dolegend = T,
                       legendpos = "topright",
                       shownregion_kb = myRow$region,
                       #shownregion_kb = 250,
                       maintitle = mainText,
                       subtitle = subText,
                       weakR2 = 0.1,
                       cex_genname=1, title_size=1.5, subtitle_size = 0.8)
  
  points(myDat$position[myDat$snp %in% myRow$otherSignals], 
         -log10(myDat$pvalue)[myDat$snp %in% myRow$otherSignals], 
         pch=21,cex=2.5,lwd=3, col="red")
  
}
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
