#' ---
#' title: "Supplemental Figures: Regional Association Plots (2)"
#' subtitle: "PCSK9 GWAMA stratified"
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
source(path_RAPlot_function)

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/03_GCTA_COJO.RData")
IndepSignals = IndepSignals[candidateGene!="PCSK9",]
IndepSignals[,dumID := paste(candidateGene,pheno,sep="_")]

load("../temp/05_GWAS_allLoci.RData")
erg = copy(data_GWAS)
erg = erg[candidateGene != "PCSK9"]
erg = erg[!is.na(pval)]
names(erg)[c(1:4,6,10:12,16,17)]
erg = erg[,c(1:4,6,10:12,16,17)]
erg[,dumID := paste(candidateGene,phenotype,sep="_")]
erg = erg[dumID %in% IndepSignals$dumID]
erg[,table(phenotype,candidateGene)]
setorder(erg,pval)

erg = erg[markername != "12:21117156:GAA:GA",]
todofile = erg[, .SD[1], by=candidateGene]
todofile[,gene := candidateGene]
todofile[,lipids := c(T,T,T,F,F, F,T,T,F,F)]
todofile[,cyto := c("2p24.1","19p13.11","11q12.2","12q24.22","6q11.1","10q11.21",
                    "16q22.2","10q21.3","12p12.2","7q36.1")]
todofile[,rsID := gsub(":.*","",markername)]
todofile[,region := 250]
todofile[,resultdir := paste0("../temp/input_RA_otherLoci/",gene)]  
todofile[,maf := EAF]  
todofile[EAF>0.5,maf := 1-EAF]  
setorder(todofile,-lipids)
todofile

#' # RA Plotting ####
#' ***
#' * PDF 1: lipid loci
#' * PDF 2: other loci
#' 
filename1 = "../figures/SupplementalFigure_RAPlots_lipidLoci.pdf"
filename2 = "../figures/SupplementalFigure_RAPlots_otherLoci.pdf"

todofile1 = copy(todofile)
todofile1 = todofile1[lipids == T,]

todofile2 = copy(todofile)
todofile2 = todofile2[lipids == F,]

#' ## Lipid loci ####
pdf(filename1,8,12)
par(mfrow=c(3,2))

dumTab<-foreach(i = 1:dim(todofile1)[1])%do%{
  #i=5
  message("Working on candidate gene ",todofile1$candidateGene[i])
  myRow<-todofile1[i,]
  
  # prep data
  myDat = copy(erg)
  myDat = myDat[candidateGene == myRow$gene]
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
  dummy2 = input$myGenes
  dummy2 = dummy2[dummy2$biotype=="protein_coding",]
  input$myGenes = dummy2
  
  if(i %in% c(1,4,5)){
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = F,
                         legendpos = "topleft",
                         shownregion_kb = myRow$region,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.6, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 2){
    dummy = input$myLocus
    dummy = dummy[dummy$POS<=19843906,]
    dummy = dummy[dummy$POS>=19179794,]
    input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = F,
                         legendpos = "topleft",
                         #shownregion_kb = myRow$region,
                         #shownregion_kb = 400,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.4, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 3){
    dummy = input$myLocus
    dummy = dummy[dummy$POS>=61426711,]
    dummy = dummy[dummy$POS<=61686500,]
    input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = F,
                         legendpos = "topleft",
                         #shownregion_kb = myRow$region,
                         #shownregion_kb = 400,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.6, title_size=1.5, subtitle_size = 0.8)
  }
}
dev.off()


#' ## Other loci ####
pdf(filename2,8,12)
par(mfrow=c(3,2))

dumTab<-foreach(i = 1:dim(todofile2)[1])%do%{
  #i=5
  message("Working on candidate gene ",todofile2$candidateGene[i])
  myRow<-todofile2[i,]
  
  # prep data
  myDat = copy(erg)
  myDat = myDat[candidateGene == myRow$gene]
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
  dummy2 = input$myGenes
  dummy2 = dummy2[dummy2$biotype=="protein_coding",]
  input$myGenes = dummy2
  
  customASplot_woBROAD(locus = input$myLocus,
                       lead_snp = input$leadsnp, 
                       map = input$myMap,
                       genes = input$myGenes,
                       center_lead_snp = F,
                       dolegend = F,
                       legendpos = "topleft",
                       shownregion_kb = myRow$region,
                       maintitle = mainText,
                       subtitle = subText,
                       weakR2 = 0.1,
                       cex_genname=0.6, title_size=1.5, subtitle_size = 0.8)
  
  
}
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
