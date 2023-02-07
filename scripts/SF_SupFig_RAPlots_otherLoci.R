#' ---
#' title: "Supplemental Figures: Regional Association Plots (2)"
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
source(paste0(basicpath,"/07_programme/rtools/1404_regional_assoc_plot_WO_broad/RegAssocPlot_hg19_EnsemblGenes_noBROAD_9.10_variablegengroesse.R"))

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
load("../results/05_GCTA_COJO.RData")
IndepSignals = IndepSignals[candidateGene!="PCSK9",]
IndepSignals[,dumID := paste(candidateGene,pheno,sep="_")]

load("../temp/06_GWAS_allLoci.RData")
erg = copy(data_GWAS)
erg = erg[candidateGene != "PCSK9"]
erg = erg[!is.na(pval)]
names(erg)[c(1:4,6,10:12,16,17)]
erg = erg[,c(1:4,6,10:12,16,17)]
erg[,dumID := paste(candidateGene,phenotype,sep="_")]
erg = erg[dumID %in% IndepSignals$dumID]
erg[,table(phenotype,candidateGene)]
erg[candidateGene == "gene desert",candidateGene := "geneDesert"]
setorder(erg,pval)

todofile = erg[, .SD[1], by=candidateGene]
todofile[,gene := candidateGene]
todofile[,lipids := c(T,T,F,F,T,F,T,F,F,T,T,F)]
todofile[,cyto := c("2p24.1","19p13.11","6q11.1","6q24.3","2p23.2","18q23",
                    "10q21.3","2p24.3","7q36.1","12p12.2","12q24.22","20p12.1")]

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
pdf("../figures/SupplementalFigure_RAPlots_otherLoci.pdf",8,12)
par(mfrow=c(3,2))

dumTab<-foreach(i = 1:dim(todofile)[1])%do%{
  #i=7
  message("Working on candidate gene ",todofile$candidateGene[i])
  myRow<-todofile[i,]
  
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
  
  if(i %nin% c(2,4,5,7,10)){
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = T,
                         legendpos = "topleft",
                         shownregion_kb = myRow$region,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=1, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 2){
    dummy2 = input$myGenes
    dummy2 = dummy2[dummy2$biotype=="protein_coding",]
    input$myGenes = dummy2
    dummy = input$myLocus
    dummy = dummy[dummy$POS<=19843906,]
    dummy = dummy[dummy$POS>=19179794,]
    input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = T,
                         legendpos = "topleft",
                         #shownregion_kb = myRow$region,
                         #shownregion_kb = 400,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.4, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 4){
    dummy2 = input$myGenes
    dummy2 = dummy2[dummy2$biotype=="protein_coding",]
    input$myGenes = dummy2
    dummy = input$myLocus
    dummy = dummy[dummy$POS<=65521463,]
    input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = T,
                         legendpos = "topleft",
                         #shownregion_kb = myRow$region,
                         #shownregion_kb = 400,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.6, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 5){
    # dummy2 = input$myGenes
    # dummy2 = dummy2[dummy2$biotype=="protein_coding",]
    # input$myGenes = dummy2
    dummy = input$myLocus
    dummy = dummy[dummy$POS>=20829710,]
    input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = T,
                         legendpos = "topleft",
                         #shownregion_kb = myRow$region,
                         #shownregion_kb = 400,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.6, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 7){
    # dummy2 = input$myGenes
    # dummy2 = dummy2[dummy2$biotype=="protein_coding",]
    # input$myGenes = dummy2
    # dummy = input$myLocus
    # dummy = dummy[dummy$POS>=20829710,]
    # input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = T,
                         legendpos = "topleft",
                         #shownregion_kb = myRow$region,
                         #shownregion_kb = 400,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=0.6, title_size=1.5, subtitle_size = 0.8)
  }else if(i == 10){
    dummy2 = input$myGenes
    dummy2 = dummy2[dummy2$biotype=="protein_coding",]
    input$myGenes = dummy2
    # dummy = input$myLocus
    # dummy = dummy[dummy$POS>=20829710,]
    # input$myLocus = dummy
    customASplot_woBROAD(locus = input$myLocus,
                         lead_snp = input$leadsnp, 
                         map = input$myMap,
                         genes = input$myGenes,
                         center_lead_snp = F,
                         dolegend = T,
                         legendpos = "topleft",
                         shownregion_kb = myRow$region,
                         maintitle = mainText,
                         subtitle = subText,
                         weakR2 = 0.1,
                         cex_genname=1, title_size=1.5, subtitle_size = 0.8)
  }

}
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
