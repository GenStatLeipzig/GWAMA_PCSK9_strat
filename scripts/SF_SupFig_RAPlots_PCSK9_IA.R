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
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")
source(paste0(basicpath,"/07_programme/rtools/1404_regional_assoc_plot_WO_broad/RegAssocPlot_hg19_EnsemblGenes_noBROAD_9.10_variablegengroesse.R"))

setwd(paste0(projectpath_main,"scripts/"))

#' # Load data ####
#' ***
statistics = list.files("../data/")

dumTab = foreach(i=1:8)%do%{
  #i=1 
  data = fread(paste0("../data/",statistics[i]))
  data = data[chr==1,]
  data = data[bp_hg19<=55505647 + 500000,]
  data = data[bp_hg19>=55505647 - 500000,]
  data
}
data_GWAS = rbindlist(dumTab)

myTab = dcast(data_GWAS,formula = markername ~ phenotype, 
              value.var = c("beta", "SE", "pval","invalidAssocs"))
myTab = myTab[complete.cases(myTab), ]

#' # Interaction Tests ####
#' ***
dumTab = foreach(i=1:dim(myTab)[1])%do%{
  #i=2
  IATest1 =interactionTest(mean1 = myTab$beta_PCSK9_females[i], se1 = myTab$SE_PCSK9_females[i],
                           mean2 = myTab$beta_PCSK9_males[i], se2 = myTab$SE_PCSK9_males[i])
  IATest2 =interactionTest(mean1 = myTab$beta_PCSK9_free[i], se1 = myTab$SE_PCSK9_free[i],
                           mean2 = myTab$beta_PCSK9_treated[i], se2 = myTab$SE_PCSK9_treated[i])
  
  IATest3 =interactionTest(mean1 = myTab$beta_PCSK9_females_free[i],
                           se1 = myTab$SE_PCSK9_females_free[i],
                           mean2 = myTab$beta_PCSK9_males_free[i], 
                           se2 = myTab$SE_PCSK9_males_free[i])
  IATest4 =interactionTest(mean1 = myTab$beta_PCSK9_females_treated[i], 
                           se1 = myTab$SE_PCSK9_females_treated[i],
                           mean2 = myTab$beta_PCSK9_males_treated[i], 
                           se2 = myTab$SE_PCSK9_males_treated[i])
  IATest5 =interactionTest(mean1 = myTab$beta_PCSK9_females_free[i], 
                           se1 = myTab$SE_PCSK9_females_free[i],
                           mean2 = myTab$beta_PCSK9_females_treated[i], 
                           se2 = myTab$SE_PCSK9_females_treated[i])
  IATest6 =interactionTest(mean1 = myTab$beta_PCSK9_males_free[i], 
                           se1 = myTab$SE_PCSK9_males_free[i],
                           mean2 = myTab$beta_PCSK9_males_treated[i], 
                           se2 = myTab$SE_PCSK9_males_treated[i])
  
  IATest = rbind(IATest1,IATest2,IATest3,IATest4,IATest5,IATest6)
  IATest[,type := c("females vs males",
                    "free vs treated",
                    "females vs males (free)",
                    "females vs males (treated)",
                    "free vs treated (females)",
                    "free vs treated (males)")]
  
  IATest[,type2 := c("sex",
                    "statin",
                    "sex",
                    "sex",
                    "statin",
                    "statin")]
  
  IATest[,snp := myTab$markername[i]]
  setorder(IATest,meandiff_p)
  IATest_Dummy = IATest[!duplicated(type2)]
  IATest = IATest[1,]
  IATest[,type3 := 0]
  IATest[IATest_Dummy[type2=="statin",meandiff_p <0.05],type3 := 1]
  IATest[IATest_Dummy[type2=="sex",meandiff_p <0.05],type3 := 2]
  IATest[IATest_Dummy[type2=="sex",meandiff_p <0.05] & IATest_Dummy[type2=="statin",meandiff_p <0.05],type3 := 3]
  
  IATest
}

myBestIAs = rbindlist(dumTab)

matched = match(myBestIAs$snp,data_GWAS$markername)
myBestIAs[,position := data_GWAS[matched,bp_hg19]]
setnames(myBestIAs,"meandiff_p","pvalue")
myBestIAs[,chr :=1]
setorder(myBestIAs,position)
table(myBestIAs$type,myBestIAs$pvalue<0.05)
table(myBestIAs$type3,myBestIAs$pvalue<0.05)

#' # RA Plotting ####
#' ***
# create input
myBestIAs[pvalue == min(pvalue),]

input = createInputfilesRegAssocPlot1kg19Ensembl(myBestIAs, 
                                                 doplink=T, 
                                                 r_on_server = T,
                                                 gene_validation_level = c("KNOWN"),
                                                 path_ldreference = path_1000Genomes, 
                                                 path_ldreference_snps = path_1000Genomes_snps_fn,
                                                 leadsnp = "rs28385704:55506103:C:G",
                                                 resultdir = "../temp/input_RA_PCSK9_IA/")

par(mfrow=c(1,1))
dummy = input$myLocus
dummy = dummy[dummy$POS>=55443166,]
dummy$RSQR = 0.05
input$myLocus = dummy

customASplot_woBROAD(locus = input$myLocus,
                     lead_snp = input$leadsnp, 
                     map = input$myMap,
                     genes = input$myGenes,
                     center_lead_snp = F,
                     dolegend = F, 
                     #legendpos = "topright",
                     #shownregion_kb = myRow$region,
                     shownregion_kb = 250,
                     maintitle = "Minimal interaction pvalue per SNP",
                     subtitle = "best SNP: rs28385704 (MAF: 0.126)",
                     weakR2 = 0.1,col_recombirate="red",plot_recrate_ontop=T,
                     cex_genname=1, title_size=1.5, subtitle_size = 0.8)

points(myBestIAs[type3==1 & position %in% dummy$POS, position],
       myBestIAs[type3==1 & position %in% dummy$POS, -log10(pvalue)],
       pch=21,cex=1.1,lwd=2, col="darkorange")
points(myBestIAs[type3==2 & position %in% dummy$POS, position],
       myBestIAs[type3==2 & position %in% dummy$POS, -log10(pvalue)],
       pch=21,cex=1.1,lwd=2, col="darkseagreen")
points(myBestIAs[type3==3 & position %in% dummy$POS, position],
       myBestIAs[type3==3 & position %in% dummy$POS, -log10(pvalue)],
       pch=21,cex=1.1,lwd=2, col="blue3")


tiff(filename = "../figures/SupplementalFigure_RAPlots_PCSK9_Interaction.tiff", 
     width = 1800, height = 1200, res = 200, compression = 'lzw')

par(mfrow=c(1,1))
customASplot_woBROAD(locus = input$myLocus,
                     lead_snp = input$leadsnp, 
                     map = input$myMap,
                     genes = input$myGenes,
                     center_lead_snp = F,
                     dolegend = F, 
                     #legendpos = "topright",
                     #shownregion_kb = myRow$region,
                     shownregion_kb = 250,
                     maintitle = "Minimal interaction pvalue per SNP",
                     subtitle = "best SNP: rs28385704 (MAF: 0.126)",
                     weakR2 = 0.1,col_recombirate="red",plot_recrate_ontop=T,
                     cex_genname=1, title_size=1.5, subtitle_size = 0.8)

points(myBestIAs[type3==1 & position %in% dummy$POS, position],
       myBestIAs[type3==1 & position %in% dummy$POS, -log10(pvalue)],
       pch=21,cex=1.1,lwd=2, col="darkorange")
points(myBestIAs[type3==2 & position %in% dummy$POS, position],
       myBestIAs[type3==2 & position %in% dummy$POS, -log10(pvalue)],
       pch=21,cex=1.1,lwd=2, col="darkseagreen")
points(myBestIAs[type3==3 & position %in% dummy$POS, position],
       myBestIAs[type3==3 & position %in% dummy$POS, -log10(pvalue)],
       pch=21,cex=1.1,lwd=2, col="blue3")


dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
