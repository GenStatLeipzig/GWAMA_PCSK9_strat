#' ---
#' title: "Get Supplemental Tables"
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
#' **Supplemental Tables**
#' 
#' 1) Description of Studies (as received from participating studies) **--> not included here**
#' 2) Sample Sizes, SNP Numbers, genomic inflation factor $\lambda$, and LDSC heritability results per phenotype
#' 3) Overview of independent SNPs
#' 4) Interaction analyses for the lead SNPs
#' 5) Genetic correlation
#' 6) Step80: summary of all SNPs with at least p<1e-6 in one trait 
#'    a) toploci
#'    b) GWAS Catalog annotation
#'    c) eQTL annotation
#'    d) Gene annotation
#' 7) Step80: CS 99 annotation for PCSK9 locus only!
#' 8) Coloc results
#'    a) males vs females and statin vs no statin (including cond stats for PCSK9 locus)
#'    b) PCSK9 vs eQTLs (including cond stats for PCSK9 locus)
#'    c) PCSK9 vs lipids / other traits
#' 9) MR: Bidirectional results of PCSK9 <--> LDLC (sex-strat) (allele scores?)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(1:9)),
                    Title = c("Study descriptions",
                              "Sample sizes, SNP Numbers, inflation factor and LDSC results per phenotype and setting",
                              "Interaction test results comparing SNP effects of males and females, and of statin-free and statin-treated individuals",
                              "Genetic correlation estimated by LDSC and GLGC data",
                              "Annotation of all SNPs reaching at least suggestive significance in one trait",
                              "Overview of independent SNPs",
                              "Summary of all SNPs within the 99% Credible Sets at the PCSK9 gene locus",
                              "Co-localization results",
                              "Mendelian Randomization results"),
                    Source = c("done manually",
                               "step15 of GWAS pipeline|../results/08_LDSC_Heritability.RData",
                               "../results/03_InteractionTests.RData",
                               "../results/04_LDSC_GenCor.RData",
                               "step80 of GWAS pipeline",
                               "../results/02_LociOverallPhenotypes.RData",
                               "step80 of GWAS pipeline for CS99 SNPs only",
                               "../results/06_",
                               "../results/07_"
                    ))
  
  tab0
  
}

#' # Get Sup Tab 2 ####
#' ***
#' Sample Sizes, SNP Numbers, genomic inflation factor $\lambda$, and LDSC heritability results per phenotype
#' 
#' 
{
  tab2 = fread("../../2210_GWAMA/06_Annotation/results/synopsis/step15_snp_statistics_incl_lambdas.txt")
  tab2 = tab2[,c(1,9,11:14)]
  names(tab2) = c("phenotype","nSNPs","nStudies_max","nSamples_min","nSamples_max","LambdaGC_all")
  
  load("../results/04_LDSC_Heritability.RData")
  heritab = heritab[grepl("PCSK9",trait)]
  heritab = heritab[,c(1,4,5,2:3,10,6:9)]
  matched = match(tab2$phenotype,heritab$trait)
  heritab = heritab[matched]
  stopifnot(heritab$trait == tab2$phenotype)
  tab2 = cbind(tab2, heritab[,2:10])
  names(tab2)[7] = "LambdaGC_LDSC"
  tab2
  tab2_annot = data.table(column = names(tab2),
                          description = c("Analyzed phenotype and setting",
                                          "Number of SNPs after SNP-QC (number of valid SNPs after meta-analysis)",
                                          "Maximal number of participating studies",
                                          "Minimal number of individuals",
                                          "Maximal number of individuals",
                                          "Genomic inflation factor lambda using all valid SNPs (median((beta/se)^2)/qchisq(0.5, 1))",
                                          "Genomic inflation factor lambda using only SNPs used in LDSC approach",
                                          "Mean chi-square statistic in LDSC approach",
                                          "Heritability estimate in LDSC approach (slope)",
                                          "Standard error of heritability estimate",
                                          "P-value of heritability estimate",
                                          "LDSC intercept",
                                          "Standard error of LDSC intercept",
                                          "Proportion of the inflation in the mean chi^2 that the LDSC intercept ascribes to causes other than polygenic heritability (Ratio = (intercept-1)/(mean(chi^2)-1))",
                                          "Standard error of ratio"))
  
  
}

#' # Get Sup Tab 3 ####
#' ***
#' Interaction test results
#' 
{
  loaded1 = load("../results/03_InteractionTests_3way.RData")
  loaded2 = load("../results/03_InteractionTests_2way.RData")
  
  names(IATab_2way)
  names(IATab_3way)
  
  tab3 = copy(IATab_2way)
  tab3 = tab3[,c(1:5,7,8, 10,11,9, 24:30, 12:23)]
  
  matched = match(IATab_3way$markername,tab3$markername)
  table(is.na(matched))
  IATab_3way[, bestPheno:= tab3[matched, bestPheno]]
  IATab_3way[,type := "3way"]
  IATab_3way[,fix := "none"]
  IATab_3way[,cor := NA]
  
  tab3 = rbind(tab3,IATab_3way,use.names = T, fill = T)
  setnames(tab3, "type","test")
  tab3[test == "sexIA",test := "2way_sex"]
  tab3[test == "statinIA",test := "2way_statin"]
  setnames(tab3, "type2","IA_type")
  
  tab3[IA_hierarch_fdr5proz==T & test == "2way_sex" & trait1_pval<1e-6 & trait2_pval>0.01, IA_type:="male-specific"]
  tab3[IA_hierarch_fdr5proz==T & test == "2way_sex" & trait2_pval<1e-6 & trait1_pval>0.01, IA_type:="female-specific"]
  
  tab3[IA_hierarch_fdr5proz==T & test == "2way_statin" & trait1_pval<1e-6 & trait2_pval>0.01, IA_type:="treated-specific"]
  tab3[IA_hierarch_fdr5proz==T & test == "2way_statin" & trait2_pval<1e-6 & trait1_pval>0.01, IA_type:="free-specific"]
  
  tab3[IA_hierarch_fdr5proz==T & test == "3way", IA_type:="female-treated-specific"]
  
  tab3[IA_hierarch_fdr5proz==F & IA_pval<0.05, IA_type:="nom. sig. interaction"]
  
  names(tab3)
  tab3_annot = data.table(column = gsub("trait1","traitX",names(tab3)[1:23]),
                          description = c("SNP ID, naming according to 1000 Genomes Phase 3",
                                          "Chromosome of SNP",
                                          "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                          "Effect allele",
                                          "Other allele",
                                          "Proposed canidate gene",
                                          "Best-associated PCSK9 trait",
                                          "Type of interaction test, either 3-way interaction of both sex and statin, or 2-way interaction for sex and statin",
                                          "Fixed parameter in the 2-way interaction",
                                          "Pearsons correlation between the two tested traits in the 2-way interaction test",
                                          "Difference of effect estimates. In 3-way interaction: diff = (trait2_beta - trait1_beta) - (trait4_beta - trait3_beta). In 2-way interaction: trait2 - trait1.",
                                          "Standard error of difference, defined as sqaured root of standard errors of all included traits.",
                                          "Z-Statistic of difference, follow a normal distribution under null hypothesis of no interaction",
                                          "P-value of difference",
                                          "Adjusted p-value, using hierarchical FDR correction",
                                          "TRUE/FALSE indicator for significance after hier. FDR correction",
                                          "Type of interaction",
                                          "Name of trait X",
                                          "Effect allele frequency of trait X",
                                          "Beta estimate of trait X",
                                          "Standard error of trait X",
                                          "P-value of effect of trait X",
                                          "Sample size for trait X"))
}


#' # Get Sup Tab 4 ####
#' ***
#' Genetic correlation estimated by LDSC
#' 
{
  tab4 = load("../results/04_LDSC_GeneticCorrelation.RData")
  tab4 = get(tab4)
  
  tab4_annot = data.table(column = names(tab4),
                          description = c("First phenotype",
                                          "Second phenotype",
                                          "Genetic correlation",
                                          "Standard error of genetic correlation",
                                          "Z-Score of genetic correlation",
                                          "P-value of genetic correlation"))

}

#' # Get Sup Tab 5 ####
#' ***
#' Step80: summary of all SNPs with at least p<1e-6 in one trait
#' 
#' tab5_a: toplist
#' 
#' tab5_b: GWAS catalog look-up
#' 
#' tab5_c: eQTL look-up
#' 
#' tab5_d: Nearby genes

{
  tab5_a = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/topliste_2023-01-26_PCSK9_sex_strat_v2.txt")
  names(tab5_a)
  
  col1 = grep("score_PCSK9",names(tab5_a))
  col2 = grep("I2_PCSK9",names(tab5_a))
  
  myNames1 = c("markername","cyto","chr","pos","effect_allele","other_allele","tagger","r2_tagger","tagsnp",
               "nearestgenes","gene_biotype","nearestgene","CADD_scaled","corinfo_gwas2","cisgene","transgene", "KEGG", 
               "toppheno","topmaf","topeaf","topinfo","topn","numberStudies",names(tab5_a)[col1],names(tab5_a)[col2])
  myNames1 = myNames1[c(1:26,48,27:29,49,30:32,50,33:35,51,36:38,52,39:41,53,42:44,54,45:47,55)]
  
  colsOut<-setdiff(colnames(tab5_a),myNames1)
  tab5_a[,get("colsOut"):=NULL]
  setcolorder(tab5_a,myNames1)
  
  names(tab5_a) = gsub("_score_","_",names(tab5_a))
  
  tab5a_annot = data.table(column = gsub("_PCSK9_.*","",names(tab5_a)[1:27]),
                          description = c("SNP ID, naming according to 1000 Genomes Phase 3",
                                          "Cytoband with 850 resolution",
                                          "Chromosome of SNP",
                                          "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                          "Allele for which effect sizes were calculated",
                                          "Other allele",
                                          "SNP that tags the marker in column markername (LD r^2 >=0.1)",
                                          "LD r-square value between markername and tagSNP",
                                          "TRUE/FALSE flag indicating if SNP is identical with markername",
                                          "HGNC (Ensembl) symbols of the nearest genes as specified in the settings with information on distance to the marker in SNP including functional relevance if within or within the flanking 5 kb of a gene",
                                          "Functional relevance of the gene and its validation level  for proximate genes(Ensembl)",
                                          "HGNC (Ensembl) full name of nearest gene that actually has got a full name",
                                          "PHRED-like (-10*log10(rank/total)) scaled C-score ranking a variant relative to all possible substitutions of the human genome (8.6x10^9). Scaled C-score of greater or equal 10 indicates that these are predicted to be the 10% most deleterious substitutions that you can do to the human genome, a score of greater or equal 20 indicates the 1% most deleterious and so on",
                                          "r-square value between marker in markername and a hit from GWAS catalog for named phenotype",
                                          "cis-eQTL genes for which reported eQTL-SNP has specified min R2 with marker in markername",
                                          "trans-eQTL genes for which reported eQTL-SNP has specified min R2 with marker in markername",
                                          "Analysis for nominally significant enrichment of genes in  nearest genes, cis-eQTL genes and trans-eQTL genes according to pathway enrichment specific settings in KEGG pathways",
                                          "Best associated phenotype",
                                          "Minor allele frequency for best associated phenotype (n-weighted)",
                                          "Effect allele frequency for best associated phenotype (n-weighted)",
                                          "Imputation info score (n-weighted)", 
                                          "Sample size for the best associated phenotype",
                                          "Number of studies included for this SNP",
                                          "Effect estimate",
                                          "Standard error for effect estimate",
                                          "log10-transformed p-value",
                                          "Heterogenetity I-squared"))

}

{
  tab5_b = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/gwasinfo_2023-01-26_PCSK9_sex_strat_v2.txt")
  names(tab5_b)
  tab5_b = tab5_b[snps %in% tab5_a[tagsnp==T,markername]]
  length(unique(tab5_b$pmid_gwas))
  names(tab5_b)[1] = "markername"
  names(tab5_b)[10] = "phenotype_gwas_details"
  tab5_b = tab5_b[,snp_vereinzelt_gwas := NULL]
  
  tab5b_annot = data.table(column = names(tab5_b),
                           description = c("SNP ID, naming according to 1000 Genomes Phase 3",
                                           "Chromosome of SNP",
                                           "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                           "Cytoband with 850 resolution",
                                           "LD r-square value between markername and GWAS Catalog SNP",
                                           "Distance between markername and GWAS Catalog SNP",
                                           "Reference used for LD-calculation",
                                           "GWAS Catalog SNP",
                                           "Phenotype as reported in GWAS Catalog",
                                           "Details on phenotype as reported in GWAS Catalog (in case of multiple traits in one publication)",
                                           "P-value  as reported in GWAS Catalog",
                                           "Reported genes in GWAS Catalog",
                                           "First author of publication",
                                           "Date of publication",
                                           "PMID of publication"))
}

{
  tab5_c = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/eqtlinfo_2023-01-26_PCSK9_sex_strat_v2.txt")
  names(tab5_c)
  tab5_c = tab5_c[snps %in% tab5_a[tagsnp==T,markername]]
  length(unique(tab5_c$PMID_study))
  names(tab5_c)[1] = "markername"
  
  tab5_c[,study:=NULL]
  tab5_c[,top_r2_or_p:=NULL]
  
  tab5c_annot = data.table(column = names(tab5_c),
                           description = c("SNP ID, naming according to 1000 Genomes Phase 3",
                                           "Chromosome of SNP",
                                           "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                           "LD r-square value between markername and eQTL-SNP",
                                           "Distance between markername and eQTL-SNP",
                                           "Reference used for LD-calculation",
                                           "eQTL-SNP",
                                           "Gene associating with the SNP in the eQTL study",
                                           "Probe ID of the gene associating with the SNP in the eQTL study (if available)",
                                           "Local (cis) or distant (trans) association of eQTL-SNP and eQTL-gene",
                                           "Tissue analysed in eQTL study",
                                           "Reported p-value or reported q-value",
                                           "Either p-value or q-value, relates to column p_or_q_value",
                                           "Chromosome of eQTL-SNP",
                                           "Position of eQTL-SNP in basepairs on the chromosome according to genome build hg19",
                                           "Chromosome of eQTL-gene",
                                           "Start-Position of eQTL-gene in basepairs on the chromosome according to genome build hg19",
                                           "Stop-Position of eQTL-gene in basepairs on the chromosome according to genome build hg19",
                                           "Distance between eQTL-SNP and eQTL-gene",
                                           "First author of publication",
                                           "PMID of publication"))
}

{
  tab5_d = fread("../../2210_GWAMA/06_Annotation/results/synopsis/topliste_tabdelim/proximate_genes_2023-01-26_PCSK9_sex_strat_v2.txt")
  names(tab5_d)
  tab5_d = tab5_d[markername %in% tab5_a[tagsnp==T,markername]]
  length(unique(tab5_d$genename))
  tab5_d = tab5_d[,1:17]
  tab5_d = distinct(tab5_d)
  tab5_d[,snptype:=NULL]
  tab5_d[,formateddist:=NULL]
  
  tab5d_annot = data.table(column = names(tab5_d),
                           description = c("SNP ID, naming according to 1000 Genomes Phase 3",
                                           "Chromosome of SNP",
                                           "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                           "Distance between markername and gene",
                                           "Abbreviated genename",
                                           "Functional relevance of the gene and its validation level  for proximate genes(Ensembl)",
                                           "Entrez-gene ID",
                                           "ENSG-gene ID",
                                           "Chromosome of gene",
                                           "Start-Position of gene in basepairs on the chromosome according to genome build hg19",
                                           "Stop-Position of gene in basepairs on the chromosome according to genome build hg19",
                                           "Strand of gene",
                                           "Cytoband",
                                           "Biotype of gene",
                                           "Gene Validation level. Can be either KNOWN, NOVEL, or PREDICTED"))
}

#' # Get Sup Tab 6 ####
#' ***
#' Overview of independent SNPs
#' 
{
  tab6 = load("../results/02_LociOverallPhenotypes.RData")
  tab6 = get(tab6)
  load("../results/05_GCTA_COJO.RData")
  
  table(tab6$markername %in% tab5_a$markername)
  
  matched = match(tab6$markername, tab5_a$markername)
  tab6[,cytoband := tab5_a[matched,cyto]]
  tab6[,explVar := beta^2/(beta^2 + nSamples*SE^2)]
  
  IndepSignals[,matchSNP := SNP]
  IndepSignals[SNP=="rs73487385:62616515:T:G",matchSNP := "rs3076276:62650147:CAA:CA"]
  matched = match(tab6$markername,IndepSignals$matchSNP)
  table(is.na(matched))
  tab6[,candidateGene := IndepSignals[matched,candidateGene]]
  
  names(tab6)
  
  tab6 = tab6[,c(7,2,3,20,4:6,8:14,16:19,21,15,22)]
  
  tab6_annot = data.table(column = names(tab6),
                           description = c("Chromosome of locus",
                                           "Start-Position of locus (hg19)",
                                           "Stop-Position of locus (hg19)",
                                           "Cytoband with 850 resolution",
                                           "Best-associated phenotype",
                                           "Other phenotypes associated with at least suggestive significance",
                                           "SNP ID, naming according to 1000 Genomes Phase 3",
                                           "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                           "Effect allele",
                                           "Other allele",
                                           "Effect allele frequency",
                                           "Imputation info score",
                                           "Samples size regarding top-phenotype",
                                           "Study number regarding top-phenotype",
                                           "Beta estimate regarding top-phenotype",
                                           "Standard error regarding top-phenotype",
                                           "P-value regarding top-phenotype",
                                           "Heterogeneity I^2 regarding top-phenotype",
                                           "Explained variance (SNP-wise heritability) regarding top-phenotype",
                                           "Number of additional SNPs associated at this locus with at least suggestive significance",
                                           "Proposed candidate gene. Not reported for loci with 1 or less additional SNPs"))
  
}

#' # Get Sup Tab 7 ####
#' ***
#' Step80: CS 99 annotation for PCSK9 locus only!
#' 
{
  tab7_0 = fread("../../2210_GWAMA/07_Annotation_CS99_PCSK9/results/synopsis/topliste_tabdelim/topliste_2023-01-26_PCSK9_sex_strat_CredSet99.txt")
  tab7_0 = tab7_0[chr==1,]
  tab7 = distinct(tab7_0)
  tab7 = copy(tab7)
  
  dumTab = data.table(groups = unique(tab7$group))
  dumTab[grepl("rs11591147",groups),SNPregion := "A"]
  dumTab[grepl("rs28385704",groups),SNPregion := "B"]
  dumTab[is.na(SNPregion),SNPregion := "C"]
  dumTab[,num := gsub("::.*","",groups)]
  matched = match(dumTab$num,IndepSignals$num)
  dumTab[,phenotype := IndepSignals[matched,pheno]]
  dumTab[,phenotype2 := gsub("PCSK9_","",phenotype)]
  dumTab[,group2 := paste0("Region ",SNPregion,": ",phenotype2)]
  
  matched = match(tab7$group,dumTab$groups)
  table(tab7$group == dumTab[matched,groups])
  tab7[,indepSNP := gsub(".*::","",group)]
  tab7[,group := dumTab[matched,group2]]
  
  names(tab7)
  
  tab7_help = foreach(i = 1:dim(dumTab)[1])%do%{
    #i=1
    myRegion = dumTab[i,]
    tab7_2 = copy(tab7)
    tab7_2 = tab7_2[group==myRegion$group2,]
    
    myNames1 = names(tab7_2)[c(2,127,9,12,15,61,62,3:4,24)]
    myNames2 = paste0(c("beta_score_","se_score_","logp_score_","I2_"),myRegion$phenotype)
    myNames = c(myNames1,myNames2)
   
    colsOut<-setdiff(colnames(tab7_2),myNames)
    tab7_2[,get("colsOut"):=NULL]
    setcolorder(tab7_2,myNames1)
    myNames3 = c("group","indepSNP","markername","chr", "pos","effect_allele","other_allele","PostProb","SumProb","CADD_scaled",
                 "beta_GWAMA","se_GWAMA","logp_GWAMA","I2_GWAMA")
    names(tab7_2) = myNames3
    
    SNP = unique(tab7_2$indepSNP)
    if(grepl("treated",myRegion$phenotype)==F){
      cojoTab = fread(paste0("../results/05_GCTA_COJO_cond/PCSK9::",myRegion$phenotype,"_signal_",SNP,".cma.cojo"))
      matched = match(tab7_2$markername,cojoTab$SNP)
      stopifnot(is.na(matched)==0)
      cojoTab = cojoTab[matched]
      tab7_2[,beta_COJO := cojoTab$bC]
      tab7_2[,se_COJO := cojoTab$bC_se]
      tab7_2[,logp_COJO := -log10(cojoTab$pC)]
    }
    tab7_2
  }
  tab7_help = rbindlist(tab7_help,use.names = T,fill=T)
  tab7 = tab7_help
  
  tab7_annot = data.table(column = names(tab7),
                          description = c("Region of the PCSK9 locus, defined by independent SNPs across the associated traits",
                                          "rsID of the independent SNP",
                                          "SNP ID, naming according to 1000 Genomes Phase 3",
                                          "Chromosome of SNP",
                                          "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                          "Effect allele",
                                          "Other allele",
                                          "Posterior probability of conditioned statistics",
                                          "Cumulative sum of posterior probability per region and trait",
                                          "PHRED-like (-10*log10(rank/total)) scaled C-score ranking a variant relative to all possible substitutions of the human genome (8.6x10^9). Scaled C-score of greater or equal 10 indicates that these are predicted to be the 10% most deleterious substitutions that you can do to the human genome, a score of greater or equal 20 indicates the 1% most deleterious and so on",
                                          "Beta estimate of GWAMA regarding phenotype defined in column group",
                                          "Standard error of GWAMA regarding phenotype defined in column group",
                                          "-log10 transformed p-value of GWAMA regarding phenotype defined in column group",
                                          "Heterogeneity I^2 of GWAMA regarding phenotype defined in column group",
                                          "Beta estimate of conditional analyses regarding phenotype defined in column group",
                                          "Standard error of conditional analyses regarding phenotype defined in column group",
                                          "-log10 transformed p-value of conditional analyses  regarding phenotype defined in column group"))
}

#' # Get Sup Tab 8 ####
#' ***
#' Co-localization results
#' 
#'    a) males vs females and statin vs no statin (including cond stats for PCSK9 locus)
#'    b) PCSK9 vs eQTLs (including cond stats for PCSK9 locus)
#'    c) PCSK9 vs lipids / other traits
{
  loaded = load("../results/06_2_coloc_eQTLs.RData")
  coloc_eQTLs_uncond = get(loaded[1])
  loaded = load("../results/06_3_coloc_eQTLs_cond.RData")
  coloc_eQTLs_cond = get(loaded[1])
  
  loaded = load("../results/06_4_coloc_withinPCSK9.RData")
  coloc_within_uncond = get(loaded[1])
  loaded = load("../results/06_5_coloc_withinPCSK9_cond.RData")
  coloc_within_cond = get(loaded[1])
  
  loaded = load("../results/06_6_coloc_otherGWAS.RData")
  coloc_otherTraits = get(loaded[1])
  loaded = load("../results/06_7_coloc_otherGWAScond.RData")
  coloc_otherTraits_cond = get(loaded[1])
  
  # tab 8 a: tests within our data
  names(coloc_within_uncond)
  names(coloc_within_cond)
  coloc_within_cond[,NR:=NULL]
  matched = match(coloc_within_uncond$gene,tab6$candidateGene)
  table(is.na(matched))
  coloc_within_uncond[,locus:=tab6[matched,cytoband]]
  coloc_within_uncond[,sex1:=unlist(strsplit(trait1,"_"))[seq(2,90,3)]]
  coloc_within_uncond[,sex2:=unlist(strsplit(trait2,"_"))[seq(2,90,3)]]
  coloc_within_uncond[,statin1:=unlist(strsplit(trait1,"_"))[seq(3,90,3)]]
  coloc_within_uncond[,statin2:=unlist(strsplit(trait2,"_"))[seq(3,90,3)]]
  
  coloc_within_uncond[sex1 == sex2,fixed := "sex"]
  coloc_within_uncond[statin1 == statin2,fixed := "statin"]
  
  coloc_within_cond[1:8,fixed := "sex"]
  coloc_within_cond[9:12,fixed := "statin"]
  
  coloc_within_uncond[,sex1:=NULL]
  coloc_within_uncond[,sex2:=NULL]
  coloc_within_uncond[,statin1:=NULL]
  coloc_within_uncond[,statin2:=NULL]
  
  tab8_a = rbind(coloc_within_cond,coloc_within_uncond)
  
  # tab 8 b: tests with eQTLs from GTEx
  names(coloc_eQTLs_uncond)
  names(coloc_eQTLs_cond)
  setnames(coloc_eQTLs_uncond,"cytoband","locus")
  
  tab8_b = rbind(coloc_eQTLs_cond,coloc_eQTLs_uncond)
  setnames(tab8_b,"gene","trait2_gene")
  setnames(tab8_b,"trait2","trait2_tissue")
  tab8_b[,trait2_tissue := gsub("GE in ","",trait2_tissue)]
  
  # tab 8 c: tests with other GWAS
  names(coloc_otherTraits)
  names(coloc_otherTraits_cond)
  
  tab8_c = rbind(coloc_otherTraits_cond,coloc_otherTraits)
  names(tab8_c)
  matched = match(tab8_c$gene,tab6$candidateGene)
  table(is.na(matched))
  tab8_c[,locus:=tab6[matched,cytoband]]
  tab8_c = tab8_c[,c(10,1:9)]
  
  # checks
  names(tab8_a)
  names(tab8_b)
  names(tab8_c)
  
  # annotation 
  myNames = unique(c(names(tab8_a),names(tab8_b),names(tab8_c)))
  myNames = myNames[c(1:5,12,13,6:11)]
  tab8_annot = data.table(column = myNames,
                           description = c("Genomic cytoband of index SNP",
                                           "Proposed canidate gene",
                                           "In Table 8a: fixed strata for this test",
                                           "First trait to be tested, always one of our PCSK9 traits",
                                           "Second trait to be tested, either one of our PCSK9 traits (Tab 8a) or other GWAS trait",
                                           "In Table 8b: Second trait to be tested, respective gene",
                                           "In Table 8b: Second trait to be tested, respective tissue",
                                           "Number of SNPs included in co-localization analysis per test",
                                           "Posterior probability for hypothesis 0: neither trait associated",
                                           "Posterior probability for hypothesis 1: only trait 1 associated",
                                           "Posterior probability for hypothesis 2: only trait 2 associated",
                                           "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                           "Posterior probability for hypothesis 4: both trait associated, shared signal"))
}

#' # Get Sup Tab 9 ####
#' ***
#' Bidirectional Mendelian Randomization results
{
  loaded = load(file="../results/07_MR_PCSK9_LDLC_SNPwise.RData")
  MR_PCSK9_LDLC_SNPwise = get(loaded)
  loaded = load(file="../results/07_MR_PCSK9_LDLC_summary.RData")
  MR_PCSK9_LDLC_summary = get(loaded)
  
  loaded = load(file="../results/07_MR_LDLC_PCSK9_SNPwise.RData")
  MR_LDLC_PCSK9_SNPwise = get(loaded)
  loaded = load(file="../results/07_MR_LDLC_PCSK9_summary.RData")
  MR_LDLC_PCSK9_summary = get(loaded)

  # tab 9a: MR results per SNP
  names(MR_PCSK9_LDLC_SNPwise)
  names(MR_LDLC_PCSK9_SNPwise)
  
  MR_PCSK9_LDLC_SNPwise[,P.LDLC := beta.LDLC/SE.LDLC]
  MR_LDLC_PCSK9_SNPwise[,logP.LDLC := beta.LDLC/SE.LDLC]
  setnames(MR_PCSK9_LDLC_SNPwise,"P.LDLC","ZScore.LDLC")
  setnames(MR_LDLC_PCSK9_SNPwise,"logP.LDLC","ZScore.LDLC")
  
  MR_PCSK9_LDLC_SNPwise[,exposure := phenotype]
  MR_PCSK9_LDLC_SNPwise[,outcome := "LDLC_all"]
  MR_PCSK9_LDLC_SNPwise[grepl("_fe",phenotype),outcome := "LDLC_female"]
  MR_PCSK9_LDLC_SNPwise[grepl("_ma",phenotype),outcome := "LDLC_male"]
  MR_PCSK9_LDLC_SNPwise[,phenotype := NULL]
  MR_PCSK9_LDLC_SNPwise[,EA := NULL]
  MR_PCSK9_LDLC_SNPwise = MR_PCSK9_LDLC_SNPwise[,c(19,20,1:18)]
  
  MR_LDLC_PCSK9_SNPwise[,exposure := "LDLC_all"]
  MR_LDLC_PCSK9_SNPwise[grepl("_fe",phenotype),exposure := "LDLC_female"]
  MR_LDLC_PCSK9_SNPwise[grepl("_ma",phenotype),exposure := "LDLC_male"]
  MR_LDLC_PCSK9_SNPwise[,outcome := phenotype]
  MR_LDLC_PCSK9_SNPwise[,phenotype := NULL]
  MR_LDLC_PCSK9_SNPwise[,EA.LDLC := NULL]
  MR_LDLC_PCSK9_SNPwise = MR_LDLC_PCSK9_SNPwise[,c(19,20,1:18)]
  
  tab9_a = rbind(MR_PCSK9_LDLC_SNPwise,MR_LDLC_PCSK9_SNPwise)
  tab9_a
  
  # tab 9b: MR results combined
  names(MR_PCSK9_LDLC_summary)
  names(MR_LDLC_PCSK9_summary)
  
  MR_PCSK9_LDLC_summary[,exposure := phenotype]
  MR_PCSK9_LDLC_summary[,outcome := "LDLC_all"]
  MR_PCSK9_LDLC_summary[grepl("_fe",phenotype),outcome := "LDLC_female"]
  MR_PCSK9_LDLC_summary[grepl("_ma",phenotype),outcome := "LDLC_male"]
  MR_PCSK9_LDLC_summary[,phenotype := NULL]
  MR_PCSK9_LDLC_summary = MR_PCSK9_LDLC_summary[,c(9,2,1,3:8)]
  
  MR_LDLC_PCSK9_summary[,outcome := phenotype]
  MR_LDLC_PCSK9_summary[,exposure := "LDLC_all"]
  MR_LDLC_PCSK9_summary[grepl("_fe",phenotype),exposure := "LDLC_female"]
  MR_LDLC_PCSK9_summary[grepl("_ma",phenotype),exposure := "LDLC_male"]
  MR_LDLC_PCSK9_summary[,phenotype := NULL]
  MR_LDLC_PCSK9_summary = MR_LDLC_PCSK9_summary[,c(2,9,1,3:8)]
  
  tab9_b = rbind(MR_PCSK9_LDLC_summary,MR_LDLC_PCSK9_summary)
  tab9_b
  
  # annotation 
  tab9a_annot = data.table(column = names(tab9_a),
                          description = c("Exposure used in MR",
                                          "Outcome used in MR",
                                          "SNP ID, naming according to 1000 Genomes Phase 3",
                                          "Chromosome of SNP",
                                          "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                          "Effect allele frequency in PCSK9",
                                          "Beta estimate in PCSK9",
                                          "Standard error in PCSK9",
                                          "P-value in PCSK9",
                                          "Sample size in PCSK9",
                                          "Effect allele frequency in LDLC",
                                          "Beta estimate in LDLC",
                                          "Standard error in LDLC",
                                          "-log10 transformed p-value in LDLC",
                                          "Sample size in LDLC",
                                          "Ratio of beta estimates (outcome / exposure)",
                                          "Standard error of ratio, first term of delta method expansion",
                                          "Standard error of ratio, first two terms of delta method expansion",
                                          "P-value for ratio, using 'SEst.Ratio'",
                                          "P-value for ratio, using 'SEnd.Ratio'"))
  tab9b_annot = data.table(column = names(tab9_b),
                           description = c("Exposure used in MR",
                                           "Outcome used in MR",
                                           "Number of SNPs used in the inverse-variance weighted MR",
                                           "Beta estimate of inverse-variance weighted MR",
                                           "Standard error of 'beta_IV'",
                                           "P-value of 'beta_IV'",
                                           "Cochrans Q of 'beta_IV'",
                                           "P-value of Cochrans Q",
                                           "Comment to SNP selection"))
}


#' # Save tables ###
#' ***
tosave4 = data.table(data = c("tab0", "tab2","tab6","tab3","tab4","tab5_a","tab5_b","tab5_c","tab5_d",
                              "tab7","tab8_a","tab8_b","tab8_c","tab9_a","tab9_b"), 
                     SheetNames = c("Content","TableS2","TableS3","TableS4","TableS5","TableS6a","TableS6b","TableS6c",
                                    "TableS6d","TableS7","TableS8a","TableS8b","TableS8c","TableS9a","TableS9b"))
excel_fn = "../tables/SupplementalTables.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

tosave4 = data.table(data = c("tab2_annot","tab6_annot","tab3_annot","tab4_annot","tab5a_annot","tab5b_annot",
                              "tab5c_annot","tab5d_annot","tab7_annot","tab8_annot","tab9a_annot","tab9b_annot"), 
                     SheetNames = c("TableS2_annot","TableS3_annot","TableS4_annot","TableS5_annot","TableS6a_annot",
                                    "TableS6b_annot","TableS6c_annot","TableS6d_annot","TableS7_annot","TableS8_annot",
                                    "TableS9a_annot","TableS9b_annot"))
excel_fn = "../tables/SupplementalTables_Annotation.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab0, tab2,tab3,tab4,tab5_a,tab5_b,tab5_c,tab5_d,tab6,tab7,tab8_a,tab8_b,tab8_c,tab9_a,tab9_b,
     file = "../tables/SupplementalTables.RData")
save(tab2_annot,tab3_annot,tab4_annot,tab5a_annot,tab5b_annot,tab5c_annot,tab5d_annot,tab6_annot, 
     tab7_annot,tab8_annot,tab9a_annot,tab9b_annot,
     file = "../tables/SupplementalTables_annot.RData")

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))

