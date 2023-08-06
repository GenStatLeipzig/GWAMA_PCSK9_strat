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
#' 2) Sample Sizes, SNP Numbers, genomic inflation factor $\lambda$
#' 3) Overview of independent SNPs
#' 4) Step80: summary of all SNPs with at least p<1e-6 in one trait & at valid loci
#'    a) toploci
#'    b) GWAS Catalog annotation
#'    c) eQTL annotation
#'    d) Gene annotation
#' 5) GCTA COJO-slct results (PCSK9 only)
#' 6) Coloc results
#'    a) PCSK9 vs eQTLs (including cond stats for PCSK9 locus)
#'    b) PCSK9 vs other traits (including cond stats for PCSK9 locus)
#'    c) within PCSK9 (including cond stats for PCSK9 locus)
#' 7) Interaction results
#'    a) Lead SNPs
#'    b) PCSK9 special
#' 8) MR results
#'    a) MR of PCSK9 on LDLC
#'    b) Interaction tests of causal estimates
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(1:8)),
                    Title = c("Study descriptions",
                              "Sample sizes, SNP Numbers, and inflation factor per subgroup",
                              "Overview of all associated loci",
                              "Annotation of GWAMA results",
                              "GCTA COJO select results",
                              "Colocalization results",
                              "Interaction results",
                              "Mendelian Randomization results"),
                    Source = c("done manually",
                               "step15 of GWAS pipeline",
                               "../results/02_LociOverallPhenotypes.RData",
                               "step80 of GWAS pipeline",
                               "../results/03_GCTA_COJO.RData",
                               "../results/05_",
                               "../results/04_",
                               "../results/06_"))
  
  tab0
  
}

#' # Get Sup Tab 2 ####
#' ***
#' Sample Sizes, SNP Numbers, genomic inflation factor $\lambda$
#' 
#' 
{
  tab2 = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/step15_snp_statistics_incl_lambdas.txt")
  tab2 = tab2[,c(1,9,11:14)]
  names(tab2) = c("phenotype","nSNPs","nStudies_max","nSamples_min","nSamples_max","LambdaGC_all")
  
  tab2_annot = data.table(column = names(tab2),
                          description = c("Analyzed subgroups",
                                          "Number of SNPs after SNP-QC (number of valid SNPs after meta-analysis)",
                                          "Maximal number of participating studies",
                                          "Minimal number of individuals",
                                          "Maximal number of individuals",
                                          "Genomic inflation factor lambda using all valid SNPs (median((beta/se)^2)/qchisq(0.5, 1))"))
  
  
}

#' # Get Sup Tab 3 ####
#' ***
#' Overview all loci
{
  load("../results/02_LociOverallPhenotypes.RData")
  tab3 = copy(result.4)
  
  tab3_annot = data.table(column = names(tab3),
                          description = c("Numbered locus per chromosome",
                                          "First base position of this locus",
                                          "Last base position of this locus",
                                          "Best associated subgroup",
                                          "Other associated subgroups",
                                          "SNP ID of best associated SNP at this locus (lead SNP)",
                                          "Chromosome number",
                                          "Base position of lead SNP",
                                          "Effect alllele of lead SNP",
                                          "Other allele of lead SNP",
                                          "Effect allele frequency of lead SNP in best subgroup",
                                          "Impuation info score of lead SNP in best subgroup",
                                          "Sample size for lead SNP in best subgroup",
                                          "Number of studies available for lead SNP in best subgroup",
                                          "Number of SNPs associated at this locus with at least suggestive significance, summed up over all associated subgroups (SNPs can be counted multiple times)",
                                          "Beta estimate of lead SNP in best subgroup",
                                          "Standard error of lead SNP effect in best subgroup",
                                          "Association p-value of lead SNP in best subgroup",
                                          "Heterogeneity I2 parameter of lead SNP in best subgroup"))
  
}

#' # Get Sup Tab 4 ####
#' ***
#' Step80: summary of all SNPs with at least p<1e-6 in one trait
#' 
#' tab4_a: toplist (all SNPs at valid loci!)
#' 
#' tab4_b: GWAS catalog look-up (only for independent signals)
#' 
#' tab4_c: eQTL look-up (only for independent signals)
#' 
#' tab4_d: Nearby genes (only for independent signals)

{
  tab4a = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/topliste_2023-07-26_PCSK9_strat.txt")
  
  # check valid loci
  invalidSNPs = tab3[NR_SNPs <2,markername]
  table(is.element(tab4a$markername,invalidSNPs))
  table(is.element(tab4a$tagger,invalidSNPs))
  tab4a = tab4a[!is.element(tagger,invalidSNPs)]
  names(tab4a)
  
  col1 = grep("score_PCSK9",names(tab4a))
  col2 = grep("I2_PCSK9",names(tab4a))
  
  myNames1 = c("markername","cyto","chr","pos","effect_allele","other_allele","tagger","r2_tagger","tagsnp",
               "nearestgenes","gene_biotype","nearestgene","CADD_scaled","corinfo_gwas2","cisgene","transgene", "KEGG", 
               "toppheno","topmaf","topeaf","topinfo","topn","numberStudies",names(tab4a)[col1],names(tab4a)[col2])
  myNames1 = myNames1[c(1:26,48,27:29,49,30:32,50,33:35,51,36:38,52,39:41,53,42:44,54,45:47,55)]
  
  colsOut<-setdiff(colnames(tab4a),myNames1)
  tab4a[,get("colsOut"):=NULL]
  setcolorder(tab4a,myNames1)
  
  names(tab4a) = gsub("_score_","_",names(tab4a))
  
  tab4a_annot = data.table(column = gsub("_PCSK9_.*","",names(tab4a)[1:27]),
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
  tab4b = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/gwasinfo_2023-07-26_PCSK9_strat.txt")
  names(tab4b)
  tab4b = tab4b[snps %in% tab3[NR_SNPs>2,markername]]
  length(unique(tab4b$pmid_gwas))
  names(tab4b)[1] = "markername"
  names(tab4b)[10] = "phenotype_gwas_details"
  tab4b = tab4b[,snp_vereinzelt_gwas := NULL]
  
  tab4b_annot = data.table(column = names(tab4b),
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
  tab4c = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/eqtlinfo_2023-07-26_PCSK9_strat.txt")
  names(tab4c)
  tab4c = tab4c[snps %in% tab3[NR_SNPs>2,markername]]
  length(unique(tab4c$PMID_study))
  names(tab4c)[1] = "markername"
  
  tab4c[,study:=NULL]
  tab4c[,top_r2_or_p:=NULL]
  
  tab4c_annot = data.table(column = names(tab4c),
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
  tab4d = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/proximate_genes_2023-07-26_PCSK9_strat.txt")
  names(tab4d)
  tab4d = tab4d[markername %in% tab3[NR_SNPs>2,markername]]
  length(unique(tab4d$genename))
  tab4d = tab4d[,1:17]
  tab4d = distinct(tab4d)
  tab4d[,snptype:=NULL]
  tab4d[,formateddist:=NULL]
  
  tab4d_annot = data.table(column = names(tab4d),
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

#' # Get Sup Tab 5 ####
#' ***
#' GCTA COJO results
{
  load("../results/03_GCTA_COJO.RData")
  tab5 = copy(IndepSignals)
  names(tab5)
  tab5 = tab5[,c(15,21,26,2,1,3:14)]
  
  tab5_annot = data.table(column = names(tab5),
                          description = c("Analyzed subgroup (multiple columns per subgroup and locus possible)",
                                          "Candidate gene",
                                          "RS ID of independent signal",
                                          "SNP ID of independent signal",
                                          "Chromosome number",
                                          "Base position (hg19)",
                                          "Effect allele",
                                          "Effect allele frequency",
                                          "Beta estimate of SNP in GWAMA",
                                          "Standard error of SNP effect in GWAMA",
                                          "Association p-value of SNP in GWAMA",
                                          "Estimated ffective sample size",
                                          "Frequency of the effect allele in LIFE",
                                          "Effect size from joint analysis",
                                          "Standard error from joint analysis",
                                          "P-value from joint analysis",
                                          "LD correlation between the SNP i and SNP i + 1 for the SNPs on the list per subgroup"))
  
}

#' # Get Sup Tab 6 ####
#' ***
#' Co-localization results
#' 
#'  a) PCSK9 vs eQTLs (including cond stats for PCSK9 locus)
#'  b) PCSK9 vs lipids / other traits (including cond stats for PCSK9 locus)
#'  c) males vs females and statin vs no statin (including cond stats for PCSK9 locus)
{
  loaded = load("../results/05_2_coloc_eQTLs.RData")
  coloc_eQTLs_uncond = get(loaded[1])
  loaded = load("../results/05_3_coloc_eQTLs_cond.RData")
  coloc_eQTLs_cond = get(loaded[1])
  
  loaded = load("../results/05_4_coloc_withinPCSK9.RData")
  coloc_within_uncond = get(loaded[1])
  loaded = load("../results/05_5_coloc_withinPCSK9_cond.RData")
  coloc_within_cond = get(loaded[1])
  
  loaded = load("../results/05_6_coloc_otherGWAS.RData")
  coloc_otherTraits = get(loaded[1])
  loaded = load("../results/05_7_coloc_otherGWAScond.RData")
  coloc_otherTraits_cond = get(loaded[1])
  
  # tab 6 a: tests with eQTLs from GTEx
  names(coloc_eQTLs_uncond)
  names(coloc_eQTLs_cond)
  setnames(coloc_eQTLs_cond,"locus","cytoband")
  coloc_eQTLs_cond[,cytoband := "01p32.3"]
  coloc_eQTLs_cond[,indepSignal := gsub(".*__","",trait1)]
  coloc_eQTLs_cond[,trait1 := gsub("__.*","",trait1)]
  
  tab6a = rbind(coloc_eQTLs_cond,coloc_eQTLs_uncond,fill=TRUE)
  setnames(tab6a,"gene","trait2_gene")
  setnames(tab6a,"trait2","trait2_tissue")

  # tab 6 b: tests with other GWAS
  names(coloc_otherTraits)
  names(coloc_otherTraits_cond)
  coloc_otherTraits_cond[,gene := "PCSK9"]
  
  tab6b = rbind(coloc_otherTraits_cond,coloc_otherTraits,fill=T)
  names(tab6b)
  matched = match(tab6b$gene,tab6a$trait2_gene)
  table(is.na(matched))
  tab6b[,cytoband:=tab6a[matched,cytoband]]
  tab6b = tab6b[,c(11,2,10,3:9,1)]
  setnames(tab6b,"SNP","indepSignal")
  
  # tab 6 c: tests within our data
  names(coloc_within_uncond)
  names(coloc_within_cond)
  matched = match(coloc_within_uncond$gene,tab6a$trait2_gene)
  table(is.na(matched))
  coloc_within_uncond[,cytoband:=tab6a[matched,cytoband]]
  matched = match(coloc_within_cond$gene,tab6a$trait2_gene)
  table(is.na(matched))
  coloc_within_cond[,cytoband:=tab6a[matched,cytoband]]
  coloc_within_cond[,locus := NULL]
  
  coloc_within_uncond[,sex1:=unlist(strsplit(trait1,"_"))[seq(2,78,3)]]
  coloc_within_uncond[,sex2:=unlist(strsplit(trait2,"_"))[seq(2,78,3)]]
  coloc_within_uncond[,statin1:=unlist(strsplit(trait1,"_"))[seq(3,78,3)]]
  coloc_within_uncond[,statin2:=unlist(strsplit(trait2,"_"))[seq(3,78,3)]]
  
  coloc_within_uncond[sex1 == sex2,fixed := "sex"]
  coloc_within_uncond[statin1 == statin2,fixed := "statin"]
  
  coloc_within_cond[1:4,fixed := "statin"]
  coloc_within_cond[5:8,fixed := "sex"]
  coloc_within_cond[9:12,fixed := "statin"]
  coloc_within_cond[13:16,fixed := "statin"]
  coloc_within_cond[17:20,fixed := "sex"]
  coloc_within_cond[21:24,fixed := "sex"]
  
  coloc_within_uncond[,sex1:=NULL]
  coloc_within_uncond[,sex2:=NULL]
  coloc_within_uncond[,statin1:=NULL]
  coloc_within_uncond[,statin2:=NULL]
  
  tab6c = rbind(coloc_within_cond,coloc_within_uncond,fill=T)
  tab6c[!grepl("PCSK9",trait1),trait1 := paste0("PCSK9_",trait1)]
  tab6c[!grepl("PCSK9",trait2),trait2 := paste0("PCSK9_",trait2)]
  tab6c = tab6c[,c(12,4,3,2,5:11,1)]
  setnames(tab6c,"SNP","indepSignal")
  
  # checks
  names(tab6a)
  names(tab6b)
  names(tab6c)
  
  # annotation 
  myNames = unique(c(names(tab6a),names(tab6b),names(tab6c)))
  myNames = myNames[c(1,12,2,14,13,3,4,5:11)]
  tab6_annot = data.table(column = myNames,
                          description = c("Genomic cytoband of locus",
                                          "Proposed canidate gene",
                                          "First trait to be tested, always one of our PCSK9 traits",
                                          "In Table 6c: fixed strata for this test",
                                          "Second trait to be tested, either one of our PCSK9 traits (Tab 6c) or other GWAS trait (6c)",
                                          "In Table 6b: Second trait to be tested, respective gene",
                                          "In Table 6b: Second trait to be tested, respective tissue",
                                          "Number of SNPs included in co-localization analysis per test",
                                          "Posterior probability for hypothesis 0: neither trait associated",
                                          "Posterior probability for hypothesis 1: only trait 1 associated",
                                          "Posterior probability for hypothesis 2: only trait 2 associated",
                                          "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                          "Posterior probability for hypothesis 4: both trait associated, shared signal",
                                          "RS ID of independent signal, in case of conditional statistics used for PCSK9 trait"))

}

#' # Get Sup Tab 7 ####
#' ***
#' Interaction test results
#' 
{
  loaded1 = load("../results/04_IATest_complete.RData")
  tab7a = copy(IATab)
  loaded2 = load("../results/04_IATest_PCSK9Special_complete.RData")
  tab7b = copy(IATab)
  
  tab7a[,fix := NULL]
  tab7b[,fix := NULL]
  
  table(names(tab7a)==names(tab7b))
  myNames = names(tab7a)
  myNames = myNames[1:20]
  myNames = gsub("trait1","traitX",myNames)
  
  tab7_annot = data.table(column = myNames,
                          description = c("SNP ID, naming according to 1000 Genomes Phase 3",
                                          "Chromosome of SNP",
                                          "Position of SNP in basepairs on the chromosome according to genome build hg19",
                                          "Effect allele",
                                          "Other allele",
                                          "Proposed canidate gene",
                                          "Best-associated PCSK9 trait",
                                          "Type of interaction test, either 3-way interaction of both sex and statin, or 2-way interaction for sex and statin",
                                          "Difference of effect estimates. In 3-way interaction: diff = (trait2_beta - trait1_beta) - (trait4_beta - trait3_beta). In 2-way interaction: trait2 - trait1.",
                                          "Standard error of difference, defined as sqaured root of standard errors of all included traits.",
                                          "Z-Statistic of difference, follow a normal distribution under null hypothesis of no interaction",
                                          "P-value of difference",
                                          "Adjusted p-value, using hierarchical FDR correction",
                                          "TRUE/FALSE indicator for significance after hier. FDR correction",
                                          "Name of trait X",
                                          "Effect allele frequency of trait X",
                                          "Beta estimate of trait X",
                                          "Standard error of trait X",
                                          "P-value of effect of trait X",
                                          "Sample size for trait X"))
}

#' # Get Sup Tab 8 ####
#' ***
#' Mendelian Randomization results
{
  loaded = load(file="../results/06_MRresults.RData")
  tab8a = get(loaded)
  loaded = load(file="../results/06_MR_interaction.RData")
  tab8b = get(loaded)
  
  # tab 9a: MR results per SNP
  names(tab8a)
  names(tab8b)

  # annotation 
  tab8a_annot = data.table(column = names(tab8a),
                          description = c("Exposure used in MR",
                                          "Outcome used in MR",
                                          "SNP selection model in MR",
                                          "Number of SNPs used in the MR according to setting",
                                          "Beta estimate of inverse-variance weighted MR",
                                          "Standard error of 'beta_IVW'",
                                          "P-value of 'beta_IVW'",
                                          "Cochrans Q of 'beta_IVW'",
                                          "P-value of Cochrans Q in IVW",
                                          "Beta estimate of MR-Egger",
                                          "Standard error of 'beta_egger'",
                                          "P-value of 'beta_egger'",
                                          "Beta estimate of MR-Egger intercept",
                                          "Standard error of 'beta_egger_int'",
                                          "P-value of 'beta_egger_int'",
                                          "Cochrans Q of 'beta_egger'",
                                          "P-value of Cochrans Q in Egger"))
  
  myNames = names(tab8b)
  myNames = myNames[1:12]
  myNames = gsub("trait1","traitX",myNames)
  
  tab8b_annot = data.table(column = myNames,
                          description = c("Type of interaction test, either 3-way interaction of both sex and statin, or 2-way interaction for sex and statin",
                                          "SNP selection model in MR",
                                          "Difference of effect estimates. In 3-way interaction: diff = (trait2_beta - trait1_beta) - (trait4_beta - trait3_beta). In 2-way interaction: trait2 - trait1.",
                                          "Standard error of difference, defined as sqaured root of standard errors of all included traits.",
                                          "Z-Statistic of difference, follow a normal distribution under null hypothesis of no interaction",
                                          "P-value of difference",
                                          "Adjusted p-value, using hierarchical FDR correction",
                                          "TRUE/FALSE indicator for significance after hier. FDR correction",
                                          "Name of trait X",
                                          "Beta estimate of trait X",
                                          "Standard error of trait X",
                                          "P-value of effect of trait X"))
  
  
}


#' # Save tables ###
#' ***
tosave4 = data.table(data = c("tab0", "tab2","tab3","tab4a","tab4b","tab4c","tab4d","tab5",
                              "tab6a","tab6b","tab6c","tab7a","tab7b","tab8a","tab8b"), 
                     SheetNames = c("Content","TableS2","TableS3","TableS4a","TableS4b","TableS4c","TableS4d","TableS5",
                                    "TableS6a","TableS6b","TableS6c","TableS7a","TableS7b","TableS8a","TableS8b"))
excel_fn = "../tables/SupplementalTables_230806.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

tosave4 = data.table(data = c("tab2_annot","tab3_annot","tab4a_annot","tab4b_annot","tab4c_annot","tab4d_annot",
                              "tab5_annot","tab6_annot","tab7_annot","tab8a_annot","tab8b_annot"), 
                     SheetNames = c("TableS2_annot","TableS3_annot","TableS4a_annot","TableS4b_annot","TableS4c_annot",
                                    "TableS4d_annot","TableS5_annot","TableS6_annot","TableS7_annot","TableS8a_annot",
                                    "TableS8b_annot"))
excel_fn = "../tables/SupplementalTables_Annotation_230806.xlsx"
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab0, tab2,tab3,tab4a,tab4b,tab4c,tab4d,tab5,tab6a,tab6b,tab6c,tab7a,tab7b,tab8a,tab8b,
     file = "../tables/SupplementalTables_230806.RData")
save(tab2_annot,tab3_annot,tab4a_annot,tab4b_annot,tab4c_annot,tab4d_annot,tab5_annot,tab6_annot, 
     tab7_annot,tab8a_annot,tab8b_annot,
     file = "../tables/SupplementalTables_annot_230806.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0, tz,units = "mins"),2))

