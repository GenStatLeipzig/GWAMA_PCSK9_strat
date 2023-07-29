#' ---
#' title: "Co-localization Part 1: get eQTL data"
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
#' *Co-localization analysis of gene-expression quantitative trait loci*: 
#' 
#' Used eQTL databases:
#' 
#' * GTEx Analysis V8, 49 tissues, hg38-> liftover using the GTEx SNP Annotation file
#' 
#' Gene selection: 
#' 
#' * nearby genes (+/- 250 kb)
#' * cis eQTL genes (LD r^2>=0.3)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))

#' # Step 1: Get genes per region ####
#' ***
load("../results/03_GCTA_COJO.RData")

myPhenos = unique(IndepSignals$pheno)
mySNPs = unique(IndepSignals$SNP)

genes = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/proximate_genes_2023-07-26_PCSK9_strat.txt")
eQTLs = fread("../../2307_GWAMA/06_Annotation2/results/synopsis/topliste_tabdelim/eqtlinfo_2023-07-26_PCSK9_strat.txt")

genes = genes[markername %in% mySNPs,]
eQTLs = eQTLs[snps %in% mySNPs,]
genes = genes[!duplicated(genename),]
eQTLs = eQTLs[!duplicated(genesymbol),]
eQTLs = eQTLs[cistrans == "cis",]

res = data.table(markername = c(eQTLs$snps,genes$markername),
                 genes = c(eQTLs$genesymbol,genes$genename),
                 source = c(rep("eQTL",dim(eQTLs)[1]),rep("proxGene",dim(genes)[1])),
                 chr =c(eQTLs$chr,genes$chr))

res[genes == "MARCH8", genes:="MARCHF8"]
res[genes == "ANUBL1", genes:="ZFAND4"]
res[genes == "FAM21C", genes:="WASHC2C"]
res[genes == "FTHL16", genes:="FTH1P16"]
res[genes == "FTSJD1", genes:="CMTR2"]
res[genes == "MEF2BNB", genes:="BORCS8"]
res[genes == "LST3", genes:="SLCO1B3"]
res[genes == "KIAA1683", genes:="IQCN"]
res[genes == "KIAA0892", genes:="MAU2"]
res[genes == "C2orf43", genes:="LDAH"]
res[genes == "C19orf60", genes:="REX1BD"]
res[genes == "C1orf177", genes:="LEXM"]

res[,dumID := paste(genes,markername, sep="__")]
table(duplicated(res$dumID))
dups = res[duplicated(dumID),dumID]
res[dumID %in% dups,source := "eQTL and proxGene"]
res = res[!duplicated(dumID),]
res

candidateGenes = unique(res$genes)
candidateGenes = candidateGenes[candidateGenes!=""]
myGenTab<-data.table(genename=candidateGenes)

genes38 = fread("../temp/05_HGNC_Download_221125.txt")
table(is.element(myGenTab$genename, genes38$symbol))
myGenTab[!is.element(genename,genes38$symbol),]

myGenTab = myGenTab[is.element(genename,genes38$symbol),]
m1 <- match(myGenTab$gene, genes38$symbol)
genes38 = genes38[m1,]

myGenTab[, `:=`(
  ensg = genes38[, ensembl_gene_id],
  entrez = genes38[, entrez_id],
  hgnc = genes38[, hgnc_id],
  description = genes38[, name],
  type = genes38[,locus_group],
  cytoband = genes38[,location_sortable ]
)]
myGenTab

table(is.na(myGenTab$ensg))
table(is.na(myGenTab$entrez))
setorder(myGenTab,cytoband)

#' # Get GTEx v8 eQTLs ####
#' ***
#' 
myeQTLs<-dir(path = path_GTExv8, pattern = "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8")

gtex8annot <- fread(paste0(path_GTExv8,"GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"))
gtex8annot

dumTab1<-foreach(i=c(1:length(myeQTLs)))%do%{
  #i=1
  myTissue<-gsub(".allpairs.txt.gz","",myeQTLs[i])
  myTissue<-gsub("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_","",myTissue)
  time1 = Sys.time()

  # Step 1: Loading data
  message("Loading eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  data0<-fread(paste0(path_GTExv8, myeQTLs[i]))
  time2 = Sys.time()
  message("          Finished loading eQTL data in ",round(difftime(time2, time1, tz,units = "min"),2)," minutes")

  # Step 2: Change GeneID and filter for candidate genes
  data0[,ENSG:=gsub(gene_id, pattern = "\\..*", replacement = "")]
  # data1 = copy(data0)
  # data0 = copy(data1)
  message("          Filtering eQTL data ",myTissue," number ",i," of ",length(myeQTLs))
  filt<-is.element(data0$ENSG,myGenTab$ensg)
  data0<-data0[filt,]
  time3 = Sys.time()
  message("          Finished filtering eQTL data in ",round(difftime(time3, time2, tz,units = "min"),2)," minutes")

  # Step 3: Get chr, pos and alleles
  message("          Harmonizing column names for ",myTissue," number ",i," of ",length(myeQTLs))
  dummy<-unlist(strsplit(data0$variant_id,"_"))
  chr<-dummy[seq(1,length(dummy),by=5)]
  chr<-gsub("chr","",chr)
  pos<-as.numeric(dummy[seq(2,length(dummy),by=5)])
  other_allele<-dummy[seq(3,length(dummy),by=5)]
  effect_allele<-dummy[seq(4,length(dummy),by=5)]

  data0[,chr:=chr,]
  data0[,pos_b38:=pos,]
  data0[,other_allele:=other_allele,]
  data0[,effect_allele:=effect_allele,]
  data0[,n_samples:=round(ma_count/maf/2,1)]
  setnames(data0,"slope_se","se")
  setnames(data0,"pval_nominal","pval")
  setnames(data0,"slope","beta")

  # Step 4: Get SNP ID of hg19
  data0[,variant_id_b38:= variant_id]
  data0[,variant_id:= NULL]
  matched<-match(data0$variant_id_b38,gtex8annot$variant_id)
  table(is.na(matched))
  table(data0$variant_id_b38==gtex8annot$variant_id[matched])
  data0$variant_id_b37 <- gtex8annot[matched, variant_id_b37]
  table(is.na(data0$variant_id_b37))
  data0<-data0[!is.na(variant_id_b37)]

  dummy<-unlist(strsplit(data0$variant_id_b37,"_"))
  chr_b37<-dummy[seq(1,length(dummy),by=5)]
  filt = chr_b37 == data0$chr
  data0<-data0[filt,]
  dummy<-unlist(strsplit(data0$variant_id_b37,"_"))
  pos<-as.numeric(dummy[seq(2,length(dummy),by=5)])
  data0[,pos_b37:=pos,]
  data0[,chrPos_b37 := paste(chr,pos_b37,sep=":")]
  data0[,chrPos_b37 := paste0("chr",chrPos_b37)]
  data0[,chrPos_b38 := paste(chr,pos_b38,sep=":")]
  data0[,chrPos_b38 := paste0("chr",chrPos_b38)]

  # Step 5: Add gene name and cytoband
  matched<-match(data0$ENSG,myGenTab$ensg)
  data0[,gene:=myGenTab[matched,genename]]
  data0[,cyto:=myGenTab[matched,cytoband]]
  data0[,tissue:=myTissue]
  data0[,SNPGene:=paste(chrPos_b37,gene,sep=":")]

  # Step 6: Check for duplicates (tri allelic positions) and filter for MAF>=1%
  data0[,table(duplicated(SNPGene),duplicated(chrPos_b37))]
  dups = data0[duplicated(data0$SNPGene),]
  data0 = data0[!is.element(SNPGene,dups$SNPGene),]
  data0 = data0[maf>=0.01,]
  x1 = dim(data0)[1]
  message("          Saving ",x1," filtered SNPs")

  # Step 7: Save eqtl data
  if(dir.exists("../temp/05_coloc/")==F) dir.create("../temp/05_coloc/") 
  outfn1<-paste0("../temp/05_coloc/GTEx_v8_filtered_",myTissue,".RData")
  save(data0,file=outfn1)
  # load(outfn1)
  
  # Step 8: Return something
  dummy = data0[,.N,by=gene]
  res = data.table(tissue = myTissue,
                   genes = dummy$gene,
                   n_eQTLs = dummy$N)
  res
}
dumTab1 = rbindlist(dumTab1)
dumTab1
dumTab2<-dcast(dumTab1,
                   formula = genes ~ tissue,
                   value.var = c("n_eQTLs"),
                   sep = "_")
matched = match(myGenTab$genename,dumTab2$genes)
myGenTab = cbind(myGenTab,dumTab2[matched,])
myGenTab[,genes := NULL]

x = dim(myGenTab)[2]
mySums = rowSums(myGenTab[,c(9:x),with=F], na.rm=TRUE)
table(mySums == 0)
myGenTab = myGenTab[mySums != 0,]

save(myGenTab,file="../results/05_1_usedGenes.RData")

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

