#' ---
#' title: "Co-localization Part 2: Run coloc eQTLs"
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
source("../helperFunctions/colocFunction_jp.R")

#' # Prep data ####
#' ***
load("../results/05_1_usedGenes.RData")
load("../results/03_GCTA_COJO.RData")
table(is.element(IndepSignals$candidateGene,myGenTab$genename))
matched = match(IndepSignals$candidateGene,myGenTab$genename)
IndepSignals[,cytoband := myGenTab[matched,cytoband]]

myGenTab[,cytoband2 := cytoband]
myGenTab[cytoband == "10q11.21-q11.22",cytoband2 := "10q11.21"]
myGenTab[cytoband == "10q11.22",cytoband2 := "10q11.21"]
myGenTab[cytoband == "11q12.2-q12.3",cytoband2 := "11q12.2"]
myGenTab[cytoband == "11q12.3",cytoband2 := "11q12.2"]
myGenTab[cytoband == "11q14.1",cytoband2 := "11q12.2"]
myGenTab[cytoband == "12p12.1",cytoband2 := "12p12.2"]
myGenTab[cytoband == "12q24.22-q24.23",cytoband2 := "12q24.22"]
myGenTab[cytoband == "16q22.2-q22.3",cytoband2 := "16q22.2"]
myGenTab[cytoband == "19p12",cytoband2 := "19p13.11"]
myGenTab[cytoband == "19p13",cytoband2 := "19p13.11"]

ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]


#' # Run Coloc ####
#' ***
#' Loop 1: PCSK9 Phenotype (k=1-8)
#' 
#' Loop 2: Tissue (i=1:49 for GTEx tissues)
#' 
#' Loop 3: Genes (all available gene per tissue, independent of top phenotype)
#' 

registerDoMC(cores=20)

dumTab = foreach(k=1:dim(ToDoList)[1])%do%{
  #k=1
  myPheno = ToDoList[k,pheno]
  # myPheno2 = gsub(myPheno,pattern="\\_.*",replacement = "")
  message("Working on phenotype ",myPheno)
  myFileName = ToDoList[k,statistic_path]
  
  data_GWAS = fread(myFileName)
  data_GWAS = data_GWAS[invalidAssocs==F,]
  data_GWAS[,chrPos_b37 := paste0("chr",chr,":",bp_hg19)]
  setnames(data_GWAS,"bp_hg19","pos_b37")
  setnames(data_GWAS,"SE","se")
  setnames(data_GWAS,"nSamples","n_samples")
  setnames(data_GWAS,"EA","effect_allele")
  setnames(data_GWAS,"OA","other_allele")
  data_GWAS[,maf := EAF]
  data_GWAS[EAF>0.5,maf := 1-EAF]
  
  # relevant loci
  # uniqueCytos1 = IndepSignals[pheno == myPheno, unique(cytoband)]
  # uniqueCytos2 = myGenTab[cytoband2 %in% uniqueCytos1, unique(cytoband)]
  uniqueCytos2 = myGenTab[, unique(cytoband)]
  
  myeQTLs2<-dir(path = "../temp/05_coloc/",pattern = ".RData")
  
  dumTab2 = foreach(i=c(1:length(myeQTLs2)))%dopar%{
    #i=4
    loaded = load(paste0("../temp/05_coloc/",myeQTLs2[i]))
    data_eQTLs = get(loaded)
    myTissue2 = gsub("GTEx_v8_filtered_","",myeQTLs2[i])
    myTissue2 = gsub(".RData","",myTissue2)
    message("Working on tissue ",myTissue2)
    data_eQTLs = data_eQTLs[cyto %in% uniqueCytos2,]
    
    # get To Do list
    dummy = data_eQTLs[,.N,by=gene]
    ToDoList2 = copy(myGenTab)
    dummy = dummy[gene %in% ToDoList2[,genename]]
    
    dumTab3 = foreach(j=c(1:dim(dummy)[1]))%do%{
      #j=10
      myRow = dummy[j,]
      moreInfo = copy(ToDoList2)
      moreInfo = moreInfo[genename == myRow$gene,]
      
      data_eQTL1 = copy(data_eQTLs)
      data_eQTL1 = data_eQTL1[gene ==myRow$gene, ]
      dup2 = data_eQTL1[duplicated(chrPos_b37),]
      data_eQTL1 = data_eQTL1[!is.element(chrPos_b37,dup2$chrPos_b37),]
      
      data_GWAS1 = copy(data_GWAS)
      data_GWAS1 = data_GWAS1[chrPos_b37 %in% data_eQTL1$chrPos_b37,]
      dups = data_GWAS1[duplicated(chrPos_b37),]
      data_GWAS1 = data_GWAS1[!is.element(chrPos_b37,dups$chrPos_b37),]
      
      data_eQTL1 = data_eQTL1[chrPos_b37 %in% data_GWAS1$chrPos_b37,]
      
      setorder(data_eQTL1,pos_b37)
      setorder(data_GWAS1,pos_b37)
      stopifnot(data_eQTL1$chrPos_b37 == data_GWAS1$chrPos_b37)
      
      res = colocFunction_jp(tab1 = data_GWAS1,tab2 = data_eQTL1,trait1 = myPheno,trait2 = myTissue2,
                             locus = moreInfo$cytoband,locus_name = myRow$gene,plotting = F,
                             col_SNPID = "chrPos_b37", col_pos = "pos_b37",
                             col_beta = "beta",col_se = "se",col_P = "pval",
                             col_N = "n_samples",col_MAF = "maf",
                             col_effect = "effect_allele",col_other="other_allele")
      x2<-as.data.table(res)
      x3<-t(x2)
      x4<-as.data.table(x3)
      names(x4)<-names(res)
      x4[,cytoband:=moreInfo$cytoband]
      x4[,gene:= myRow$gene]
      x4[,trait1:=myPheno]
      x4[,trait2:=paste0("GE in ",myTissue2)]
      x4
      
    }
    dumTab3 = rbindlist(dumTab3)
    dumTab3
  }
  ColocTable = rbindlist(dumTab2)
  ColocTable
}

ColocTable = rbindlist(dumTab)

ColocTable[,table(PP.H4.abf>=0.75)]
ColocTable[,table(PP.H3.abf>=0.75)]
ColocTable[,table(PP.H2.abf>=0.75)]
ColocTable[,table(PP.H1.abf>=0.75)]
ColocTable[,table(PP.H0.abf>=0.75)]

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H3.abf>=0.75,]

#' # Save results ####
#' ***
ColocTable = ColocTable[,c(7,9,8,10,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Genomic cytoband of index SNP",
                                        "Tested PCSK9 phenotype",
                                        "Tested gene",
                                        "Tested tissue",
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/05_2_coloc_eQTLs.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/05_2_coloc_eQTLs.xlsx"

WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)



#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

