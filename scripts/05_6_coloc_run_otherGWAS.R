#' ---
#' title: "Co-localization Part 3: Run coloc against other GWAS"
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
#' *Co-localization analysis between PCSK9 and other GWAS traits*: 
#' 
#' Used traits:
#' 
#' - lipids
#' - CAD
#' - sleep duration
#' - bilirubin
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath_main,"scripts/"))
source("../helperFunctions/colocFunction_jp.R")

#' # Get ToDoList ####
#' ***
load("../results/03_GCTA_COJO_slct.RData")
load("../results/02_LociOverallPhenotypes_filtered.RData")
IndepSignals = copy(result.5)

#' # Load and filter data ####
#' ***
load("../temp/05_GWAS_allLoci.RData")

#' # Load and filter other GWAS data ####
#' ***
#' ## CAD data ####
myCAD = fread(paste0(path_CAD,"CAD_META.gz"))
setnames(myCAD,"MarkerName","SNP")
setnames(myCAD,"CHR","chr")
setnames(myCAD,"BP","bp_hg19")
setnames(myCAD,"Allele1","EA")
setnames(myCAD,"Allele2","OA")
setnames(myCAD,"Freq1","EAF")
setnames(myCAD,"Effect","beta")
setnames(myCAD,"StdErr","SE")
setnames(myCAD,"P-value","pval")
myCAD[,nSamples := 43159+126268]
myCAD[,MAF := EAF]
myCAD[EAF>0.5,MAF := 1-EAF]
myCAD = myCAD[MAF>=0.01,]
myCAD[,phenotype := "CAD"]
myCAD[,candidateGene := "PCSK9"]
myCAD[,EA:=toupper(EA)]
myCAD[,OA:=toupper(OA)]
myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(myCAD),myNames)
myCAD[,get("colsOut"):=NULL]
setcolorder(myCAD,myNames)

dumTab2 = foreach(j = 1:dim(IndepSignals)[1])%do%{
  # j=1
  myRow = copy(IndepSignals)
  myRow = myRow[j,]
  message("     Working on gene ",myRow$candidateGene)
  
  data_GWAS2 = copy(myCAD)
  data_GWAS2 = data_GWAS2[chr == myRow$chr,]
  data_GWAS2 = data_GWAS2[bp_hg19 >= myRow$region_start,]
  data_GWAS2 = data_GWAS2[bp_hg19 <= myRow$region_end,]
  data_GWAS2 = data_GWAS2[EAF >= 0.01,]
  data_GWAS2[,candidateGene := myRow$candidateGene]
  data_GWAS2
  
}
myCAD = rbindlist(dumTab2)
myCAD[,table(candidateGene)]


#' ## Lipids data ####
myLipids = list.files(path = paste0(path_lipids),pattern = "with_N_1.gz",recursive = T)
myLipids = myLipids[!grepl("tbi",myLipids)]
myLipids

dumTab2 = foreach(i = 1:length(myLipids))%do%{
  #i=1
  myLipid = myLipids[i]
  myLipid = gsub(".*_AFR_EAS_EUR_HIS_SAS_","",myLipid)
  myLipid = gsub("_with_N_1.gz","",myLipid)

  message("Working on lipid ",myLipid)
  
  lipid = fread(paste0(path_lipids,myLipids[i]))
  
  names(lipid)
  lipid = lipid[CHROM %in% IndepSignals$chr,]
  
  dumTab2 = foreach(j = 1:dim(IndepSignals)[1])%do%{
    # j=1
    myRow = copy(IndepSignals)
    myRow = myRow[j,]
    message("     Working on gene ",myRow$candidateGene)
    
    data_GWAS2 = copy(lipid)
    data_GWAS2 = data_GWAS2[CHROM == myRow$chr,]
    data_GWAS2 = data_GWAS2[POS_b37 >= myRow$region_start,]
    data_GWAS2 = data_GWAS2[POS_b37 <= myRow$region_end,]
    data_GWAS2 = data_GWAS2[POOLED_ALT_AF >= 0.01,]
    data_GWAS2[,candidateGene := myRow$candidateGene]
    data_GWAS2
    
  }
  dumTab2 = rbindlist(dumTab2)
  dumTab2[,phenotype := myLipid]
  dumTab2[,METAL_Pvalue := as.numeric(METAL_Pvalue)]
  dumTab2[,pvalue := as.numeric(pvalue)]
  dumTab2[,pvalue_GC := as.numeric(pvalue_GC)]
  dumTab2[,MAF := POOLED_ALT_AF]
  dumTab2[POOLED_ALT_AF>0.5,MAF := 1-POOLED_ALT_AF]
  
  setnames(dumTab2,"CHROM","chr")
  setnames(dumTab2,"POS_b37","bp_hg19")
  setnames(dumTab2,"ALT","EA")
  setnames(dumTab2,"REF","OA")
  setnames(dumTab2,"METAL_Effect","beta")
  setnames(dumTab2,"METAL_StdErr","SE")
  setnames(dumTab2,"METAL_Pvalue","pval")
  setnames(dumTab2,"N","nSamples")
  setnames(dumTab2,"rsID","SNP")
  setnames(dumTab2,"POOLED_ALT_AF","EAF")
  
  myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
  colsOut<-setdiff(colnames(dumTab2),myNames)
  dumTab2[,get("colsOut"):=NULL]
  setcolorder(dumTab2,myNames)
  dumTab2
}

myLipids = rbindlist(dumTab2)
myLipids[,table(candidateGene,phenotype)]

#' ## Sleep duration data ####
mySleep = fread(paste0(path_otherGWAS,"sleepdurationsumstats.txt"))
setnames(mySleep,"CHR","chr")
setnames(mySleep,"BP","bp_hg19")
setnames(mySleep,"ALLELE1","EA")
setnames(mySleep,"ALLELE0","OA")
setnames(mySleep,"A1FREQ","EAF")
setnames(mySleep,"BETA_SLEEPDURATION","beta")
setnames(mySleep,"SE_SLEEPDURATION","SE")
setnames(mySleep,"P_SLEEPDURATION","pval")
mySleep[,nSamples := 446118]
mySleep[,MAF := EAF]
mySleep[EAF>0.5,MAF := 1-EAF]
mySleep = mySleep[MAF>=0.01,]
mySleep[,phenotype := "sleep_duration"]
mySleep[,candidateGene := "PCSK9"]
myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(mySleep),myNames)
mySleep[,get("colsOut"):=NULL]
setcolorder(mySleep,myNames)

dumTab2 = foreach(j = 1:dim(IndepSignals)[1])%do%{
  # j=1
  myRow = copy(IndepSignals)
  myRow = myRow[j,]
  message("     Working on gene ",myRow$candidateGene)
  
  data_GWAS2 = copy(mySleep)
  data_GWAS2 = data_GWAS2[chr == myRow$chr,]
  data_GWAS2 = data_GWAS2[bp_hg19 >= myRow$region_start,]
  data_GWAS2 = data_GWAS2[bp_hg19 <= myRow$region_end,]
  data_GWAS2 = data_GWAS2[EAF >= 0.01,]
  data_GWAS2[,candidateGene := myRow$candidateGene]
  data_GWAS2
  
}
mySleep = rbindlist(dumTab2)
mySleep[,table(candidateGene)]

myBili = fread(paste0(path_otherGWAS,"GCST90019521_buildGRCh37.tsv"))
setnames(myBili,"chromosome","chr")
setnames(myBili,"base_pair_location","bp_hg19")
setnames(myBili,"effect_allele","EA")
setnames(myBili,"other_allele","OA")
setnames(myBili,"standard_error","SE")
setnames(myBili,"p_value","pval")
setnames(myBili,"variant_id","SNP")
myBili[,nSamples := 354368]
myBili[,phenotype := "bilirubin"]
myBili[,EAF := NA]
myBili[,MAF := NA]
myBili[,pval := as.numeric(pval)]
myBili[,EA:=toupper(EA)]
myBili[,OA:=toupper(OA)]
myBili[,candidateGene := "PCSK9"]
myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(myBili),myNames)
myBili[,get("colsOut"):=NULL]
setcolorder(myBili,myNames)

dumTab2 = foreach(j = 1:dim(IndepSignals)[1])%do%{
  # j=1
  myRow = copy(IndepSignals)
  myRow = myRow[j,]
  message("     Working on gene ",myRow$candidateGene)
  
  data_GWAS2 = copy(myBili)
  data_GWAS2 = data_GWAS2[chr == myRow$chr,]
  data_GWAS2 = data_GWAS2[bp_hg19 >= myRow$region_start,]
  data_GWAS2 = data_GWAS2[bp_hg19 <= myRow$region_end,]
  #data_GWAS2 = data_GWAS2[EAF >= 0.01,]
  data_GWAS2[,candidateGene := myRow$candidateGene]
  data_GWAS2
  
}
myBili = rbindlist(dumTab2)
myBili[,table(candidateGene)]

#' ## Save other GWAS data ####
myOtherGWAS = rbind(myCAD,myLipids,mySleep,myBili)
save(myOtherGWAS, file="../temp/05_OtherGWASs.RData")
load("../temp/05_OtherGWASs.RData")

#' # Run Coloc ####
#' ***
data_GWAS[,ChrPos := paste(chr,bp_hg19,sep=":")]
data_GWAS[,dumID := paste(candidateGene,phenotype,sep="__")]
myLoci = unique(data_GWAS$dumID)

dumTab = foreach(i=1:dim(IndepSignals)[1])%do%{
  #i=1
  myRow = IndepSignals[i,]
  otherPhenos = unlist(strsplit(myRow$other_phenotypes, "[|]"))
  otherPhenos = gsub(" ","",otherPhenos)
  phenos = c(myRow$phenotype,otherPhenos)
  
  dumTab2 = foreach(j=1:length(phenos))%do%{
    #j=1
    myRow2 = copy(myRow)
    myRow2[,phenotype := phenos[j]]
    myRow2
  }
  res = rbindlist(dumTab2)
  res
}
IndepSignals = rbindlist(dumTab)

IndepSignals[,dumID := paste(candidateGene,phenotype,sep="__")]
filt = myLoci %in% IndepSignals$dumID
table(filt)
myLoci = myLoci[filt]
data_GWAS = data_GWAS[dumID %in% myLoci,]

myOtherGWAS[,ChrPos := paste(chr,bp_hg19,sep=":")]
myOtherGWAS[,phenotype := gsub("MALE","MALES",phenotype)]

dumTab1 = foreach(i=1:length(myLoci))%do%{
  #i=16
  myLocus = myLoci[i]
  myGene = gsub("__.*","",myLocus)
  myPheno = gsub(".*__","",myLocus)
  
  myRow = copy(IndepSignals)
  myRow = myRow[dumID == myLocus,]
  sex = unlist(strsplit(myPheno,"_"))[2]
  sex = toupper(sex)
  if(sex %in% c("FREE","TREATED"))sex = "ALL"
  
  message("working on phenotype ",myPheno, " and candidate gene ",myGene)
  
  data_GWAS1 = copy(data_GWAS)
  data_GWAS1 = data_GWAS1[dumID == myLocus,]
  data_GWAS1[,chrPosEAOA := paste(chr,bp_hg19,EA,OA,sep="_")]
  data_GWAS1[,chrPosOAEA := paste(chr,bp_hg19,OA,EA,sep="_")]
  
  data_GWAS2 = copy(myOtherGWAS)
  data_GWAS2 = data_GWAS2[candidateGene == myGene,]
  data_GWAS2[,chrPosEAOA := paste(chr,bp_hg19,EA,OA,sep="_")]
  data_GWAS2[,chrPosOAEA := paste(chr,bp_hg19,OA,EA,sep="_")]
  
  otherPhenos = unique(data_GWAS2$phenotype)
  if(sex == "FEMALES") otherPhenos = otherPhenos[c(1,2,4,6,8,10,17,18)]
  if(sex == "MALES") otherPhenos = otherPhenos[c(1,3,5,7,9,11,17,18)]
  if(sex == "ALL") otherPhenos = otherPhenos[c(1,12:16,17,18)]
  data_GWAS2 = data_GWAS2[phenotype %in% otherPhenos,]
  
  dumTab2 = foreach(j=1:length(otherPhenos))%do%{
    #j=8
    myOtherPheno = otherPhenos[j]
    message("    working on other GWAS trait ",myOtherPheno)
    
    data_GWAS3 = copy(data_GWAS2)
    data_GWAS3 = data_GWAS3[phenotype == myOtherPheno,]
    
    data_GWAS4 = copy(data_GWAS1)
    data_GWAS4 = data_GWAS4[chrPosEAOA %in% data_GWAS3$chrPosEAOA | chrPosOAEA %in% data_GWAS3$chrPosEAOA,]
    data_GWAS3 = data_GWAS3[chrPosEAOA %in% data_GWAS4$chrPosEAOA | chrPosOAEA %in% data_GWAS4$chrPosEAOA,]
    data_GWAS3 = data_GWAS3[!duplicated(chrPosEAOA),]
    setorder(data_GWAS3,bp_hg19)
    setorder(data_GWAS4,bp_hg19)
    if(myOtherPheno == "bilirubin")data_GWAS3[,MAF := data_GWAS4$MAF]
    
    res = colocFunction_jp(tab1 = data_GWAS4,tab2 = data_GWAS3,
                           trait1 = myRow$pheno,trait2 = myLipidTrait,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "ChrPos", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myGene]
    x4[,trait1:= myPheno]
    x4[,trait2:=myOtherPheno]
    x4
  }
  
  dumTab2 = rbindlist(dumTab2)
  dumTab2
}

ColocTable = rbindlist(dumTab1)

ColocTable[,table(PP.H4.abf>=0.75)]
ColocTable[,table(PP.H3.abf>=0.75)]
ColocTable[,table(PP.H2.abf>=0.75)]
ColocTable[,table(PP.H1.abf>=0.75)]
ColocTable[,table(PP.H0.abf>=0.75)]

ColocTable[PP.H4.abf>=0.75,]
ColocTable[PP.H3.abf>=0.75,]
ColocTable[PP.H1.abf>=0.75,]

#' # Coloc Plot ####
#' ***
plotData = copy(ColocTable)
plotData[,trait2 := gsub("_INV.*","",trait2)]
plotData[,table(trait2)]
plotData[,dumID := paste(gene,trait1,sep="_")]

data_H4<-dcast(plotData,
               formula = dumID ~ trait2,
               value.var = "PP.H4.abf")
data_H3<-dcast(plotData,
               formula = dumID ~ trait2,
               value.var = "PP.H3.abf")

M4<-as.matrix(data_H4[,-1])
M3<-as.matrix(data_H3[,-1])

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

rownames(M)<-data_H3$dumID
colnames(M)<-names(data_H3)[-1]

# Filter M for top phenos / relevant phenos
rownames(M)
M = M[c(1,2,4,5,9:14,17,18,22,23,27),]

# reorder genes so that it is easier to read
rownames(M)
M = M[c(10:12,2:4,7,5,6,15,8,13,1,14,9),]

# reorder traits so that it is easier to read
colnames(M)
M = M[,c(4,7,3,2,6,5,1,8)]

col_fun = colorRamp2(c(-1,0,1), c("coral","white","steelblue"))

plot1 = Heatmap(M, name = "PP",col=col_fun, 
                cluster_rows = F, cluster_columns = F,
                show_column_dend = FALSE, show_row_dend = FALSE,
                show_heatmap_legend = F,show_column_names = FALSE, 
                bottom_annotation = HeatmapAnnotation(text = anno_text(colnames(M), 
                                                                       rot = 45, 
                                                                       location = unit(1, "npc"), 
                                                                       just = "right"),
                                                      annotation_height = max_text_width(colnames(M))))
plot1
tiff(filename = "../figures/SugFigure_ColocPlot_otherGWAS.tiff", 
     width = 2540, height = 1920, res=300, compression = 'lzw')
plot1
dev.off()

plot2 = Heatmap(M, name = "PP",col=col_fun, 
                cluster_rows = T, cluster_columns = T,
                show_column_dend = FALSE, show_row_dend = FALSE,
                show_heatmap_legend = F,show_column_names = FALSE, 
                bottom_annotation = HeatmapAnnotation(text = anno_text(colnames(M), 
                                                                       rot = 45, 
                                                                       location = unit(1, "npc"), 
                                                                       just = "right"),
                                                      annotation_height = max_text_width(colnames(M))))
plot2
tiff(filename = "../figures/SugFigure_ColocPlot_otherGWAS_clustered.tiff", 
     width = 2540, height = 1920, res=300, compression = 'lzw')
plot2
dev.off()

#' # Save results ####
#' ***
ColocTable = ColocTable[,c(7:9,1:6)]

description = data.table(column = names(ColocTable),
                         description = c("Candidate Gene of tested region", 
                                         "Tested PCSK9 subgroup 1 ",
                                        "Tested PCSK9 subgroup 1 ",
                                        "Number of SNPs included in co-localization analysis per test",
                                        "Posterior probability for hypothesis 0: neither trait associated",
                                        "Posterior probability for hypothesis 1: only trait 1 associated (CKDGen trait)",
                                        "Posterior probability for hypothesis 2: only trait 2 associated (GE trait)",
                                        "Posterior probability for hypothesis 3: both trait associated, but different signals",
                                        "Posterior probability for hypothesis 4: both trait associated, shared signal"))

save(ColocTable, description,file="../results/05_6_coloc_otherGWAS.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/05_6_coloc_otherGWAS.xlsx"

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

