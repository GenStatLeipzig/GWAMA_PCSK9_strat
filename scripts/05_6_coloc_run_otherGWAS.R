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
#' + Lipids (sex-stratified input from GLGC for HDL, LDL, TC, TG, nonHDL) --> all loci!
#' + Total bilirubin levels (Sinnott-Armstrong N, GWAS catalog study accession GCST90019521) --> *SLCO1B1*
#' + Sleep duration (Dashti HS, GWAS catalog study accession GCST007561) --> *NOS1*
#' + Systolic blood pressure, Pulse pressure, Medication use --> *PLB1* 
#' 
#' For locus *PCSK9*, I will test all 6 strata with all 5 lipid traits
#'  
#' For the other loci, I will test the best associated phenotype only. 
#' 
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
load("../results/05_GCTA_COJO.RData")

ToDoList = copy(IndepSignals)
ToDoList = ToDoList[,c(1:5,15,19:21,26)]
ToDoList = ToDoList[c(1,3,5,7,9,11,12,13:24)]
ToDoList[,otherGWAS := "lipids"]

ToDoList2 = copy(ToDoList)
ToDoList2 = ToDoList2[candidateGene == "PCSK9",]
ToDoList2[,otherGWAS := "CAD"]

ToDoList3 = copy(ToDoList)
ToDoList3 = ToDoList3[candidateGene == "SLCO1B3",]
ToDoList3[,otherGWAS := "bilirubin"]

ToDoList4 = copy(ToDoList)
ToDoList4 = ToDoList4[candidateGene == "NOS1",]
ToDoList4[,otherGWAS := "sleep"]

ToDoList5 = copy(ToDoList)
ToDoList5 = ToDoList5[candidateGene == "FOSL2",]
ToDoList5 = ToDoList5[c(1,1,1),]
ToDoList5[,otherGWAS := c("SBP","PP","medication")]


ToDoList6 = copy(ToDoList)
ToDoList6 = ToDoList6[candidateGene == "KHDRBS2",]
ToDoList6[,otherGWAS := c("chronotype")]

ToDoList = rbind(ToDoList,ToDoList2,ToDoList3,ToDoList4,ToDoList5,ToDoList6)
ToDoList

#' # Load and filter data ####
#' ***
load("../temp/06_GWAS_allLoci.RData")

#' # Load and filter other GWAS data ####
#' ***
#' ## CAD data ####
myCAD = fread(paste0(path_CAD,"CAD_META.gz"))
myCAD = myCAD[CHR == 1]
myCAD = myCAD[BP <= ToDoList[1,region_end]]
myCAD = myCAD[BP >= ToDoList[1,region_start]]
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
myCAD[,candidateGene := ToDoList[1,candidateGene]]
myCAD[,EA:=toupper(EA)]
myCAD[,OA:=toupper(OA)]
myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(myCAD),myNames)
myCAD[,get("colsOut"):=NULL]
setcolorder(myCAD,myNames)

#' ## Lipids data ####
ToDoList2 = copy(ToDoList)
ToDoList2 = ToDoList2[otherGWAS=="lipids",]
ToDoList2 = ToDoList2[!duplicated(candidateGene),]

myLipids = list.files(path = paste0(path_lipids),pattern = "with_N_1.gz",recursive = T)
myLipids = myLipids[!grepl("tbi",myLipids)]
myLipids

dumTab2 = foreach(i = 1:length(myLipids))%do%{
  #i=15
  myLipid = myLipids[i]
  myLipid = gsub(".*_AFR_EAS_EUR_HIS_SAS_","",myLipid)
  myLipid = gsub("_with_N_1.gz","",myLipid)

  message("Working on lipid ",myLipid)
  
  lipid = fread(paste0(path_lipids,myLipids[i]))
  
  names(lipid)
  lipid = lipid[CHROM %in% ToDoList2$Chr,]
  
  dumTab2 = foreach(j = 1:dim(ToDoList2)[1])%do%{
    # j=1
    myRow = copy(ToDoList2)
    myRow = myRow[j,]
    message("     Working on gene ",myRow$candidateGene)
    
    data_GWAS2 = copy(lipid)
    data_GWAS2 = data_GWAS2[CHROM == myRow$Chr,]
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

#' ## Sleep & bilirubin data ####
ToDoList3 = copy(ToDoList)
ToDoList3 = ToDoList3[otherGWAS!="lipids",]

mySleep = fread(paste0(path_otherGWAS,"sleepdurationsumstats.txt"))
mySleep = mySleep[CHR == 12]
mySleep = mySleep[BP <= ToDoList3[otherGWAS=="sleep",region_end]]
mySleep = mySleep[BP >= ToDoList3[otherGWAS=="sleep",region_start]]
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
mySleep[,phenotype := ToDoList3[otherGWAS=="sleep",otherGWAS]]
mySleep[,candidateGene := ToDoList3[otherGWAS=="sleep",candidateGene]]

myBili = fread(paste0(path_otherGWAS,"GCST90019521_buildGRCh37.tsv"))
myBili = myBili[chromosome == 12]
myBili = myBili[base_pair_location <= ToDoList3[otherGWAS=="bilirubin",region_end]]
myBili = myBili[base_pair_location >= ToDoList3[otherGWAS=="bilirubin",region_start]]
setnames(myBili,"chromosome","chr")
setnames(myBili,"base_pair_location","bp_hg19")
setnames(myBili,"effect_allele","EA")
setnames(myBili,"other_allele","OA")
setnames(myBili,"standard_error","SE")
setnames(myBili,"p_value","pval")
setnames(myBili,"variant_id","SNP")
myBili[,nSamples := 354368]
myBili[,phenotype := ToDoList3[otherGWAS=="bilirubin",otherGWAS]]
myBili[,EAF := NA]
myBili[,MAF := NA]
myBili[,pval := as.numeric(pval)]
myBili[,EA:=toupper(EA)]
myBili[,OA:=toupper(OA)]
myBili[,candidateGene := ToDoList3[otherGWAS=="bilirubin",candidateGene]]

myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(mySleep),myNames)
mySleep[,get("colsOut"):=NULL]
setcolorder(mySleep,myNames)

colsOut<-setdiff(colnames(myBili),myNames)
myBili[,get("colsOut"):=NULL]
setcolorder(myBili,myNames)

#' ## Hypertension data ####
mySBP = fread(paste0(path_otherGWAS,"34594039-GCST90018972-EFO_0006335-Build37.f.tsv.gz"))
myPP = fread(paste0(path_otherGWAS,"34594039-GCST90018970-EFO_0005763-Build37.f.tsv.gz"))
myMed = fread(paste0(path_otherGWAS,"34594039-GCST90018988-EFO_0009931-Build37.f.tsv.gz"))

mySBP = mySBP[chromosome == ToDoList3[otherGWAS=="SBP",Chr]]
mySBP = mySBP[base_pair_location <= ToDoList3[otherGWAS=="SBP",region_end]]
mySBP = mySBP[base_pair_location >= ToDoList3[otherGWAS=="SBP",region_start]]
mySBP[,phenotype := ToDoList3[otherGWAS=="SBP",otherGWAS]]
mySBP[,nSamples := 340159 ]

myPP = myPP[chromosome == ToDoList3[otherGWAS=="PP",Chr]]
myPP = myPP[base_pair_location <= ToDoList3[otherGWAS=="PP",region_end]]
myPP = myPP[base_pair_location >= ToDoList3[otherGWAS=="PP",region_start]]
myPP[,phenotype := ToDoList3[otherGWAS=="PP",otherGWAS]]
myPP[,nSamples := 360863  ]

myMed = myMed[chromosome == ToDoList3[otherGWAS=="medication",Chr]]
myMed = myMed[base_pair_location <= ToDoList3[otherGWAS=="medication",region_end]]
myMed = myMed[base_pair_location >= ToDoList3[otherGWAS=="medication",region_start]]
myMed[,phenotype := ToDoList3[otherGWAS=="medication",otherGWAS]]
myMed[,nSamples := 237530 ]

myBP = rbind(mySBP,myPP,myMed)
setnames(myBP,"chromosome","chr")
setnames(myBP,"base_pair_location","bp_hg19")
setnames(myBP,"effect_allele","EA")
setnames(myBP,"other_allele","OA")
setnames(myBP,"standard_error","SE")
setnames(myBP,"p_value","pval")
setnames(myBP,"variant_id","SNP")
setnames(myBP,"effect_allele_frequency","EAF")
myBP[,MAF := EAF]
myBP[EAF>0.5,MAF := 1-EAF]
myBP[,pval := as.numeric(pval)]
myBP[,SNP := paste(chr,bp_hg19,sep=":")]
myBP[,candidateGene := "FOSL2"]

myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(myBP),myNames)
myBP[,get("colsOut"):=NULL]
setcolorder(myBP,myNames)

#' ## Chronotype data ####
myChrono = fread(paste0(path_otherGWAS,"morning_person_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3_logORs.txt.gz"))
myChrono = myChrono[CHR == ToDoList[otherGWAS=="chronotype",Chr]]
myChrono = myChrono[BP <= ToDoList[otherGWAS=="chronotype",region_end]]
myChrono = myChrono[BP >= ToDoList[otherGWAS=="chronotype",region_start]]
myChrono[,phenotype := ToDoList[otherGWAS=="chronotype",otherGWAS]]
myChrono[,candidateGene := ToDoList[otherGWAS=="chronotype",candidateGene]]
myChrono[,nSamples := 403195 ]
setnames(myChrono,"CHR","chr")
setnames(myChrono,"BP","bp_hg19")
setnames(myChrono,"ALLELE1","EA")
setnames(myChrono,"ALLELE0","OA")
setnames(myChrono,"A1FREQ","EAF")
setnames(myChrono,"LOGOR","beta")
setnames(myChrono,"LOGOR_SE","SE")
setnames(myChrono,"P_BOLT_LMM","pval")
myChrono[,MAF := EAF ]
myChrono[EAF>0.5,MAF := 1-EAF ]
myNames = c("SNP","chr","bp_hg19","EA","OA","EAF","MAF","beta","SE","pval","nSamples","candidateGene","phenotype")
colsOut<-setdiff(colnames(myChrono),myNames)
myChrono[,get("colsOut"):=NULL]
setcolorder(myChrono,myNames)

#' ## Save other GWAS data ####
myOtherGWAS = rbind(myCAD,myLipids,mySleep,myBili,myBP,myChrono)
save(myOtherGWAS, file="../temp/06_OtherGWASs.RData")
load("../temp/06_OtherGWASs.RData")

#' # Run Coloc ####
#' ***
data_GWAS[,ChrPos := paste(chr,bp_hg19,sep=":")]
myOtherGWAS[,ChrPos := paste(chr,bp_hg19,sep=":")]
myOtherGWAS[,phenotype := gsub("MALE","MALES",phenotype)]

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=32
  myRow = copy(ToDoList)
  myRow = myRow[i,]
  sex = unlist(strsplit(myRow$pheno,"_"))[2]
  sex = toupper(sex)
  if(sex %in% c("FREE","TREATED"))sex = "ALL"
  
  message("working on phenotype ",myRow$pheno, " and candidate gene ",myRow$candidateGene)
  
  data_GWAS1 = copy(data_GWAS)
  data_GWAS1 = data_GWAS1[candidateGene == myRow$candidateGene,]
  data_GWAS1 = data_GWAS1[phenotype == myRow$pheno ,]
  
  data_GWAS2 = copy(myOtherGWAS)
  data_GWAS2 = data_GWAS2[candidateGene == myRow$candidateGene,]
  
  if(myRow$otherGWAS == "lipids"){
    myLipidTraits = unique(data_GWAS2$phenotype)
    myLipidTraits = myLipidTraits[grepl(paste0("_",sex),myLipidTraits)]
    
    dumTab2 = foreach(j=1:length(myLipidTraits))%do%{
      #j=1
      myLipidTrait = myLipidTraits[j]
      message("    working on lipid ",myLipidTrait)
      
      data_GWAS3 = copy(data_GWAS2)
      data_GWAS3 = data_GWAS3[phenotype == myLipidTrait,]
      
      data_GWAS4 = copy(data_GWAS1)
      data_GWAS4 = data_GWAS4[ChrPos %in% data_GWAS3$ChrPos,]
      data_GWAS3 = data_GWAS3[ChrPos %in% data_GWAS4$ChrPos,]
      
      setorder(data_GWAS3,bp_hg19)
      setorder(data_GWAS4,bp_hg19)
      
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
      x4[,gene:= myRow$candidateGene]
      x4[,trait1:= myRow$pheno]
      x4[,trait2:=myLipidTrait]
      x4
    }
    dumTab2 = rbindlist(dumTab2)
  }else{
    message("    working on trait ",myRow$otherGWAS)
    
    data_GWAS3 = copy(data_GWAS2)
    data_GWAS3 = data_GWAS3[phenotype == myRow$otherGWAS]
    trait2 = unique(data_GWAS3$phenotype)
    
    data_GWAS4 = copy(data_GWAS1)
    data_GWAS4 = data_GWAS4[ChrPos %in% data_GWAS3$ChrPos,]
    data_GWAS3 = data_GWAS3[ChrPos %in% data_GWAS4$ChrPos,]
    
    setorder(data_GWAS3,bp_hg19)
    setorder(data_GWAS4,bp_hg19)
    
    if(sum(!is.na(data_GWAS3$MAF))==0){
      data_GWAS3[,MAF := NULL]
    }

    res = colocFunction_jp(tab1 = data_GWAS4,tab2 = data_GWAS3,
                           trait1 = myRow$pheno,trait2 = trait2,
                           locus = myRow$rsID,locus_name = myRow$candidateGene,plotting = F,
                           col_SNPID = "ChrPos", col_pos = "bp_hg19",
                           col_beta = "beta",col_se = "SE",col_P = "pval",
                           col_N = "nSamples",col_MAF = "MAF",
                           col_effect = "EA",col_other="OA")
    x2<-as.data.table(res)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(res)
    x4[,gene:= myRow$candidateGene]
    x4[,trait1:= myRow$pheno]
    x4[,trait2:=trait2]
    dumTab2 = x4
  }
  dumTab2
  ColocTable = copy(dumTab2)
  ColocTable
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

#' **Summary**: 
#' 
#' * PCSK9 males: shared with HDL! 
#' * NOS1 males: shared with sleep duration!
#' * PLB1, CYP51A1P3, & SLCO1B1 have only signal in females (not in males)
#' * MYNC has signal only in males & only in free setting
#' * PRKAG2 & NOS1 have only signal in males (not in females)
#' * TM6SF2 has signal only in free setting
#' 
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

save(ColocTable, description,file="../results/06_6_coloc_otherGWAS.RData")

tosave4 = data.table(data = c("ColocTable", "description"), 
                     SheetNames = c("ColocTable", "Description"))
excel_fn = "../results/06_6_coloc_otherGWAS.xlsx"

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

