#' ---
#' title: "GCTA COJO"
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
#' For each loci, I want to check for multiple associations using GCTA - conditional joint (COJO) analyses. 
#' 
#' * Step 1: get GCTA input (.ma files) of each phenotype
#' * Step 2: perform GCTA COJO select to identify independent signals per loci
#' * Step 3: in case of multiple independent signals at a locus, perform GCTA COJO conditional to estimate the conditional statistics at this locus for each independent signal. 
#' 
#' Reference data set is the combination of LIFE-Adult and LIFE-Heart.
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

#' # Loop 1: get GCTA file format ####
#' ***
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  
  # load data
  data = fread(ToDoList[i,statistic_path])
  x1 = dim(data)[1]
  
  # check & filter SNPs
  data<-data[invalidAssocs==F,]
  stopifnot(min(data$EAF)>=0.01)
  stopifnot(max(data$EAF)<=0.99)
  stopifnot(min(data$info)>=0.8)
  stopifnot(max(data$I2)<=0.9)
  stopifnot(min(data$nStudies) >=2)
  stopifnot(sum(is.na(data$beta_score))==0)
  x2 = dim(data)[1]
  
  # filter columns
  myNames1<-c("markername","EA","OA","EAF","beta","SE","pval","nSamples")
  stopifnot(sum(is.element(colnames(data),myNames1))==8)
  colsOut<-setdiff(colnames(data),myNames1)
  data[,get("colsOut"):=NULL]
  setcolorder(data,myNames1)
  myNames2<-c("SNP","A1","A2","freq","b","se","p","N")
  names(data)<-myNames2
  
  # save as .ma file
  if(dir.exists("../temp/05_GCTA_input/")==F) dir.create("../temp/05_GCTA_input/") 
  outfn = paste0("../temp/05_GCTA_input/",ToDoList[i, pheno],".ma")
  fwrite(data,file=outfn,sep = "\t")
  
  # return some numbers
  res = data.table(pheno = ToDoList[i, pheno],
                   nSNPs_unfiltered = x1,
                   nSNPs_filtered = x2)
  res
  
}
dumTab = rbindlist(dumTab)
dumTab

#' # Loop 2: perform GCTA COJO select ####
#' ***
#' For the *PCSK9* locus on chromosome 1, I want to test each trait for independent signal (main locus). For the other loci, I only want to test those with 2 or more supportive variants with p<1e-6 and only the best setting. 
#' 
load("../results/02_LociOverallPhenotypes_filtered.RData")
load("../results/02_LociPerPhenotype.RData")

result.3 = result.3[markername %in% result.5$markername,]
matched = match(result.3$markername,result.5$markername)
result.3[,candidateGene := result.5[matched,candidateGene]]
result.3[,region := result.5[matched,region]]
result.3[,input := paste0("../temp/05_GCTA_input/",phenotype,".ma")]
result.3[,cutoff := 5e-8]
result.3[pval>cutoff,cutoff := 1e-6]
result.3[candidateGene == "KHDRBS2",cutoff := 1e-6]
result.3[candidateGene == "KHDRBS2",markername := "rs112875382:62600689:G:GT"]
result.3[candidateGene == "SLCO1B3",markername := "rs4762806:21067768:T:C"]
result.3 = result.3[-17,]

SLCTtab<-foreach(i = 1:dim(result.3)[1])%do%{
  #i=1
  myRow = result.3[i,]
  myPheno<-myRow[1,phenotype]
  mySNP<-myRow[1,markername]
  myChr<-myRow[1,chr]
  myNR<-myRow[1,region]
  myInput<-myRow[1,input]
  myCut<-myRow[1,cutoff]
  myOutput = paste(i,myRow[1,candidateGene],myRow[1,phenotype],sep="_")
  if(dir.exists("../results/05_GCTA_COJO_slct/")==F) dir.create("../results/05_GCTA_COJO_slct/") 
  
  mycall<-paste0(path_gcta,
                 " --bfile ", path_gcta_RefData,myChr,
                 " --chr ",myChr,
                 " --maf 0.01",
                 " --cojo-file ",myInput,
                 " --cojo-p ", myCut,
                 " --cojo-slct",
                 " --extract-region-snp ",mySNP," 500",
                 " --out ../results/05_GCTA_COJO_slct/",myOutput)
  mycall
  system(mycall)
  
  res<-fread(paste0("../results/05_GCTA_COJO_slct/",myOutput,".jma.cojo"))
  res[,pheno:=myPheno]
  res[,TopSNP:=mySNP]
  res[,num:=i]
  res[,input_slct:= paste0("../results/05_GCTA_COJO_slct/",myOutput,".jma.cojo")]
  res[,region_start := myRow$region_start]
  res[,region_end := myRow$region_end]
  res[,candidateGene := myRow$candidateGene]
  res[,EA := myRow$EA]
  res[,OA := myRow$OA]
  res[,EAF := myRow$EAF]
  res[,info := myRow$info]
  res
  
}
SLCTtab<-rbindlist(SLCTtab)
table(SLCTtab$SNP == SLCTtab$TopSNP)
SLCTtab[grepl("PCSK9",candidateGene)]
SLCTtab[!grepl("PCSK9",candidateGene)]
SLCTtab[,rsID := gsub(":.*","",SNP)]


#' **Summary of the PCSK9 locus**
#' 
#' Multiple hits!
#' 
#' **Summary of the other loci**
#' 
#' All other loci have only one signal!
#' 
SLCTtab[, multipleSignals := F]
SLCTtab[num<6, multipleSignals := T]

SLCTtab[multipleSignals == F, input_CS := paste0("../data/SumStat_",pheno,"_230120.txt.gz")]


#' # Loop 2: perform GCTA COJO conditional #### 
#' ***
#' I only do this for the PCSK9 locus, which has multiple signals for four traits
#' 
SLCTtab2 = copy(SLCTtab)
SLCTtab2 = SLCTtab2[multipleSignals == T,]
myPhenos = unique(SLCTtab2$pheno)

for(i in 1:length(myPhenos)){
  #i=1
  myPheno<-myPhenos[i]
  res1 = copy(SLCTtab2)
  res1<-res1[pheno == myPheno,]
  
  message("Working on trait ",myPheno,", ",i," of ",length(myPhenos),", with n=",dim(res1)[1]," independent signals ...")
  
  myChr<-res1$Chr[1]
  
  for(j in 1:dim(res1)[1]){
    #j=1
    
    # write SNP to exclude
    SNPList<-res1$SNP[-j]
    mySNP<-res1$SNP[j]
    myRSID = res1$rsID[j]
    myPheno<-res1$pheno[j]
    myInput<-result.5[phenotype==myPheno,input]
    myInput = myInput[1]
    
    if(dir.exists("../results/05_GCTA_COJO_cond/")==F) dir.create("../results/05_GCTA_COJO_cond/") 
    write.table(SNPList,file=paste0("../results/05_GCTA_COJO_cond/PCSK9::",myPheno,"_signal_",myRSID,".snplist"),
                col.names = F,row.names = F,quote = F)
    
    # run GCTA cond
    mycall<-paste0(path_gcta,
                   " --bfile ", path_gcta_RefData,myChr,
                   " --chr ",myChr,
                   " --maf 0.01",
                   " --cojo-file ",myInput,
                   " --cojo-cond ../results/05_GCTA_COJO_cond/PCSK9::",myPheno,"_signal_",myRSID,".snplist",
                   " --extract-region-snp ",mySNP," 500",
                   " --out ../results/05_GCTA_COJO_cond/PCSK9::",myPheno,"_signal_",myRSID)
    
    mycall
    system(mycall)
    
    # add input
    SLCTtab[pheno == myPheno & rsID == myRSID,
            input_CS := paste0("../results/05_GCTA_COJO_cond/PCSK9::",myPheno,"_signal_",myRSID,".cma.cojo")]
  }
  
}


#' # Save results #### 
#' ***
#' I want to save the select table

IndepSignals = copy(SLCTtab)
save(IndepSignals,file="../results/05_GCTA_COJO.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

