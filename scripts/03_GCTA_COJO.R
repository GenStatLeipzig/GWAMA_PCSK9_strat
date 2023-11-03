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
#' For *PCSK9*, I want to check for multiple associations using GCTA - conditional joint (COJO) analyses.
#' 
#' - Step 1: get GCTA input (.ma files) of each subgroup
#' - Step 2: perform GCTA COJO select to identify independent signals per subgroup
#' - Step 3: test pairwise LD of all selected SNPs (done in online tool)
#' - Step 4: perform GCTA COJO joint to estimate the joint effects of the independent SNPs in each subgroup
#' - Step 5: perform GCTA COJO conditional to estimate the conditional statistics at this locus for each independent signal. 
#' 
#' Reference data set is the combination of LIFE-Adult and LIFE-Heart.
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

#' # Step 1: get GCTA file format ####
#' ***
load("../results/02_LociOverallPhenotypes.RData")
result.4 = result.4[markername == "rs11591147:55505647:G:T",]

ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]

if(dir.exists("../temp/03_GCTA_input/")==F) dir.create("../temp/03_GCTA_input/") 

dumTab1 = foreach(i = 1:dim(ToDoList)[1])%do%{
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
  data2 = copy(data)
  
  # filter columns
  myNames1<-c("markername","EA","OA","EAF","beta","SE","pval","nSamples")
  stopifnot(sum(is.element(colnames(data),myNames1))==8)
  colsOut<-setdiff(colnames(data),myNames1)
  data[,get("colsOut"):=NULL]
  setcolorder(data,myNames1)
  myNames2<-c("SNP","A1","A2","freq","b","se","p","N")
  names(data)<-myNames2
  
  # save as .ma file
  outfn = paste0("../temp/03_GCTA_input/",ToDoList[i, pheno],".ma")
  fwrite(data,file=outfn,sep = "\t")
  
  # return PCSK9 locus
  data2 = data2[chr == result.4$chr,]
  data2 = data2[bp_hg19 >= result.4$region_start,]
  data2 = data2[bp_hg19 <= result.4$region_end,]
  data2
  
}
INPUTtab = rbindlist(dumTab1)
INPUTtab

SNPList = unique(INPUTtab$markername)
fn0 = "../temp/03_GCTA_input/PCSK9.snplist"
write.table(unique(INPUTtab$markername),
            file=fn0,col.names = F,row.names = F,quote = F)

#' # Step 2: perform GCTA COJO select ####
#' ***
load("../results/02_LociPerPhenotype.RData")
result.3 = result.3[markername == "rs11591147:55505647:G:T",]
result.3[,input := paste0("../temp/03_GCTA_input/",phenotype,".ma")]
if(dir.exists("../results/03_GCTA_COJO_slct/")==F) dir.create("../results/03_GCTA_COJO_slct/") 

dumTab2<-foreach(i = 1:dim(result.3)[1])%do%{
  #i=1
  myRow = result.3[i,]
  myPheno<-myRow[1,phenotype]
  mySNP<-myRow[1,markername]
  myChr<-myRow[1,chr]
  myNR<-myRow[1,region]
  myInput<-myRow[1,input]
  myCut<-5e-8
  myOutput = paste(i,myRow[1,phenotype],sep="_")
  
  mycall<-paste0(path_gcta,
                 " --bfile ", path_gcta_RefData,myChr,
                 " --chr ",myChr,
                 " --maf 0.01",
                 " --cojo-file ",myInput,
                 " --cojo-p ", myCut,
                 " --cojo-slct",
                 " --extract ", fn0,
                 " --out ../results/03_GCTA_COJO_slct/",myOutput)
  mycall
  system(mycall)
  
  res<-fread(paste0("../results/03_GCTA_COJO_slct/",myOutput,".jma.cojo"))
  res[,pheno:=myPheno]
  res[,TopSNP:=mySNP]
  res[,num:=i]
  res[,input_slct:= paste0("../results/03_GCTA_COJO_slct/",myOutput,".jma.cojo")]
  res[,region_start := myRow$region_start]
  res[,region_end := myRow$region_end]
  
  # match to input to get SNP infos
  data2 = copy(INPUTtab)
  data2 = data2[phenotype == myPheno]
  data2 = data2[markername %in% res$SNP,]
  matched = match(res$SNP,data2$markername)
  data2 = data2[matched,]
  res[,EA := data2$EA]
  res[,OA := data2$OA]
  res[,EAF := data2$EAF]
  res[,info := data2$info]
  res[,I2 := data2$I2]
  res
  
}
SLCTtab<-rbindlist(dumTab2)
table(SLCTtab$SNP == SLCTtab$TopSNP)
SLCTtab[,rsID := gsub(":.*","",SNP)]

save(SLCTtab, file = "../results/03_GCTA_COJO_slct.RData")

#' # Step 3: LD look-up ####
#' ***
PCSK9_SNPs = SLCTtab[,unique(rsID)]
PCSK9_SNPs

#' Check SNPs in [LDlink](https://ldlink.nih.gov/?tab=ldmatrix) (using EUR population)
#' 
#' - Cluster 1: rs2495491
#' - Cluster 2: rs11591147
#' - Cluster 3: rs11583680 and rs28385704
#' - Cluster 4: rs2495477, rs472495, and rs693668
#' 
SLCTtab[rsID %in% PCSK9_SNPs[c(5,2)],c(15,26,5:8,11:13)]

#' strongest signal in males for rs11583680 --> keep this SNP as proxy for all
#' 
SLCTtab[rsID %in% PCSK9_SNPs[c(6,3,7)],c(15,26,5:8,11:13)]

#' strongest signal in males for rs693668 --> keep this SNP as proxy for all
#' 
#' However, this SNP is not available in statin-free individuals (filtered due to I2>90%). Hence, for this subgroup I will generate joint and conditional statistics using rs2495477 instead.
#' 
#' In the interaction test & MR, I want a shared SNP set. There I will use rs693668 for all subgroups. 
#' 
#' # Step 4: perform GCTA COJO joint ####
#' ***
if(dir.exists("../results/03_GCTA_COJO_joint/")==F) dir.create("../results/03_GCTA_COJO_joint/") 

fn1 = paste0("../results/03_GCTA_COJO_joint/PCSK9.snplist")
write.table(unique(SLCTtab[rsID %in% PCSK9_SNPs[c(1,4,5,7)],SNP]),
            file=fn1,
            col.names = F,row.names = F,quote = F)
fn2 = paste0("../results/03_GCTA_COJO_joint/PCSK9_free.snplist")
write.table(unique(SLCTtab[rsID %in% PCSK9_SNPs[c(1,4,5,6)],SNP]),
            file=fn2,
            col.names = F,row.names = F,quote = F)

dumTab3 = foreach(i = 1:dim(result.3)[1])%do%{
  #i=1
  myPheno<-result.3[i,phenotype]
  myChr<-1
  myInput<-result.3[i,unique(input)]
  
  fn = fn1
  if(i==4) fn = fn2
  
  message("Working on trait ",myPheno,", ",i," of ",dim(result.3)[1],", with n=4 independent signals ...")
  
  # run GCTA joint
  mycall<-paste0(path_gcta,
                 " --bfile ", path_gcta_RefData,myChr,
                 " --chr ",myChr,
                 " --maf 0.01",
                 " --cojo-file ",myInput,
                 " --extract ", fn,
                 " --cojo-joint ",
                 " --out ../results/03_GCTA_COJO_joint/",myPheno)
  
  mycall
  system(mycall)
  
  #load result
  res<-fread(paste0("../results/03_GCTA_COJO_joint/",myPheno,".jma.cojo"))
  res[,pheno:=myPheno]
  res
}
JOINTtab = rbindlist(dumTab3)  
JOINTtab[,Tstat := b/se]
JOINTtab[,Fstat := (b/se)^2]
JOINTtab[,TstatJ := bJ/bJ_se]
JOINTtab[,FstatJ := (bJ/bJ_se)^2]

JOINTtab[Fstat <10,]
JOINTtab[FstatJ<10,]

save(JOINTtab, file = "../results/03_GCTA_COJO_joint.RData")

#' # Step 5: perform GCTA COJO conditional #### 
#' ***

if(dir.exists("../results/03_GCTA_COJO_cond/")==F) dir.create("../results/03_GCTA_COJO_cond/") 
dummy = c(1,4,5,7,6)

for(i in 1:5){
  #i=3
  j = dummy[i]
  fn = paste0("../results/03_GCTA_COJO_cond/PCSK9_",PCSK9_SNPs[j],".snplist")
  if(i!=5){
    PCSK9_SNPs_filtered = unique(SLCTtab[rsID %in% PCSK9_SNPs[c(1,4,5,7)],SNP])
    PCSK9_SNPs_filtered = PCSK9_SNPs_filtered[-i]
  }else{
    PCSK9_SNPs_filtered = unique(SLCTtab[rsID %in% PCSK9_SNPs[c(1,4,5)],SNP])
  }
  write.table(PCSK9_SNPs_filtered,
              file=fn,
              col.names = F,row.names = F,quote = F)
}

for(i in 1:dim(result.3)[1]){
  #i=1
  myPheno<-result.3[i, phenotype]
  myChr<-1
  myInput<-result.3[i,input]
  
  message("Working on trait ",myPheno,", ",i," of ",dim(result.3)[1],", with n=4 independent signals ...")
  
  PCSK9_indepSNPs_rsID = PCSK9_SNPs[c(1,4,5,7)]
  if(myPheno == "PCSK9_free") PCSK9_indepSNPs_rsID = PCSK9_SNPs[c(1,4,5,6)]

  for(j in 1:4){
    #j=1
    fn = paste0("../results/03_GCTA_COJO_cond/PCSK9_",PCSK9_indepSNPs_rsID[j],".snplist")
    mySNP = "rs11591147:55505647:G:T"
    myRSID = PCSK9_indepSNPs_rsID[j]
    
    # run GCTA cond
    mycall<-paste0(path_gcta,
                   " --bfile ", path_gcta_RefData,myChr,
                   " --chr ",myChr,
                   " --maf 0.01",
                   " --cojo-file ",myInput,
                   " --cojo-cond ",fn,
                   " --extract ", fn0,
                   " --out ../results/03_GCTA_COJO_cond/",myPheno,"_signal_",myRSID)
    
    mycall
    system(mycall)
    
  }

}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

