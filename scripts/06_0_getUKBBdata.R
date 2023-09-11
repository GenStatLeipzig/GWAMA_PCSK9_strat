#' ---
#' title: "Prep Cis-MR: Getting LDLC summary statistics from UKBB"
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
#' I want to extract from the UKBB phenotype file all necessary columns for my cis-MR using PCSK9 (stratified) as exposure and LDLC as outcome. In addition to LDLC, I want to test for an direct (not mediated by LDLC) effect on CAD. 
#' 
#' I want the following columns
#' 
#' - ID
#' - Sex
#' - Age
#' - Medication
#' - LDLC levels
#' - Genetic principal components (the first 10)
#' - Kinship
#' - Ethnic background
#' - Diagnoses: acute myocardial infarction
#' - Diagnoses: subsequent myocardial infarction
#' - Diagnoses: certain current complications following acute MI
#' - Diagnoses: other acture ischaemic heart diseases
#' - Diagnoses: chronic ischemic heart diseases
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

.libPaths()
library(data.table)
setDTthreads(2)
library(readxl)
library(foreach)
library(meta)

#' # Load data ####
#' ***
myTab_20 = fread(paste0(path_UKBB,"phenotypes/ukb672224.tab"), 
                 header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

#' # Get all columns I want ####
#' ***
myAnnot[colNm == "f.eid", parameter := "ID"]
myAnnot[colNm == "f.eid", comment := "pseudoID"]

myAnnot[colNm == "f.31.0.0", parameter := "Sex"]
myAnnot[colNm == "f.31.0.0", comment := "coding females as 0 and males as 1"]

myAnnot[colNm %in% c("f.21022.0.0"), parameter := "Age at recuitment (years)"]
myAnnot[colNm %in% c("f.21022.0.0"), comment := c("initial assessment")]

myAnnot[grepl("f.20003.0.",colNm), parameter := paste("Medication",1:48,sep="_")]
myAnnot[grepl("f.20003.0.",colNm), comment := c("List of medications, up to 48 possible")]

myAnnot[colNm %in% c("f.30780.0.0"), parameter := "LDL direct (mmol/L)"]
myAnnot[colNm %in% c("f.30780.0.0"), comment := c("initial assessment")]

myAnnot[colNm %in% c(paste0("f.22009.0.",1:10)), parameter := "Genetical Principal Component"]
myAnnot[colNm %in% c(paste0("f.22009.0.",1:10)), comment := c("first 10 (out of 40)")]

myAnnot[colNm %in% c(paste0("f.22021.0.",0)), parameter := "Kinship"]
myAnnot[colNm %in% c(paste0("f.22021.0.",0)), comment := "0 - no kinship; 1 - at least one kinship; 10 - more than 10"]

myAnnot[colNm %in% c("f.21000.0.0"), parameter := "Ethnic background"]
myAnnot[colNm %in% c("f.21000.0.0"), comment := c("1001 codes for White British")]

myAnnot[grepl("f.41202.0.",colNm), parameter := paste("Diagnoses_Main_ICD10",1:80,sep="_")]
myAnnot[grepl("f.41202.0.",colNm), comment := c("List of main diagnoses, up to 80 possible")]

myAnnot[grepl("f.41204.0.",colNm), parameter := paste("Diagnoses_Sec_ICD10",1:210,sep="_")]
myAnnot[grepl("f.41204.0.",colNm), comment := c("List of secondary diagnoses, up to 210 possible")]

myAnnot[grepl("f.41270.0.",colNm), parameter := paste("Diagnoses_ICD10",1:259,sep="_")]
myAnnot[grepl("f.41270.0.",colNm), comment := c("List of any diagnoses, up to 259 possible")]


#' # Load all samples ####
#' ***
x = myAnnot[!is.na(comment),colNR]
myTab <- fread(paste0(path_UKBB,"phenotypes/ukb672224.tab"), 
               header=TRUE, sep="\t",select = x)
names(myTab)
names(myTab) = c("ID","sex",paste("medication",1:48,sep="_"),"ethnic_background","age",
                 paste("PC",1:10,sep="_"),"kinship","LDLC",
                 paste("diagnoses_main",1:80,sep="_"),paste("diagnoses_sec",1:210,sep="_"),
                 paste("diagnoses_any",1:259,sep="_"))

save(myTab, file = "../temp/06_UKBB_complete.RData")

#' # Filter samples ####
#' ***
#' Filter for white British people with no kinship
myTab = myTab[ethnic_background == 1001,]
myTab = myTab[kinship == 0,]

#' Filter for consent
ToExclude = fread(paste0(path_UKBB,"phenotypes/withdraw98032_19.txt"))
myTab = myTab[!is.element(ID,ToExclude$V1),]

#' Get statin strata
codingTable = data.table(read_xlsx("../temp/Wu_2019_SupplementalData1_modified.xlsx",sheet=1))

names(myTab)
myMeds = names(myTab)[grep("medication_",names(myTab))]
statins = codingTable[grepl("C10AA",ATC_Code),Coding]
lipidLowering = codingTable[grepl("C10A",ATC_Code),Coding]
myTab[,statin := 0]
myTab[,lipidLow := 0]

for(i in 1:length(myMeds)){
  #i=1
  myTab[get(myMeds[i]) %in% statins,statin := 1]
  myTab[get(myMeds[i]) %in% lipidLowering,lipidLow := 1]
}
myTab[,get("myMeds"):=NULL]
dim(myTab)
table(myTab$lipidLow,myTab$statin)

#' Get coronary artery disease
codingTable = fread("coding19.tsv")

myDiag1 = names(myTab)[grep("diagnoses_main",names(myTab))]
myDiag2 = names(myTab)[grep("diagnoses_sec",names(myTab))]
myDiag3 = names(myTab)[grep("diagnoses_any",names(myTab))]
myDiag = c(myDiag1,myDiag2,myDiag3)

codingTable = codingTable[grepl("I20",meaning) | grepl("I21",meaning) | grepl("I22",meaning) | 
                            grepl("I23",meaning) | grepl("I24",meaning) | grepl("I25",meaning),]
codingTable = codingTable[selectable == "Y",]
CAD_codes = codingTable[,coding]

myTab[,CAD := 1]

for(i in 1:length(myDiag)){
  #i=1
  myTab[get(myDiag[i]) %in% CAD_codes,CAD := 2]
  
}
myTab[,get("myDiag"):=NULL]
dim(myTab)

table(myTab$CAD)

#' Filter sex mismatches
sampleFile = fread(paste0(path_UKBB,"genotypes/ukb22418_c1_b0_v2_s488131.fam"))
sampleFile = sampleFile[V1 %in% myTab$ID,]
matched = match(sampleFile$V1,myTab$ID)
myTab = myTab[matched,]
table(myTab$ID == sampleFile$V1)
table(myTab$sex, sampleFile$V5)
myTab[,sex2 := sampleFile$V5]
myTab = myTab[(sex==1 & sex2==1) | (sex==0 & sex2==2),]
myTab[,sex:= NULL]
setnames(myTab,"sex2","sex")

#' Check missingness: sex, age, statins, genetics, LDLC, CAD
table(is.na(myTab$sex))
table(is.na(myTab$age))
table(is.na(myTab$statin))
table(is.na(myTab$PC_1))
table(is.na(myTab$LDLC))
table(is.na(myTab$CAD))
table(is.na(myTab$LDLC),myTab$CAD)

#' Save filtered data 
save(myTab, file = "../temp/06_UKBB_filtered.RData")

#' # Get files for GWAS ####
#' ***
#' Get LDLC traits
myTab[sex == 2 & lipidLow == 0, ldlc_fem_free := LDLC]
myTab[sex == 1 & lipidLow == 0, ldlc_mal_free := LDLC]
myTab[sex == 2 & lipidLow == 1, ldlc_fem_treated := LDLC]
myTab[sex == 1 & lipidLow == 1, ldlc_mal_treated := LDLC]

#' Get CAD traits
myTab[sex == 2 & lipidLow == 0, cad_fem_free := CAD]
myTab[sex == 1 & lipidLow == 0, cad_mal_free := CAD]
myTab[sex == 2 & lipidLow == 1, cad_fem_treated := CAD]
myTab[sex == 1 & lipidLow == 1, cad_mal_treated := CAD]

#' Check traits 
myTraits1 = names(myTab)[grep("ldlc_",names(myTab))]
myTraits2 = names(myTab)[grep("cad_",names(myTab))]

dumTab = foreach(i = 1:length(myTraits1))%do%{
  #i=5
  myTrait = myTraits1[i]
  myTab[,trait := get(myTrait)]
  x1 = summary(myTab$trait)
  x2 = sum(!is.na(myTab$trait))
  x3 = sd(myTab$trait,na.rm = T)
  res = as.numeric(x1)
  res = c(res,x2,x3)
  
  res2 = data.table(phenotype = myTrait,
                    mean = round(res[4],3),
                    sd = round(res[9],3),
                    min = round(res[1],3),
                    median = round(res[3],3),
                    max = round(res[6],3),
                    sampleSize = res[8],
                    NAs = res[7])
  res2
}
summary = rbindlist(dumTab)
summary

dumTab2 = foreach(i = 1:length(myTraits2))%do%{
  #i=1
  myTrait = myTraits2[i]
  myTab[,trait := get(myTrait)]
  x1 = table(myTab$trait)
  x2 = sum(!is.na(myTab$trait))
  res = as.numeric(x1)
  res = c(res,x2)
  
  res2 = data.table(phenotype = myTrait,
                    controls = res[1],
                    cases = res[2],
                    sampleSize = res[3])
  res2
}
summary2 = rbindlist(dumTab2)
summary2

#' Create files for PLINK
myPhenFile1 = copy(myTab)
myPhenFile1 = myPhenFile1[,c(1,1,20:27)]
names(myPhenFile1)[1:2] = c("FID","IID")

myCovarFile1 = copy(myTab)
myCovarFile1 = myCovarFile1[,c(1,1,3:13)]
names(myCovarFile1)[1:2] = c("FID","IID")

write.table(myPhenFile1,file="../temp/06_PhenoFile_UKBB_ldlc_cad.txt",
            sep=" ",col.names = T,row.names = F,quote=F)

write.table(myCovarFile1,file="../temp/06_CovarFile_UKBB_ldlc_cad.txt",
            sep=" ",col.names = T,row.names = F,quote=F)

#' Create PLINK call 

mycall1 = paste0("plink2", 
                 " --bgen ",path_UKBB,"genotypes-imputed/ukb22828_c1_b0_v3.bgen",
                 " --sample ",path_UKBB,"genotypes-imputed/ukb22828_c1_b0_v3_s487160.sample",
                 " --glm hide-covar firth-fallback",
                 " cols=chrom,pos,ref,alt,firth,test,nobs,machr2,a1freq,a1freqcc,a1countcc,orbeta,se,ci,tz,p",
                 " --pheno ../temp/06_PhenoFile_UKBB_ldlc_cad.txt", 
                 " --covar ../temp/06_CovarFile_UKBB_ldlc_cad.txt",
                 " --from-bp 55005647 --to-bp 56005647 --chr 1 --threads 20",
                 " --out ",path_UKBB,"hpc-work/230831_LDLC_CAD_stratified")
mycall1

#' Run mycall1 on HPC! Not done via R! (slurm scheduler)
#' 
#' # Load data ####
#' ***
ToDoFiles = list.files(path = paste0(path_UKBB,"hpc-work/"), 
                       pattern = "230831_LDLC_CAD_stratified.")
#' Remove log file
ToDoFiles = ToDoFiles[-9]

dumTab = foreach(i = 1:length(ToDoFiles))%do%{
  # i=1
  file = ToDoFiles[i]
  dat = fread(paste0(path_UKBB,"hpc-work/",file))
  dat[,pheno := gsub("230831_LDLC_CAD_stratified.","",file)]
  dat[,pheno := gsub(".glm.linear","",pheno)]
  dat[,pheno := gsub(".glm.logistic.hybrid","",pheno)]
  dat[,pheno := gsub("mal","males",pheno)]
  dat[,pheno := gsub("fem","females",pheno)]
  dat
}
dat = rbindlist(dumTab,fill=T)

mySNPs1 = c("rs11591147","rs11583680","rs693668","rs2495491")
mySNPs2 = c("rs9436961","rs11206510","rs79224150","rs2479409","rs11591147","rs150119739","rs693668","rs149093957","rs505151")
mySNPs = unique(c(mySNPs1,mySNPs2))

dat2 = copy(dat)
dat2 = dat2[ID %in% mySNPs,]
names(dat2)[1] = "CHR"
dat2[is.na(BETA), BETA := log(as.numeric(OR))]
dat2[is.na(SE), SE := `LOG(OR)_SE`]
dat2[is.na(T_STAT), T_STAT := BETA/SE]
dat2[, NR_cases := round(A1_CASE_CT/2/A1_CASE_FREQ,0)]
dat2[, NR_controls := round(A1_CTRL_CT/2/A1_CTRL_FREQ,0)]

#' # Meta-analysis ####
#' ***
ToDoList1 = data.table(SNP = rep(mySNPs,2),
                       trait = rep(c("ldlc","cad"),each=length(mySNPs)),
                       subgroup = "free",
                       strata1 = "females_free",
                       strata2 = "males_free")

ToDoList2 = data.table(SNP = rep(mySNPs,2),
                       trait = rep(c("ldlc","cad"),each=length(mySNPs)),
                       subgroup = "treated",
                       strata1 = "females_treated",
                       strata2 = "males_treated")

ToDoList3 = data.table(SNP = rep(mySNPs,2),
                       trait = rep(c("ldlc","cad"),each=length(mySNPs)),
                       subgroup = "females",
                       strata1 = "females_free",
                       strata2 = "females_treated")

ToDoList4 = data.table(SNP = rep(mySNPs,2),
                       trait = rep(c("ldlc","cad"),each=length(mySNPs)),
                       subgroup = "males",
                       strata1 = "males_free",
                       strata2 = "males_treated")

ToDoList = rbind(ToDoList1,ToDoList2,ToDoList3,ToDoList4)

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  myPhenos = c(paste(myRow$trait,myRow$strata1,sep="_"),paste(myRow$trait,myRow$strata2,sep="_"))
  
  dat3 = copy(dat2)
  dat3 = dat3[ID == myRow$SNP & pheno %in% myPhenos,]
  samplesize = sum(dat3$OBS_CT)
  if(myRow$trait == "cad"){
    NR_case = sum(dat3$NR_cases)
    NR_control = sum(dat3$NR_controls)    
  } 
  weight = dat3$OBS_CT/samplesize
  mean_EAF = round(sum(weight * dat3$A1_FREQ),6)
  mean_R2 = round(sum(weight * dat3$MACH_R2),6)
  
  mod = metagen(TE = dat3[,BETA], 
                seTE = dat3[,SE],
                studlab = dat3[,pheno])
  
  res = dat3[1,]
  res[,A1_FREQ := mean_EAF]
  res[,MACH_R2 := mean_R2]
  res[,OBS_CT := samplesize]
  res[,BETA := mod$TE.fixed]
  res[,SE := mod$seTE.fixed]
  res[,T_STAT := BETA/SE]
  res[,P := mod$pval.fixed]
  res[,pheno := paste(myRow$trait,myRow$subgroup,sep="_")]
  res[,I2 := mod$I2]
  if(myRow$trait == "cad"){
    res[,NR_cases := NR_case]
    res[,NR_controls := NR_control]    
  } 
  
  res
  
}
dat4 = rbindlist(dumTab)
dat4

#' # Combine and save ####
#' ***
dat5 = rbind(dat2,dat4,fill=T)
myNames = names(dat5)[c(20,1:6,9,12,15,24,25,21,22,23,19,26)]
colsOut<-setdiff(colnames(dat5),myNames)
dat5[,get("colsOut"):=NULL]
setcolorder(dat5,myNames)
setorder(dat5,POS,pheno)
dat5[,pheno := gsub("ldlc","LDLC",pheno)]
dat5[,pheno := gsub("cad","CAD",pheno)]

dat6 = copy(dat5)
dat6 = dat6[ID %in% mySNPs1,]
save(dat6,file="../temp/06_LDLC_CAD_assoc_stratified_indepSignals.RData")

dat7 = copy(dat5)
dat7 = dat7[ID %in% mySNPs1,]
save(dat7,file="../temp/06_LDLC_CAD_assoc_stratified_moreSNPs.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
