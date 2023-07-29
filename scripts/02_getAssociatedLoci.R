#' ---
#' title: "Get associated loci"
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
#' Define loci as in Wuttke publication. I modify a script of Katrin Horn to do this. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/getSmallestDist.R")

#' # Load data ####
#' ***
#' Load data and filter for suggestive significant and valid SNPs. Finally, merge all significant SNPs into one object. 
ToDoList = data.table(NR = 1:8)

ToDoList[,statistic := list.files(path = "../data/",pattern = "PCSK9")]
ToDoList[,statistic_path := paste0("../data/",statistic)]

ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_23.*","",pheno)]

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("Working on data ",myRow$pheno)
  data = fread(myRow$statistic_path)
  
  message("    SNPs before filter: ",dim(data)[1])
  data = data[invalidAssocs  == F & pval<=1e-6,]
  message("    SNPs after filter: ",dim(data)[1])
  
  data[,chrPosPheno := paste(chr,bp_hg19,phenotype,sep=":")]
  data
}

result.1 = rbindlist(dumTab)
table(duplicated(result.1$chrPosPheno))
table(result.1$phenotype,result.1$chr)

#' # Get Top SNP per 500 kb range ###
#' ***
#' I want to reduce the list of associated SNPs by assigning SNPs to the top SNP of the region per phenotype. To do this, I will order the SNPs by their p-value, choose the SNP with lowest p-value as lead SNP, and assign all SNPs within 500 kb around the lead SNP to this region. This will be repeated until no SNPs can be assign, as the minimal pairwise distance is larger than 500 kb.  
#' 
subset1 = unique(result.1$phenotype)
subset1

result.2 = foreach(s1 = subset1) %do% {
  # s1 = subset1[1]
  subdata = copy(result.1)
  subdata = subdata[phenotype == s1, ]
  subset2 = unique(subdata$chr)
  
  result.22 = foreach(s2 = subset2) %do% {
    # s2 = subset2[1]
    subdata2 = copy(subdata)
    subdata2 = subdata2[chr == s2, ]
    
    setkey(subdata2, bp_hg19)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, NR_SNPs := 0]
    }else{
      subdata2[, keep := NA]
      subdata2[, NR_SNPs := as.numeric(NA)]
      
      smallestDist = getSmallestDist(subdata2[, bp_hg19])
      while(smallestDist < 500000) {
        minP = min(subdata2[is.na(keep), pval])
        mybp_hg19 = subdata2[minP == pval & is.na(keep), bp_hg19]
        if(length(mybp_hg19)>1){
          mybp_hg19 = mybp_hg19[1]
        }
        subdata2[bp_hg19 == mybp_hg19, keep := T]
        
        #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
        myFilt = (subdata2[, bp_hg19] < (mybp_hg19 - 500000)) | 
          (subdata2[, bp_hg19] > (mybp_hg19 + 500000)) | 
          subdata2[, keep] 
        myFilt[is.na(myFilt)] = FALSE
        subdata2 = subdata2[myFilt == TRUE, ]
        
        subdata2[bp_hg19 == mybp_hg19, NR_SNPs := sum(myFilt==F)]
        smallestDist = getSmallestDist(subdata2[, bp_hg19])
      }
      
      #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
      subdata2[is.na(keep), NR_SNPs := 0]
      subdata2[is.na(keep), keep := TRUE]
    }
    
    subdata2
  }
  result.22 = rbindlist(result.22)
  table(result.22$phenotype)
  table(result.22$chr)
  result.22
}
result.2 = rbindlist(result.2)
table(result.2$phenotype,result.2$chr)

#add regions
result.2[, region_start := bp_hg19 - 500000]
result.2[, region_end := bp_hg19 + 500000]
result.2

#' # Locus collapsing by phenotype ###
#' ***
#' I want to collapse loci if their regions are overlapping, although their base positions differ more then 500 kb. This could be the case in large associated regions. 
#' 

result.3 = foreach(s1 = subset1) %do% {
  # s1 = subset1[1]
  subdata = copy(result.2)
  subdata = subdata[phenotype == s1, ]
  subset2 = unique(subdata$chr)
  
  result.32 = foreach(s2 = subset2) %do% {
    # s2 = subset2[1]
    subdata2 = copy(subdata)
    subdata2 = subdata2[chr == s2, ]
    
    setkey(subdata2, bp_hg19)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, region := 1]
    }else{
      subdata2[, keep := NA]
      subdata2[, region := as.numeric(NA)]
      
      subdata2[1, region := 1]
      foreach(l = c(2:nrow(subdata2))) %do% {
        if (subdata2[l, region_start] <= subdata2[l-1, region_end]) {
          subdata2[l, region := subdata2[l-1, region]]
        } else { subdata2[l, region := subdata2[l-1, region] + 1]}
      }
      
      #keep best SNPs per region
      allRegions = unique(subdata2[, region])
      
      subresult = foreach(r = allRegions) %do% {
        subdata.2 = copy(subdata2)
        subdata.2 = subdata.2[region == r, ]
        minP = min(subdata.2[, pval])
        subdata.2[pval == minP, region_start := min(subdata.2[, region_start])]
        subdata.2[pval == minP, region_end := max(subdata.2[, region_end])]
        subdata.2[pval == minP, NR_SNPs := sum(subdata.2[,NR_SNPs]) + nrow(subdata.2) - 1]
        subdata.2[pval == minP, ]
      }
      subresult = rbindlist(subresult)
      subresult[,region := paste(s2,region,sep="::")]
      subresult
    }
    
  }
  result.32 = rbindlist(result.32,fill = T)
  table(result.32$phenotype)
  table(result.32$chr)
  result.32
  
}
result.3 = rbindlist(result.3,fill = T)
table(result.3$phenotype,result.3$chr)
result.3[, region := c(1:nrow(result.3))]
result.3

names(result.3)
cols2keep = names(result.3)[c(22,20,21,16,1:9,19,10:13)]
colsOut = setdiff(colnames(result.3), cols2keep)
result.3[, get("colsOut") := NULL]
setcolorder(result.3, cols2keep)
setorder(result.3,chr,bp_hg19)
result.3
result.3[,.N,by=phenotype]

save(result.3,file="../results/02_LociPerPhenotype.RData")

#' # Locus collapsing over all phenotype ###
#' *** 
#' I want to know which loci appear in multiple settings, and how these collapsed loci would look like. 
#' 
subdata = copy(result.2)
setorder(subdata, chr, bp_hg19)
subset2 = unique(subdata$chr)

result.4 = foreach(s2 = subset2) %do% {
  # s2 = subset2[1]
  subdata2 = copy(subdata)
  subdata2 = subdata2[chr == s2, ]
  
  setkey(subdata2, bp_hg19)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, region := 1]
    subdata2[, other_phenotypes := ""]
    subdata2[, region := paste0("chr",s2,"::",region)]
    
  }else{
    subdata2[, keep := NA]
    subdata2[, region := as.numeric(NA)]
    
    subdata2[1, region := 1]
    foreach(l = c(2:nrow(subdata2))) %do% {
      if (subdata2[l, region_start] <= subdata2[l-1, region_end]) {
        subdata2[l, region := subdata2[l-1, region]]
      } else { subdata2[l, region := subdata2[l-1, region] + 1]}
    }
    
    #keep best SNPs per region
    allRegions = unique(subdata2[, region])
    
    subresult = foreach(r = allRegions) %do% {
      subdata.2 = copy(subdata2)
      subdata.2 = subdata.2[region == r, ]
      minP = min(subdata.2[, pval])
      subdata.2[pval == minP, region_start := min(subdata.2[, region_start])]
      subdata.2[pval == minP, region_end := max(subdata.2[, region_end])]
      subdata.2[pval == minP, NR_SNPs := sum(subdata.2[,NR_SNPs]) + nrow(subdata.2) - 1]
      subdata.2[pval == minP, other_phenotypes := paste(unique(subdata.2[pval != minP, phenotype]), collapse = " | ")]
      subdata.2[pval == minP, ]
    }
    subresult = rbindlist(subresult)
    subresult[,region := paste0("chr",s2,"::",region)]
    subresult
  }
}
  
result.4 = rbindlist(result.4,fill = T)
result.4[,.N,by=phenotype]
table(result.4$chr)
result.4

#' There is a total of 29 loci associated with at least one SNP. But I want to loci to be supported by at least two other variants!
#' 
result.4[NR_SNPs>1,]
#' 
#' There is a total of 11 loci associated with sufficient support. 
#' 

names(result.4)
cols2keep = names(result.4)[c(22,20,21,16,23,1:9,19,10:13)]
colsOut = setdiff(colnames(result.4), cols2keep)
result.4[, get("colsOut") := NULL]
setcolorder(result.4, cols2keep)
result.4
setorder(result.4,chr,bp_hg19)

save(result.4,file="../results/02_LociOverallPhenotypes.RData")

result.5 = copy(result.4)
result.5 = result.5[NR_SNPs >1]
result.5[,candidateGene := c("PCSK9","APOB","KHDRBS2",
                             "PRKAG2","ALOX5","JMJD1C","FADS1","SLCO1B3","NOS1",
                             "HP","TM6SF2")]
result.5
save(result.5,file="../results/02_LociOverallPhenotypes_filtered.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

