#' ---
#' title: "LD Score Regression Analyses"
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
#' I want to estimate the heritability of PCSK9 in all subsets, and then perform a correlation test within my PCSK9 traits and with lipids. 
#' 
#' Within this R script (04_1), I will extract the GWAS summary statistics and save them in the format necessary for LDSC. Then, I generate all necessary LDSC commands, which will be executed in a putty shell, not here (reticulate connection with python sometimes unstable). 
#' 
#' In script 04_2, I read in the LDSC log files for heritability and genetic correlation and generate some plots. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_angmar.R")
.libPaths()
setwd(paste0(projectpath_main,"/scripts/"))

source("../helperFunctions/LDSC_PrepInputData_PCSK9.R")
source("../helperFunctions/LDSC_PrepInputData_GLGC.R")


#' # Step 1: Munge PCSK9 data  ####
#' ***
ToDoList = data.table(statistic = list.files(path = "../data/",pattern = "txt.gz"))
ToDoList[,statistic_path := paste0("../data/",statistic)]
ToDoList[,pheno := gsub("SumStat_","",statistic)]
ToDoList[,pheno := gsub("_230120.txt.gz","",pheno)]
ToDoList

cor_comb1 = LDSC_PrepInputData_PCSK9(path_data_trait1 = ToDoList$statistic_path[1],
                                       path_data_trait2 = ToDoList$statistic_path[5],
                                       out_dir = "../temp/04_LDSC/",
                                       trait1_name = ToDoList$pheno[1],
                                       trait2_name = ToDoList$pheno[5],
                                       trait1_type = "cont",
                                       trait2_type = "cont",
                                       cortest_method = "pearson",
                                       path_ldsc_data = path_ldsc_data)

cor_free = LDSC_PrepInputData_PCSK9(path_data_trait1 = ToDoList$statistic_path[2],
                                   path_data_trait2 = ToDoList$statistic_path[6],
                                   out_dir = "../temp/04_LDSC/",
                                   trait1_name = ToDoList$pheno[2],
                                   trait2_name = ToDoList$pheno[6],
                                   trait1_type = "cont",
                                   trait2_type = "cont",
                                   cortest_method = "pearson",
                                   path_ldsc_data = path_ldsc_data)

cor_treated = LDSC_PrepInputData_PCSK9(path_data_trait1 = ToDoList$statistic_path[3],
                                      path_data_trait2 = ToDoList$statistic_path[7],
                                      out_dir = "../temp/04_LDSC/",
                                      trait1_name = ToDoList$pheno[3],
                                      trait2_name = ToDoList$pheno[7],
                                      trait1_type = "cont",
                                      trait2_type = "cont",
                                      cortest_method = "pearson",
                                      path_ldsc_data = path_ldsc_data)

cor_comb2 = LDSC_PrepInputData_PCSK9(path_data_trait1 = ToDoList$statistic_path[4],
                                       path_data_trait2 = ToDoList$statistic_path[8],
                                       out_dir = "../temp/04_LDSC/",
                                       trait1_name = ToDoList$pheno[4],
                                       trait2_name = ToDoList$pheno[8],
                                       trait1_type = "cont",
                                       trait2_type = "cont",
                                       cortest_method = "pearson",
                                       path_ldsc_data = path_ldsc_data)

cor_tabs_PCSK9 = rbind(cor_comb1,cor_comb2,cor_free,cor_treated)
cor_tabs_PCSK9
save(cor_tabs_PCSK9, file = "../results/04_PearsonsCorrelation.RData")
load("../results/04_PearsonsCorrelation.RData")

#' # Step 2: Munge lipid data  ####
#' ***
ToDoList = data.table(statistic = list.files(path = path_lipids,pattern = "_N_1.gz", recursive = TRUE))
ToDoList = ToDoList[!grepl(".tbi",statistic)]

ToDoList[,pheno := gsub(".*_AFR_EAS_EUR_HIS_SAS_","",statistic)]
ToDoList[,pheno := gsub("_with_N_1.gz","",pheno)]
ToDoList[,pheno := gsub("_INV","",pheno)]
ToDoList[,statistic_path := paste0(path_lipids,statistic)]
ToDoList

cor_HDL = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[1],
                                      path_data_trait2 = ToDoList$statistic_path[2],
                                      out_dir = "../temp/04_LDSC/",
                                      trait1_name = ToDoList$pheno[1],
                                      trait2_name = ToDoList$pheno[2],
                                      trait1_type = "cont",
                                      trait2_type = "cont",
                                      cortest_method = "pearson",
                                      path_ldsc_data = path_ldsc_data)

cor_LDL = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[3],
                                        path_data_trait2 = ToDoList$statistic_path[4],
                                        out_dir = "../temp/04_LDSC/",
                                        trait1_name = ToDoList$pheno[3],
                                        trait2_name = ToDoList$pheno[4],
                                        trait1_type = "cont",
                                        trait2_type = "cont",
                                        cortest_method = "pearson",
                                        path_ldsc_data = path_ldsc_data)

cor_logTG = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[5],
                                        path_data_trait2 = ToDoList$statistic_path[6],
                                        out_dir = "../temp/04_LDSC/",
                                        trait1_name = ToDoList$pheno[5],
                                        trait2_name = ToDoList$pheno[6],
                                        trait1_type = "cont",
                                        trait2_type = "cont",
                                        cortest_method = "pearson",
                                        path_ldsc_data = path_ldsc_data)

cor_nonHDL = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[7],
                                        path_data_trait2 = ToDoList$statistic_path[8],
                                        out_dir = "../temp/04_LDSC/",
                                        trait1_name = ToDoList$pheno[7],
                                        trait2_name = ToDoList$pheno[8],
                                        trait1_type = "cont",
                                        trait2_type = "cont",
                                        cortest_method = "pearson",
                                        path_ldsc_data = path_ldsc_data)

cor_TC = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[9],
                                        path_data_trait2 = ToDoList$statistic_path[10],
                                        out_dir = "../temp/04_LDSC/",
                                        trait1_name = ToDoList$pheno[9],
                                        trait2_name = ToDoList$pheno[10],
                                        trait1_type = "cont",
                                        trait2_type = "cont",
                                        cortest_method = "pearson",
                                        path_ldsc_data = path_ldsc_data)

cor_LDL_TC = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[12],
                                 path_data_trait2 = ToDoList$statistic_path[15],
                                 out_dir = "../temp/04_LDSC/",
                                 trait1_name = ToDoList$pheno[12],
                                 trait2_name = ToDoList$pheno[15],
                                 trait1_type = "cont",
                                 trait2_type = "cont",
                                 cortest_method = "pearson",
                                 path_ldsc_data = path_ldsc_data)

cor_HDL_nonHDL = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[11],
                                     path_data_trait2 = ToDoList$statistic_path[14],
                                     out_dir = "../temp/04_LDSC/",
                                     trait1_name = ToDoList$pheno[11],
                                     trait2_name = ToDoList$pheno[14],
                                     trait1_type = "cont",
                                     trait2_type = "cont",
                                     cortest_method = "pearson",
                                     path_ldsc_data = path_ldsc_data)

cor_TG = LDSC_PrepInputData_GLGC(path_data_trait1 = ToDoList$statistic_path[13],
                                         out_dir = "../temp/04_LDSC/",
                                         trait1_name = ToDoList$pheno[13],
                                         trait1_type = "cont",
                                         path_ldsc_data = path_ldsc_data)

cor_tabs_lipids = rbind(cor_HDL,cor_LDL,cor_logTG,cor_nonHDL,cor_TC,cor_LDL_TC,cor_HDL_nonHDL,cor_TG,use.names=T,fill=TRUE)
cor_tabs_lipids
cortabs = rbind(cor_tabs_PCSK9,cor_tabs_lipids,fill=TRUE)
save(cortabs, file = "../results/04_PearsonsCorrelation.RData")

#' # Step 3: Generate Calls ####
#' ***
#' ## Munge Calls 
#' 
ToDoList = data.table(statistic = list.files(path = "../temp/04_LDSC/",pattern = "_input", recursive = TRUE))
ToDoList[,outfn := gsub("_input","_munge",statistic)]
ToDoList[,outfn := gsub(".txt.gz","",outfn)]
ToDoList
mungeCalls = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1

  myCall = paste0("./munge_sumstats.py",
                  " --sumstats ",projectpath_main,"temp/04_LDSC/",ToDoList$statistic[i],
                  " --ignore beta",
                  " --out ",projectpath_main,"temp/04_LDSC/",ToDoList$outfn[i],
                  " --merge-alleles ",path_ldsc_data)
  myCall
}
mungeCalls = unlist(mungeCalls)

#' ## Heritability Calls 
#' 
ToDoList[,outfn_heritab := gsub("LDSC_munge_","",outfn)]
ToDoList[,outfn_heritab := paste0(outfn_heritab,"_heritab")]
ToDoList
heritabCalls = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  
  myCall = paste0("./ldsc.py",
                  " --h2 ",projectpath_main,"temp/04_LDSC/",ToDoList$outfn[i],".sumstats.gz",
                  " --ref-ld-chr ",path_ldsc_data,
                  " --w-ld-chr ",path_ldsc_data,
                  " --out ",projectpath_main,"temp/04_LDSC/",ToDoList$outfn_heritab[i])
  myCall
}
heritabCalls = unlist(heritabCalls)
heritabCalls = gsub("w_hm3.snplist","",heritabCalls)

#' ## GenCor Calls 
#' 
#' Within PCSK9: females vs males and free vs treated (2 tests)
#' 
#' With lipids: GWAMA2 PCSK9 vs lipids (30 tests)
#' 
ToDoList2 = data.table(trait1 = c(ToDoList$outfn[15],ToDoList$outfn[16],
                                  rep(ToDoList$outfn[15],5),
                                  rep(ToDoList$outfn[19],5),
                                  rep(ToDoList$outfn[16],5),
                                  rep(ToDoList$outfn[20],5)),
                       trait2 = c(ToDoList$outfn[19],ToDoList$outfn[20],
                                  ToDoList$outfn[c(2,5,8,11,22)],
                                  ToDoList$outfn[c(3,6,9,12,23)],
                                  ToDoList$outfn[c(1,4,7,10,21)],
                                  ToDoList$outfn[c(1,4,7,10,21)]))
ToDoList2[,outfn_gencor := paste0(trait1,"_VS_",trait2)]
ToDoList2[,outfn_gencor := gsub("LDSC_munge_","",outfn_gencor)]
ToDoList2[,outfn_gencor := gsub("_GLGC","",outfn_gencor)]
ToDoList2[,outfn_gencor := gsub("_female_adjusted","",outfn_gencor)]
ToDoList2[,outfn_gencor := gsub("_male_adjusted","",outfn_gencor)]
ToDoList2[1,outfn_gencor := "PCSK9_female_VS_PCSK9_male_adjusted"]
ToDoList2[2,outfn_gencor := "PCSK9_female_VS_PCSK9_male_free"]
ToDoList2[,outfn_gencor := paste0(outfn_gencor,"_gencor")]
ToDoList2

genCorCalls = foreach(i=1:dim(ToDoList2)[1])%do%{
  #i=1
  
  myCall = paste0("./ldsc.py",
                  " --rg ",projectpath_main,"temp/04_LDSC/",ToDoList2$trait1[i],".sumstats.gz,",
                  projectpath_main,"temp/04_LDSC/",ToDoList2$trait2[i],".sumstats.gz",
                  " --ref-ld-chr ",path_ldsc_data,
                  " --w-ld-chr ",path_ldsc_data,
                  " --out ",projectpath_main,"temp/04_LDSC/",ToDoList2$outfn_gencor[i])
  myCall
}
genCorCalls = unlist(genCorCalls)
genCorCalls = gsub("w_hm3.snplist","",genCorCalls)

#' ## Save Calls 
Header = c("#!/bin/bash",
           "",
           "#Activate LDSC python environment",
           paste0("cd ",path_phyton),
           "source activate ldsc",
           "cd ../../ldsc",
           "")
Header = gsub("python","",Header)

shell = c(Header,
          "# Munge Sumstats",
          mungeCalls,
          "",
          "# Estimate heritability per trait",
          heritabCalls,
          "",
          "# Estimate genetic Correlation",
          genCorCalls,
          "")

write.table(shell, "../temp/04_LDSC/LDSC_shell_commands2.sh", 
            row.names=F, quote=F, col.names=F)



#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
