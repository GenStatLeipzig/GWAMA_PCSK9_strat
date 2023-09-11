#############################
# this is a template source file
# please change all paths accordingly
#############################

#############################
# Working directory
#############################
projectpath_main = "/PATH/TO/RPROJECTS/GWAMA_PCSK9_strat/"

#############################
# R library and R packages
#############################
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(toolboxH))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(meta)) 
suppressPackageStartupMessages(library(metafor)) 
suppressPackageStartupMessages(library(dplyr)) 
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(MendelianRandomization))

#############################
# Other Tools 
#############################
# GCTA paths
path_gcta = "/PATH/TO/TOOLS/gcta/gcta_1.93.2beta/gcta64"
path_gcta_RefData = "/PATH/TO/TOOLS/gcta/REFERENCE_DATA/LIFE_COMBINED"

# PLINK2 path
path_plink2 = "/PATH/TO/TOOLS/plink2.0/20210203/unix_64/plink2"

#############################
# Downloaded data sets 
#############################
path_GTExv8 = "/PATH/TO/DATA/GTEx_v8/"
path_GTExv8_sexStrat = "/PATH/TO/DATA/GTEx_v8_sexStrat/"
path_1000Genomes = "/PATH/TO/DATA/1000genomes_phase3_vs5/plinkformat/EUR"
path_1000Genomes_snps_fn<- "/PATH/TO/DATA/1000genomes_phase3_vs5/plinkformat/EUR/autoANDchrXnonPAR._snps.txt"
path_lipids = "/PATH/TO/DATA/SummaryStatistics/GlobalLipidGeneticConsortium/sex_specific_summary_stats/"
path_otherGWAS = "/PATH/TO/DATA/SummaryStatistics/"
path_CAD ="/PATH/TO/DATA/SummaryStatistics/CAD/"

#############################
# Data from in-house pipeline or UKBB (not executable without access to UKBB or Genstat-Pipeline) 
#############################
path_UKBB = ""
path_GenStatPipeline = ""
path_RAPlot_function = ""
