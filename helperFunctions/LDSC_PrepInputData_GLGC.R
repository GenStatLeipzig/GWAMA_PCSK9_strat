#' @title Get summary Statistics in right format for LD Score regression
#' @description In this function, GWAMA summary statistics are loaded and prepped for LDSC analyses. The data is filtered for invalid associations, filtered for rsID matching the LDSC data and then saved as text file. Finally, the text file is used in the python based LDSC munge call. This can be done with two data objects at a time, allowing to estimate the raw correlation, e.g. necessary in interaction analyses of two strata.
#' @param path_data_trait1 complete path and file name of (first) GWAMA object
#' @param path_data_trait2 optional, complete path and file name of second GWAMA object, Default: 'none'
#' @param out_dir complete path to directory used to store all output files
#' @param trait1_name phenotype name of (first) trait
#' @param trait2_name optional, phenotype name of second trait, Default: 'none'
#' @param trait1_type phenotype type of (first) trait. Must be either 'cont' or 'bin' for continuous or binary traits, Default: 'cont'
#' @param trait2_type optional, phenotype type of second trait. Must be either 'cont' or 'bin' for continuous or binary traits, Default: 'cont'
#' @param cortest_method method used in correlation test. Only relevant when two traits are given as input. Must be either 'pearson', 'spearman' or 'kendall', Default: 'pearson'
#' @param path_ldsc_data path to LDSC example data
#' @return This function returns a data.table containing some trait and SNP information, and the correlation results:
#'
#' * phenoi = name of trait i
#' * nSamples_phenoi = sample size for trait i
#' * nSNPs_phenoi_total = SNP number in step10 object
#' * nSNPs_phenoi_filt1 = SNP number after excluding invalid SNPs due to MAF, info, I2 or triallelic SNPs
#' * nSNPs_phenoi_filt2 = SNP number after filtering for matching SNP IDs and alleles to LDSC SNP data
#' * cor_filt1 = correlation estimate and p-value between betas of trait 1 and 2 after excluding invalid SNPs
#' * cor_filt1_sig =  correlation estimate and p-value between betas of trait 1 and 2 after excluding invalid SNPs and filtering for significant SNPs (p<0.01 for at least one of the traits)
#' * cor_filt2 = correlation estimate and p-value between betas of trait 1 and 2 after filtering for matching SNPs
#' * cor_filt2 = correlation estimate and p-value between betas of trait 1 and 2 after filtering for matching and significant SNPs (p<0.01 for at least one of the traits)
#'
#' In addition, two files are generated in the output directory:
#'
#' * LDSC_input_trait.txt.gz: Summary statistics after filtering for matching SNPs. Containing the columns "chr","pos","snpid" (rsID),"A1" (effect allele),"A2" (other allele),"eaf","info","beta" (beta_score),"se" (se_score),"Zscore" (beta/se),"N","Pvalue" (p_score).
#' * LDSC_munge_trait.sumstats.gz and log: Summary statistics and log file resulted by the LDSC munge call. Containing the columns "SNP", "A1", "A2", "N", "Z".
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{setattr}}, \code{\link[data.table]{copy}}, \code{\link[data.table]{fread}}, \code{\link[data.table]{setcolorder}}, \code{\link[data.table]{setorder}}, \code{\link[data.table]{fwrite}}, \code{\link[data.table]{data.table-package}}
#'  \code{\link[R.utils]{compressFile}}
#'  \code{\link[reticulate]{use_python}}, \code{\link[reticulate]{conda-tools}}, \code{\link[reticulate]{py_config}}
#' @rdname step1_munge_step10object
#' @export
#' @importFrom data.table setnames copy fread setcolorder setorder fwrite data.table
#' @importFrom R.utils gzip
#' @importFrom reticulate use_python conda_list use_condaenv py_config
LDSC_PrepInputData_GLGC = function(path_data_trait1,
                                    path_data_trait2= "none",
                                    out_dir,
                                    trait1_name,
                                    trait2_name = "none",
                                    trait1_type = "cont",
                                    trait2_type = "cont",
                                    cortest_method = "pearson",
                                    path_ldsc_data){
   
  # Step 0: track time
  time0<-Sys.time()
  
  # Step 1: load data for first trait
  message("Loading data from ",path_data_trait1)
  erg1 = fread(path_data_trait1)
  n1_1 = dim(erg1)[1]
  
  # Step 2: filter invalid SNP association and triallelic SNPs
  message("Filtering data, starting with ",dim(erg1)[1]," SNPs ...")
  
  erg1 = erg1[ POOLED_ALT_AF >= 0.01,]
  message("After filtering invalid associations (maf, info, I2, k), ",dim(erg1)[1]," SNPs remain ...")
  
  erg1[,chrPos := paste(CHROM,POS_b37,sep=":")]
  dupPos = erg1[duplicated(chrPos),chrPos]
  erg1 = erg1[!is.element(chrPos,dupPos)]
  message("After filtering triallelic variants, ",dim(erg1)[1]," SNPs remain ...")
  
  # Step 3: prep colnames (beta, se, p-value, z score, N)
  data.table::setnames(erg1,"METAL_Effect","beta")
  data.table::setnames(erg1,"METAL_StdErr","se")
  data.table::setnames(erg1,"METAL_Pvalue","Pvalue")
  erg1[,pheno := trait1_name]
  erg1[,Zscore:= beta/se]
  
  # if(trait1_type == "bin"){
  #   dummy<-erg1$N
  #   dummy3<-gsub(" [(].*","",dummy)
  #   erg1[,N:=as.numeric(dummy3)]
  # }else{
  #   erg1[,N:=as.numeric(nSamples)]
  # }
  
  erg1_p1 = data.table::copy(erg1)
  nSamples_pheno1 = max(erg1_p1$N)
  n1_2 = dim(erg1_p1)[1]
  
  # Step 4: load data for second trait if there is one
  if(path_data_trait2 != "none"){
    # Step 4.1: load data for first trait
    message("Working on second phenotype ....")
    message("     Loading data from ",path_data_trait2)
    erg1 = fread(path_data_trait2)
    n2_1 = dim(erg1)[1]
    
    # Step 4.2: filter invalid SNP association and triallelic SNPs
    message("     Filtering data, starting with ",dim(erg1)[1]," SNPs ...")
    erg1 = erg1[ POOLED_ALT_AF >= 0.01,]
    message("     After filtering invalid associations (maf, info, I2, k), ",dim(erg1)[1]," SNPs remain ...")
    erg1[,chrPos := paste(CHROM,POS_b37,sep=":")]
    dupPos = erg1[duplicated(chrPos),chrPos]
    erg1 = erg1[!is.element(chrPos,dupPos)]
    message("     After filtering triallelic variants, ",dim(erg1)[1]," SNPs remain ...")
    
    # Step 4.3: prep colnames (beta, se, p-value, z score, N)
    data.table::setnames(erg1,"METAL_Effect","beta")
    data.table::setnames(erg1,"METAL_StdErr","se")
    data.table::setnames(erg1,"METAL_Pvalue","Pvalue")
    erg1[,pheno := trait2_name]
    erg1[,Zscore:= beta/se]
    
    # if(trait2_type == "bin"){
    #   dummy<-erg1$N
    #   dummy3<-gsub(" [(].*","",dummy)
    #   erg1[,N:=as.numeric(dummy3)]
    # }else{
    #   erg1[,N:=as.numeric(nSamples)]
    # }
    
    erg1_p2 = data.table::copy(erg1)
    nSamples_pheno2 = max(erg1_p2$N)
    n2_2 = dim(erg1_p2)[1]
    
    # Step 4.4: get row correlation
    dum1 = data.table::copy(erg1_p1)
    dum2 = data.table::copy(erg1_p2)
    dum1 = dum1[is.element(chrPos,dum2$chrPos),]
    dum2 = dum2[is.element(chrPos,dum1$chrPos),]
    matched = match(dum1$chrPos,dum2$chrPos)
    dum2 = dum2[matched,]
    cor1 = cor.test(x = dum1$beta, y= dum2$beta,method = cortest_method)
    filt = dum1$pvalue_neg_log10_GC>2 | dum2$pvalue_neg_log10_GC>2
    dum1 = dum1[filt,]
    dum2 = dum2[filt,]
    cor1_sig = cor.test(x = dum1$beta, y= dum2$beta,method = cortest_method)
    
  }else{
    n2_1 = 0
    nSamples_pheno2 = 0
    n2_2 = 0
  }
  
  # Step 5: filter for SNPs in LDSC SNPlist
  erg1_p1[,snpid := rsID]
  SNPlist<-data.table::fread(path_ldsc_data)
  erg1_p1 = erg1_p1[is.element(snpid,SNPlist$SNP)]
  matched = match(erg1_p1$snpid,SNPlist$SNP)
  stopifnot(sum(is.na(matched))==0)
  SNPlist = SNPlist[matched,]
  stopifnot(SNPlist$SNP == erg1_p1$snpid)
  message("After filtering for match in LDSC SNPlist, ",dim(erg1_p1)[1]," SNPs remain ...")
  
  # Step 6: check allele (same as in SNPList, turn to effect allele == minor allele == A1)
  filt1 = erg1_p1$ALT == SNPlist$A1 & erg1_p1$REF == SNPlist$A2
  filt2 = erg1_p1$REF == SNPlist$A1 & erg1_p1$ALT == SNPlist$A2
  message("Checking matching alleles: \n   ",sum(filt1)," SNPs are a direct match (same coding) \n   ",sum(filt2)," SNPs have to be transformed (reversed coding)")
  erg1_p1 = erg1_p1[filt1 | filt2, ]
  message("After filtering for matching alleles in LDSC SNPlist, ",dim(erg1_p1)[1]," SNPs remain ...")
  
  data.table::setnames(erg1_p1,"ALT","A1")
  data.table::setnames(erg1_p1,"REF","A2")
  data.table::setnames(erg1_p1,"POS_b37","pos")
  data.table::setnames(erg1_p1,"POOLED_ALT_AF","EAF")
  data.table::setnames(erg1_p1,"CHROM","chr")
  
  # Step 7: get necessary columns
  pheno = unique(erg1_p1$pheno)
  myNames = c("chr","pos","snpid","A1","A2","EAF","beta","se","Zscore","N","Pvalue")
  colsOut<-setdiff(colnames(erg1_p1),myNames)
  erg1_p1[,get("colsOut"):=NULL]
  data.table::setcolorder(erg1_p1,myNames)
  data.table::setorder(erg1_p1,chr,pos)
  message("After filtering, allele check, and column check, ",dim(erg1_p1)[1]," SNPs remain for ",pheno,": ")
  print(erg1_p1)
  n1_3 = dim(erg1_p1)[1]
  
  # Step 8: filter data for second trait if there is one
  if(path_data_trait2 != "none"){
    message("Working on second phenotype ....")
    # Step 8.5: filter for SNPs in LDSC SNPlist
    erg1_p2[,snpid := rsID]
    SNPlist<-data.table::fread(path_ldsc_data)
    erg1_p2 = erg1_p2[is.element(snpid,SNPlist$SNP)]
    matched = match(erg1_p2$snpid,SNPlist$SNP)
    stopifnot(sum(is.na(matched))==0)
    SNPlist = SNPlist[matched,]
    stopifnot(SNPlist$SNP == erg1_p2$snpid)
    message("     After filtering for match in LDSC SNPlist, ",dim(erg1_p2)[1]," SNPs remain ...")
    
    # Step 8.6: check allele (same as in SNPList, turn to effect allele == minor allele == A1)
    filt1 = erg1_p2$ALT == SNPlist$A1 & erg1_p2$REF == SNPlist$A2
    filt2 = erg1_p2$REF == SNPlist$A1 & erg1_p2$ALT == SNPlist$A2
    message("     Checking matching alleles: \n        ",sum(filt1)," SNPs are a direct match (same coding) \n        ",sum(filt2)," SNPs have to be transformed (reversed coding)")
    erg1_p2 = erg1_p2[filt1 | filt2, ]
    message("     After filtering for matching alleles in LDSC SNPlist, ",dim(erg1_p2)[1]," SNPs remain ...")
    
    data.table::setnames(erg1_p2,"ALT","A1")
    data.table::setnames(erg1_p2,"REF","A2")
    data.table::setnames(erg1_p2,"POS_b37","pos")
    data.table::setnames(erg1_p2,"POOLED_ALT_AF","EAF")
    data.table::setnames(erg1_p2,"CHROM","chr")
    
    # Step 8.7: get necessary columns
    pheno = unique(erg1_p2$pheno)
    myNames = c("chr","pos","snpid","A1","A2","EAF","beta","se","Zscore","N","Pvalue")
    colsOut<-setdiff(colnames(erg1_p2),myNames)
    erg1_p2[,get("colsOut"):=NULL]
    data.table::setcolorder(erg1_p2,myNames)
    data.table::setorder(erg1_p2,chr,pos)
    message("     After filtering, allele check, and column check, ",dim(erg1_p2)[1]," SNPs remain for ",pheno,": ")
    print(erg1_p2)
    n2_3 = dim(erg1_p2)[1]
    
    # Step 8.8: get row correlation
    dum1 = data.table::copy(erg1_p1)
    dum2 = data.table::copy(erg1_p2)
    dum1 = dum1[is.element(snpid,dum2$snpid),]
    dum2 = dum2[is.element(snpid,dum1$snpid),]
    matched = match(dum1$snpid,dum2$snpid)
    dum2 = dum2[matched,]
    cor2 = cor.test(x = dum1$beta, y= dum2$beta,method = cortest_method)
    # filt = dum1$Pvalue<=0.01 | dum2$Pvalue<=0.01
    # dum1 = dum1[filt,]
    # dum2 = dum2[filt,]
    # cor2_sig = cor.test(x = dum1$beta, y= dum2$beta,method = cortest_method)
  }else{n2_3 = 0}
  
  # Step 9: save
  if(dir.exists(out_dir)==F) dir.create(out_dir) 
  out_fn1 = paste0(out_dir, "LDSC_input_",trait1_name,".txt")
  erg1_p1[,Pvalue := as.numeric(Pvalue)]
  message("Saving SNP data to ",out_fn1)
  data.table::fwrite(erg1_p1, file = out_fn1, sep = "\t")
  R.utils::gzip(out_fn1, destname = paste0(out_fn1, ".gz"))
  
  if(path_data_trait2 != "none"){
    message("Working on second phenotype ....")
    out_fn2 = paste0(out_dir, "LDSC_input_",trait2_name,".txt")
    erg1_p2[,Pvalue := as.numeric(Pvalue)]
    message("     Saving SNP data to ",out_fn2)
    data.table::fwrite(erg1_p2, file = out_fn2, sep = "\t")
    R.utils::gzip(out_fn2, destname = paste0(out_fn2, ".gz"))
  }
  
  
  
  # Step 10: report consumed time & return something
  message("\nTOTAL TIME TO CREATE MUNGED SUMSTATS: \n     " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
  
  res = data.table::data.table(pheno1 = trait1_name,
                               nSamples_pheno1 = nSamples_pheno1,
                               nSNPs_pheno1_total = n1_1,
                               nSNPs_pheno1_filt1 = n1_2,
                               nSNPs_pheno1_filt2 = n1_3,
                               pheno2 = trait2_name,
                               nSamples_pheno2 = nSamples_pheno2,
                               nSNPs_pheno2_total = n2_1,
                               nSNPs_pheno2_filt1 = n2_2,
                               nSNPs_pheno2_filt2 = n2_3,
                               time_mins = round(difftime(Sys.time(),time0,units = "mins"),3))
  if(path_data_trait2 != "none"){
    res[,cor_filt1_cor := cor1$estimate]
    res[,cor_filt1_pval := cor1$p.value]
    res[,cor_filt1_sig_cor := cor1_sig$estimate]
    res[,cor_filt1_sig_pval := cor1_sig$p.value]
    res[,cor_filt2_cor := cor2$estimate]
    res[,cor_filt2_pval := cor2$p.value]
    # res[,cor_filt2_sig_cor := cor2_sig$estimate]
    # res[,cor_filt2_sig_pval := cor2_sig$p.value]
  }
  res
  return(res)
  
}
