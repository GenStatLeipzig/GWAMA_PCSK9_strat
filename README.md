# Genome-wide association meta-analysis (GWAMA) of PCSK9 levels

**Contributor: Janne Pott**

**Last Updated: 22/01/2024**

Supporting code for the following paper:

-   Pott J, Kheirkhah A, Gadin JR, Kleber ME, Delgado GE, Kirsten H, et al. Sex and statin-related genetic associations at the PCSK9 gene locus – results of genome-wide association meta-analysis. *Biology of Sex Differences*. 2023 (under review) 

This project is an extension to our previous publications:

-   Pott J, Burkhardt R, Beutner F, Horn K, Teren A, Kirsten H, et al. Genome-wide meta-analysis identifies novel loci of plaque burden in carotid artery. *Atherosclerosis*. 2017;259:32--40. [DOI](https://doi.org/10.1016/j.atherosclerosis.2017.02.018).
-   Pott J, Gadin J, Theusch E, Kleber ME, Delgado GE, Kirsten H, et al. Meta-GWAS of PCSK9 levels detects two novel loci at APOB and TM6SF2. *Hum Mol Genet* 2021. [DOI](https://doi.org/10.1093/hmg/ddab279).

We are providing the main scripts used in the GWAMA of PCSK9 levels in four European Cohorts (LIFE-Adult, LIFE-Heart, LURIC, TwinGene, KORA-F3, and GCKD), stratified by sex and statin treatment, to empower other researchers to reproduce our results, starting from the summary statistics. Data can be found on zenodo (LINK *tba*)

## Source File

You will need to customize a source file, indicating

-   R library and R packages: all scripts were run under [R Version 4.x](https://cran.r-project.org/). All necessary packages are listed in the source file. Additional function not published in any R package are listed in the directory 'helperFunctions'.

-   Path to tools other than R:

    -   [PLINK 2.00](https://www.cog-genomics.org/plink/2.0/)
    -   [GCTA 1.94.1](https://yanglab.westlake.edu.cn/software/gcta/#Download)

-   Path to downloaded data sets used throughout the analyses:

    -   [1000 Genomes Phase 3 EUR data](https://www.internationalgenome.org/data-portal/data-collection/phase-3)
    -   [GTEx v8 data](https://gtexportal.org/home/protectedDataAccess)
    -   [GTEx v8 data - sex stratified](https://gtexportal.org/home/protectedDataAccess)
    -   [Summary statistics for lipid - sex-stratified](http://csg.sph.umich.edu/willer/public/glgc-lipids2021/), publication: [Kanoni et al. (2022)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02837-1)
    -   [Summary statistics for lipids - sex-combined](http://csg.sph.umich.edu/willer/public/glgc-lipids2021/), publication: [Graham et al. (2021)](https://www.nature.com/articles/s41586-021-04064-3)
    -   [Summary statistics for coronary artery disease](https://data.mendeley.com/datasets/gbbsrpx6bs/1), publication: [van der Harst et al. (2018)](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.117.312086)
    -   [Summary statistics for sleep duration](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007561/), publication: [Dashti et al. (2019)](https://pubmed.ncbi.nlm.nih.gov/30846698/)

## Scripts

R scripts staring with *0x*:

1.  Get summary statistics as uploaded to zenodo (**documentary**, uses GWAS pipeline output, you will not need to rerun this when you downloaded the zenodo data)
2.  Define associated loci
3.  Fine-mapping of *PCSK9* locus (GCTA conditional joint analyses)
4.  Interaction Tests
5.  Co-localization
    1.  Preparation of data
    2.  Run within PCSK9 data
    3.  Run against eQTLs
    4.  Run against other GWAS traits
6.  Mendelian Randomization
    1.  Preparation of UKBB data (GWAS for stratified LDLC)
    2.  PCSK9 on LDL-C using UKBB data for LDLC summary statistics
    3.  PCSK9 on LDL-C using GLGC data
7.  Look-up of sex-biased gene expression and sex-biased eQTLs

## Main Figures

R scripts staring with *MF*:

1.  Heatmap of independent SNPs at *PCSK9* gene loci
2.  Interaction scatter plot of SNP estimates
3.  Forest plot of causal estimates per subgroup
4.  Interaction scatter plot of causal estimates

## Main Tables

R scripts staring with *MT*:

1.  Summary of independent SNPs in the PCSK9 GWAS
2.  Results of the SNP interaction analysis
3.  Results of the MR analysis
4.  Results of the MR interaction analysis

## Supplemental Tables

R script staring with *ST*:

1.  Description of Studies (as received from participating studies) **--\> not included here**
2.  Sample Sizes, SNP Numbers, genomic inflation factor $\lambda$, and LDSC heritability results per phenotype
3.  Overview of associated loci
4.  Annotation of associated SNPs
5.  GCTA COJO results
6.  Interaction results
7.  Colocalization results
8.  Mendelian Randomozation results
