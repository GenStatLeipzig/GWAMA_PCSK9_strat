# Genome-wide association meta-analysis (GWAMA) of PCSK9 levels

**Contributor: Janne Pott**

**Last Updated: 29/01/2023**

Supporting code for the following paper:

* Pott J, et al. Sex- and statin-stratified GWAMA of PCSK9 reveal novel hits in Europeans. 

This project is an extension to our previous publications:

* Pott J, Burkhardt R, Beutner F, Horn K, Teren A, Kirsten H, et al. Genome-wide meta-analysis identifies novel loci of plaque burden in carotid artery. _Atherosclerosis_. 2017;259:32–40. [DOI](https://doi.org/10.1016/j.atherosclerosis.2017.02.018).
* Pott J, Gadin J, Theusch E, Kleber ME, Delgado GE, Kirsten H, et al. Meta-GWAS of PCSK9 levels detects two novel loci at APOB and TM6SF2. _Hum Mol Genet_ 2021. [DOI](https://doi.org/10.1093/hmg/ddab279).

We are providing the main scripts used in the GWAMA of PCSK9 levels in four European Cohorts (LIFE-Adult, LIFE-Heart, LURIC, and TwinGene), stratified by sex and statin treatment, to empower other researchers to reproduce our results, starting from the summary statistics. Data can be found on zenodo (LINK tba)

## Source File

You will need to customize a source file, indicating

- R library and R packages: all scripts were run under [R Version 4.x](https://cran.r-project.org/). All necessary packages are listed in the source file. Additional function not published in any R package are listed in the directory 'helperFunctions'. 
- Path to tools other than R: 

    - [PLINK 2.00](https://www.cog-genomics.org/plink/2.0/)
    - [GCTA 1.94.1](https://yanglab.westlake.edu.cn/software/gcta/#Download)
    - [Python/Conda 3.6](https://www.anaconda.com/products/individual)
    - [LDSC 1.0.1](https://github.com/bulik/ldsc)

- Path to downloaded data sets used throughout the analyses:

    - [1000 Genomes Phase 3 EUR data](https://www.internationalgenome.org/data-portal/data-collection/phase-3)
    - [GTEx v8 data](https://gtexportal.org/home/protectedDataAccess)
    - [Sex-stratified summary statistics for lipids](http://csg.sph.umich.edu/willer/public/glgc-lipids2021/), publication: [Kanoni et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02837-1)
    - [Summary statistics for coronary artery disease](https://data.mendeley.com/datasets/gbbsrpx6bs/1), publication: [van der Harst et al.](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.117.312086)
    - [Summary statistics for total bilirubin levels](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90019001-GCST90020000/GCST90019521/), publication: [Sinnott-Armstrong et al.](https://pubmed.ncbi.nlm.nih.gov/33462484/)
    - [Summary statistics for sleep duration](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007561/), publication: [Dashti et al.](https://pubmed.ncbi.nlm.nih.gov/30846698/)
    - [Summary statistics for systolic blood pressure](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018972/), publication: [Sakaue et al.](https://pubmed.ncbi.nlm.nih.gov/34594039/)
    - [Summary statistics for pulse pressure](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018970/), publication: [Sakaue et al.](https://pubmed.ncbi.nlm.nih.gov/34594039/)
    - [Summary statistics for medication use (agents acting on the renin-angiotensin system)](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018988/), publication: [Sakaue et al.](https://pubmed.ncbi.nlm.nih.gov/34594039/)
 
## Scripts 

1. Get summary statistics as uploaded to zenodo (**documentary**, uses GWAS pipeline output, you will not need to rerun this when you downloaded the zenodo data)
2. Define associated loci 
3. Interaction Tests
4. LD Score regression
5. Fine-mapping
    1. GCTA conditional joint analyses
    2. Credible Sets
6. Co-localization
    1. Preparation of data
    2. Run within PCSK9 data
    3. Run against eQTLs
    4. Run against other GWAS traits
7. Mendelian Randomization
    1. Direction 1: PCSK9 on LDL-C
    2. Direction 2: LDL-C on PCSK9
    
## Figures

1. Interaction plot of 3-way interaction
2. Heatmap of independent SNPs at _PCSK9_ gene loci
3. Co-localization plot at _PCSK9_ gene loci (eQTLs + lipids + CAD)
