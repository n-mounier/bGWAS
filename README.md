
# bGWAS
[//]:========================================

:information_source: package under development :information_source:    
:warning: if you downloaded the Z-Matrix files before 25/02/2018, they are obsolete, you need to delete the old ones and download the new ones!

## Overview
[//]:-------------------------------

bGWAS is an R-package to perform a Bayesian GWAS, using summary statistics as input. Briefly, it compares the observed Z-scores from a conventional GWAS to prior Z-scores. These prior Z-scores can be provided by the user of directly calculated from publicly available GWASs (currently, a set of 58 studies, last update dd-mm-yyyy - hereinafter referred to as "prior GWASs"). In this case, only prior GWASs having a significant influence on the conventional GWAS (identified using a multivariate Mendelian Randomization (MR) approach) are used to calculate the prior Z-scores. Causal effect are estimated masking the focal chromosome to ensure independence.          
Observed and prior Z-scores are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   

The main functions are:   
-   `bGWAS()` -  core function that calculates prior Z-scores from prior GWASs, compares them to observed Z-scores and returns an object of class "bGWAS"    
<!--- returns an object of class `bGWAS-class`. See the vignette: vignette('vcf_data')
THIS USE RISK FACTORS TO CREATE THE PRIOR---> 
-   `list_priorGWASs()` directly returns information about the prior GWASs that can be used to calculate prior Z-scores   
-   `select_priorGWASs()` allows a quick selection of prior GWASs (to include/exclude specific studies when calculating prior Z-scores)   
-   `manatthan_plot_bGWAS()` Create a Manhattan Plot from bGWAS results
-   `coefficients_plot_bGWAS()` Create a Coefficients Plot (causal effect of Prior GWASs) from bGWAS results    
-   `bGWAS_fromPrior()` alternative function that compares prior Z-scores provided by the user to observed Z-scores and returns an object of class "bGWAS" # NOT IMPLEMENTED YET   



## Installation
[//]:-------------------------------

If you want to use the `bGWAS()` function to calculate prior Z-scores, you should start by downloading Z-Matrix files.   

* Download Z-Matrix files :   
These files contains the Z-scores for all prior GWASs (before and after imputation) and were created ... 
Z-scores before imputation are used for multivariate MR.  
Z-scores after imputation are used to calculate the prior Z-scores. 

`wget --no-check-certificate https://drive.switch.ch/index.php/s/BpRrDXvFPbnKCM6/download -O ZMatrices.tar.gz`

`tar xzvf ZMatrices.tar.gz`

Size ~ 2.09 GB   
  

* Install R-package
``` r
# Directly install the package from github
# install.packages("devtools")
devtools::install_github("n-mounier/bGWAS")
```


## Usage
[//]:-------------------------------

To run the analysis with `bGWAS` two inputs are needed:

1. The *GWAS* results to be tested   
Can be a regular (space/tab/comma-separated) file or a gzipped file (.gz), must contain the following columns, which can have alternative names.  
SNP-identifier:  `rs` or `rsid`, `snp`, `snpid`, `rnpid`    
Alternate allele:  `a1` or `alt`, `alts`    
Reference allele: `a2` or `a0`, `ref`    
Z-statistics: `z` or `Z`, `zscore`      
If the Z-statistics is not present, it can be calculated from effect size and standard error, in which case the following columns should be provided:
Effect-size: `b` or `beta`, `beta1`    
Standard error:  `se` or `std`     

2. Prior *GWASs* (downloaded above) - not needed for the function `bGWAS_fromPrior()`  
Matrix files, containing Z-scores for all prior GWASs should be downloaded separately and stored in `~/ZMatrices` or in the folder specified when launching the analysis.   

* Currently, this is not possible for users to directly add their own GWASs to this set of prior GWASs. Feel free to contact us for more information. *   
 
<!---  Format?
Can I add one more?--->


### Study Selection

Before running your analysis, you can select the prior GWASs you want to include. You can use the function `list_prioGWASs()` to get some information about the prior GWASs available.   
You should remove traits that by definition are not independent from your trait. For example, before analysing BMI results, make sure to exclude "Height" from prior GWASs. You can use the function `select_priorGWASs()` to automatically exclude/include some traits or some files.   
+ also check consortium (no function for selection) but if same individuals in your GWAS + prior GWASs not ok

``` r
AllStudies = list_priorGWASs()
MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body mass index", "Smoking"))
AllStudies[AllStudies$ID %in% MyStudies, ]
```

### Analysis
``` r
## Example A
# Using a GWAS from our list our prior GWASs
# Using all other (57) GWASs to built the prior
MyGWAS = 5
list_priorGWASs(MyGWAS)
A = bGWAS(name = "Test_UsingGWASfromPriorGWASs",
         GWAS = MyGWAS
         verbose=T)
         

## Example B
# Using a small GWAS (400,000 SNPs, Pilling et al data - file)
# Using only specific traits / files (resulting in 9 GWASs included)
MyGWAS = system.file("Data/SmallGWAS_Pilling2017.csv", package="bGWAS")
MyStudies = select_priorGWASs(include_traits=c("Type 2 diabetes", "Smoking"),
                          include_files=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))
list_files(MyStudies)
 
B = bGWAS(name = "Test_UsingSmallGWAS",
         GWAS = MyGWAS,
         prior_studies=MyStudies,
         verbose=T) 
print(B)
manatthan_plot_bGWAS(B)
```
    
``` r         
## Example C
# Using a small GWAS (400,000 SNPs, Pilling et al data - data.frame)
# Using only specific traits / files (resulting in 9 GWASs included)
data("SmallGWAS_Pilling2017")
MyStudies = select_priorGWASs(include_traits=c("Type 2 diabetes", "Smoking"),
                          include_files=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))

C = bGWAS(name="Test_UsingSmallDataFrame",
         GWAS = SmallGWAS_Pilling2017,
         prior_studies=MyStudies,
         verbose=T,
         save_files=T)
print(C)
         
         
# Note that B and C are using the same data (stored differently) and give the same results.

         
```

## Results
`bGWAS()` returns an object of class "bGWAS" than can be handled in `R`. \cr

```r
print(myObj)
```
-> insert figure here

Functions for graphic representations:   
`manatthan_plot_bGWAS(myObj)`   
-> insert figure here   
`coefficients_plot_bGWAS(myObj)`   
-> insert figure here  

extract significant SNPs   
extract all results   
extract significant studies   

Aditionnaly, if `save_files=T`, several files are created...   
... in your working directory :    
-   "`name`.log" - log file    
... in the folder `./name/` :   
-   "PriorGWASs.tsv" - contains Prior GWASs information (general info + status (used/removed) + MR coefficients)   
-   "CoefficientsByChromosome.csv" - contains the MR estimates when masking the focal chromosome (22 coefficients / study)    
-   "PriorBFp.csv" - contains BF and p-values estimated for all SNPs    
-   "SignificantSNPs.csv" - contains BF and p-values estimated for a subset of SNPs    


Further description of the files?   


## Runtime
[//]:-------------------------------

Analysis using all the 58 prior GWASs available, for a conventional GWAS containing ~7M SNPs in common with the prior studies ~ 145 minutes.

Analysis using 9 prior GWASs, for a conventional GWAS containing 400,000 SNPs in common with prior studies (see example B) ~ 8 minutes


## Improvements to be implemented
[//]:-------------------------------

document results files    

bGWASfromPrior()    

selection from consortium (maybe not possible - too complicated, too many cohorts in each consortium...)    
use of a subset of SNPs   
use of re-imputed studies for prior   
use of additional studies for prior   

optimize null-BF calculation




## Contact
<mounier.ninon@gmail.com>



