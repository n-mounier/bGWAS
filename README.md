
# bGWAS
[//]:========================================

:information_source: package under development

## Overview
[//]:-------------------------------

bGWAS is an R-package to perform a Bayesian GWAS, using summary statistics as input. Briefly, it compares the observed Z-score from a conventional GWAS to a prior Z-score, calculated from publicly available GWASs (currently, a set of 58 studies, last update dd-mm-yyyy - hereinafter referred to as "prior GWASs"). Only prior GWASs having a significant influence on the conventional GWAS (identified using a multivariate Mendelian Randomization (MR) approach) are used to calculate the prior Z-scores.          
Observed and prior Z-scores are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   

The main functions are:   
-   `bGWAS()` -  core function that will return a data.frame containing significant   
<!--- returns an object of class `bGWAS-class`. See the vignette: vignette('vcf_data')  ---> 
-   `availableStudies()` will directly return the available prior GWASs that can be used to calculate prior Z-scores   
-   `selectStudies()` allows a quick selection of prior GWASs (to include/exclude specific studies when calculating prior Z-scores)   
-   `bGWASfromPrior()` compare prior Z-scores pre-calculated by the user to observed Z-scores # NOT IMPLEMENTED YET   


## Installation
[//]:-------------------------------


* Download Z-Matrix files :   
These files contains the Z-scores for all prior GWASs (before and after imputation) and were created ... 
Z-scores before imputation are used for multivariate MR.  
Z-scores after imputation are used to calculate the prior Z-scores. 

`wget --no-check-certificate https://drive.switch.ch/index.php/s/pxZsWY88RSDsO8K/download -O ZMatrices.tar.gz`    

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

2. Prior *GWASs* (downloaded above)   
Matrix files, containing Z-scores for all prior GWASs should be downloaded separately and stored in `~/ZMatrices` or in the folder specified when launching the analysis.
 
<!---  Format?
Can I add one more?--->


### Study Selection
``` r
AllStudies = availableStudies()
listTraits()
MyStudies = selectStudies(includeTraits=c("Heart Rate", "Body mass index", "Smoking"))
AllStudies[AllStudies$ID %in% MyStudies, ]
```

### Analysis
``` r
## Example A
# Using a GWAS from our list our prior GWASs
# Using all other (57) GWASs to built the prior
MyGWAS = 1
listFiles(MyGWAS)
A = bGWAS(Name = "Test_UsingGWASfromList",
         GWAS = MyGWAS
         verbose=T)
         

## Example B
# Using a small GWAS (400,000 SNPs, Pilling et al data - file)
# Using only specific traits / files (resulting in 9 GWASs included)
MyGWAS = system.file("Data/SmallGWAS_Pilling2017.csv", package="bGWAS")
MyStudies = selectStudies(includeTraits=c("Type 2 diabetes", "Smoking"),
                          includeFiles=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))
listFiles(MyStudies)
 
B = bGWAS(Name = "Test_UsingSmallGWAS",
         GWAS = MyGWAS,
         PriorStudies=MyStudies,
         verbose=T) 
    
         
## Example C
# Using a small GWAS (400,000 SNPs, Pilling et al data - data.frame)
# Using only specific traits / files (resulting in 9 GWASs included)
data("SmallGWAS_Pilling2017")
MyStudies = selectStudies(includeTraits=c("Type 2 diabetes", "Smoking"),
                          includeFiles=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))

C = bGWAS(Name="Test_UsingSmallDataFrame",
         GWAS = SmallGWAS_Pilling2017,
         PriorStudies=MyStudies,
         verbose=T,
         saveFiles=T)
         
         
# Note that B and C are using the same data (stored differently) and give the same results.

         
```


## Runtime
[//]:-------------------------------

Analysis using all the 58 prior GWASs available, for a conventional GWAS containing ~7M SNPs in common with the prior studies ~ 145 minutes.

Analysis using 9 prior GWASs, for a conventional GWAS containing 400,000 SNPs in common with prior studies (see example B) ~ 8 minutes


## Improvements to be implemented
[//]:-------------------------------

document results files    

bGWASfromPrior()    

selection from consortium    
use of a subset of SNPs   
use of re-imputed studies for prior   
use of additional studies for prior   

optimize null-BF calculation




## Contact
<mounier.ninon@gmail.com>
