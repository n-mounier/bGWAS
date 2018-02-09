
# bGWAS
[//]:========================================


## Overview
[//]:-------------------------------

bGWAS is an R-package to perform a Bayesian GWAS, using summary statistics as input. Briefly, it compares the observed Z-score from a conventional GWAS to a prior z-score, calculated from prior, publicly available GWASs (currently, a set of 58 studies, last update dd-mm-yyyy).   
Observed Z-scores and priors are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   

The main functions are:

-   `bGWAS()` is core function that will return ...
-   `availableStudies()` will directly return the available GWASs that are available as priors
-   `selectStudies()` allows a quick selection of GWASs as priors
-   `bGWASfromPrior()` does ... # NOT IMPLEMENTED YET


## Installation
[//]:-------------------------------


* Download Z-Matrix (GWASs) files : 

`wget --no-check-certificate https://drive.switch.ch/index.php/s/pxZsWY88RSDsO8K/download -O ZMatrices.tar.gz`    

`tar xzvf ZMatrices.tar.gz`

Size ~ 2.09 GB   

Origin?  

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

SNP-identifier: `SNPID` or `rs`, `rsid`, `snp`, `snpid`, `rnpid`    
Alternate allele: `ALT` or `a1`, `alt`, `alts`    
Reference allele: `REF` or `a2`, `a0`, `ref`    
Z-statistics: `Z` or `z`, `Z`, `zscore`    
If the Z-statistics is not present, it can be calculated from effect size and standard error, in which case the following columns should be provided:
Effect-size: `BETA` or `b, `beta`, `beta1`    
Standard error: `SE` or `se`, `std   

2. Prior *GWASs* (downloaded above)
Matrix files, containing Z-scores for all prior GWASs should be downloaded separately and stored in `~/ZMatrices` or in the folder specified when launching the analysis.

Format?
Can I add one more?

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
# Using a small GWAS (400,000 SNPs, Pilling et al data)
# Using only specific traits / files (resulting in 9 GWASs included)
MyGWAS = system.file("Data/SmallGWAS_Pilling2017.csv", package="bGWAS")
MyStudies = selectStudies(includeTraits=c("Type 2 diabetes", "Smoking"),    
                          includeFiles=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))
listFiles(MyStudies)
 
B = bGWAS(Name = "Test_UsingSmallGWAS",
         GWAS = MyGWAS,
         PriorStudies=MyStudies,
         verbose=T)    
             
```


## Runtime
[//]:-------------------------------

Analysis using all the 58 prior GWASs available, for a conventional GWAS containing ~7M SNPs in common with the prior studies ~ 145 minutes.

Analysis using 9 prior GWASs, for a conventional GWAS containing 400,000 SNPs in common with prior studies (see example B above) ~ 8 minutes


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
