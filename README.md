
# bGWAS
[//]:========================================


## Overview
[//]:-------------------------------

bGWAS is ann R-package to perform a bayesian GWAS, using summary statistics. It compares the observed Z-score from a conventional GWAS to a prior z-score, calculated from prior GWASs publicly available (currently, a set of 58 studies, to be updated).   
Observed Z-scores and priors are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   


-   `bGWAS()` ...
-   `availableStudies()` ...
-   `selectStudies()` ...
-   `bGWASfromPrior()`  # NOT IMPLEMENTED YET


## Installation
[//]:-------------------------------

* Install R-package
``` r
# Directly install the package from github
# install.packages("devtools")
devtools::install_github("n-mounier/bGWAS")
```

* Download Z-Matrix files   
`wget --no-check-certificate https://drive.switch.ch/index.php/s/pxZsWY88RSDsO8K/download -O ZMatrices.tar.gz`    
`tar xzvf ZMatrices.tar.gz`

## Usage
[//]:-------------------------------

To run an anlysis :

- GWAS file, can be a regular (space/tab/comma-separated) file or a gzipped file (.gz), must contain the columns   
SNPID : rs, rsid, snp, snpid, rnpid    
ALT : a1, alt, alts    
REF : a2, a0, ref    
Z : z, Z, zscore    
If Z is not present, can be calculated from effect size and standard error  
BETA : b, beta, beta1    
SE : se, std   

- ZMatrices    
Matrix files, containing Z-scores for all prior GWASs should be downloaded separately and stored in "~/ZMatrices" or in the folder specified when launching the analysis.


### Study Selection
``` r
AllStudies = availableStudies()
listTraits()
MyStudies = selectStudies(includeTraits=c("Heart Rate", "Body mass index", "Smoking"))
AllStudies[AllStudies$ID %in% MyStudies, ]
```

### Analysis
``` r
MyGWAS = 1
listFiles(MyGWAS)
A = bGWAS(Name = "Test_UsingGWASfromList",
         GWAS = MyGWAS
         verbose=T)
```


## Runtime
[//]:-------------------------------

Analysis using all the 58 traits available, for a conventionnal GWAS containing ~7M SNPs in common with the prior studies : 145 minutes.


## Improvements to be implemented
[//]:-------------------------------

bGWASfromPrior()    
selection from consortium    
use of a subset of SNPs   
use of re-imputed studies for prior   
use of additional studies for prior   




## Contact
<mounier.ninon@gmail.com>
