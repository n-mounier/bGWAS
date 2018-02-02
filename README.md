
# bGWAS
[//]:========================================


## Overview
[//]:-------------------------------

bGWAS is an R-package to perform a bayesian GWAS, using summary statistics. It compares the observed Z-score from a conventional GWAS to a prior z-score, calculated from prior GWASs publicly available (currently, a set of 58 studies, to be updated).   
Observed Z-scores and priors are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   


-   `bGWAS()` ...
-   `availableStudies()` ...
-   `selectStudies()` ...
-   `bGWASfromPrior()`  # NOT IMPLEMENTED YET


## Installation
[//]:-------------------------------


* Download Z-Matrix files : 
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
## Example A
# Using a GWAS from our list our prior GWASs
# Using all other (57) GWASs to built the prior
MyGWAS = 1
listFiles(MyGWAS)
A = bGWAS(Name = "Test_UsingGWASfromList",
         GWAS = MyGWAS
         verbose=T)
         

## Example B
# Using a small GWAS (400,00 SNPs, Pilling et al data)
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

Analysis using all the 58 prior GWASs available, for a conventionnal GWAS containing ~7M SNPs in common with the prior studies ~ 145 minutes.

Analysis using 9 prior GWASs, for a conventionnal GWAS containtin 400,000 SNPs in commons with prior studies (example B) ~ 8 minutes


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
