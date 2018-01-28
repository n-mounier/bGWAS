
# bGWAS
[//]:========================================


## Overview
[//]:-------------------------------

bGWAS is ann R-package to perform a bayesian GWAS, using summary statistics. It compares the observed Z-score from a conventional GWAS to a prior z-score, calculated from prior GWASs publicly available (currently, a set of 58 studies, to be updated).   
Observed Z-scores and priors are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   


-   `bGWAS()` ...
-   `availableStudies()` ...
-   `selectStudies()` ...
-   `bGWASfromPrior()` ...


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
- GWAS file
SNPID : 
ALT :
REF :
Z or BETA & SE :
- ZMatrices





## Runtime
[//]:-------------------------------



## Improvements to be implemented
[//]:-------------------------------





## Contact
<mounier.ninon@gmail.com>
