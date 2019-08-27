
[![](https://travis-ci.org/n-mounier/bGWAS.svg?branch=master)](https://travis-ci.org/n-mounier/bGWAS)
[![](https://img.shields.io/badge/version-1.0.0-informational.svg)](https://github.com/n-mounier/bGWAS)
[![](https://img.shields.io/badge/lifecycle-maturing-9cf.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/github/last-commit/n-mounier/bGWAS.svg)](https://github.com/n-mounier/bGWAS/commits/master)
[![](https://img.shields.io/badge/license-GPL--2.0-lightgrey.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# bGWAS

<!--- <img src="doc/Figures/logo.svg" align="right" height=140/> --->

:information\_source: `bGWAS` has been updated to version 1.0.0.  
Check the [NEWS](News.md) to learn more about what has been modified\!

:warning: If you downloaded the Z-Matrix files before 20/08/2019, they
are now obsolete and you will not be able to use them with the newest
version of the package.  
Note: some Prior GWASs have been removed, you can find more details
[here](docs/ZMatrices.md).

## Overview

bGWAS is an R-package to perform a Bayesian GWAS, using summary
statistics as input. Briefly, it compares the observed Z-scores from a
conventional GWAS to prior effects. These prior Z-scores are directly
calculated from publicly available GWASs (currently, a set of 38
studies, last update 20-08-2019 - hereinafter referred to as “prior
GWASs”). Only prior GWASs having a significant influence on the
conventional GWAS (identified using a multivariate Mendelian
Randomization (MR) approach) are used to calculate the prior effects.
Causal effects are estimated masking the focal chromosome to ensure
independence.  
Observed and prior effects are compared using Bayes Factors.
Significance is assessed by calculating the probability of observing a
value larger than the observed BF (P-value) given the prior distribution
by decomposing the analytical form of the BFs and using an approximation
for most BFs to make the computation faster. Prior, posterior and direct
effects, alongside BFs and p-values are returned. Note that prior,
posterior and direct effects are estimated on the Z-score scale, but are
automatically rescaled to beta scale if possible.

The functions available are:

  - **`bGWAS()`**  
    main function that calculates prior effects from prior GWASs,
    compares them to observed Z-scores and returns an object of class
    *bGWAS*

  - **`list_priorGWASs()`** directly returns information about the prior
    GWASs that can be used to calculate prior Z-scores

  - **`select_priorGWASs()`**  
    allows a quick selection of prior GWASs (to include/exclude specific
    studies when calculating prior Z-scores)

  - **`extract_results_bGWAS()`**  
    returns results (prior, posterior and direct estimate /
    standard-error + p-value from BF for SNPs) from an object of class
    *bGWAS*

  - **`manhattan_plot_bGWAS()`**  
    creates a Manhattan Plot from an object of class *bGWAS*

  - **`extract_MRcoeffs_bGWAS()`**  
    returns multivariate MR coefficients (1 estimate using all
    chromosomes + 22 estimates with 1 chromosome masked) from an object
    of class *bGWAS*

  - **`coefficients_plot_bGWAS()`**  
    creates a Coefficients Plot (causal effect of prior GWASs) from an
    object of class *bGWAS*

  - **`heatmap_bGWAS()`**  
    creates a heatmap to represent, for each significant SNP, the
    contribution of each prior GWAS used to create the prior effect from
    an object of class *bGWAS*

All the functions available and more details about their use can be
found in the [manual]().

## Installation

You can install the current version of `bGWAS` with:

``` r
# Directly install the package from github
# install.packages("remotes")
remotes::install_github("n-mounier/bGWAS")
library(bGWAS)
```

<!--- Note: using remotes instead of devtools leads to re-build the package
and apparently, it may be a problem with R 3.4 and macOS, 
see https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option/43943631#43943631 --->

## Usage

To run the analysis with `bGWAS` two inputs are needed:

#### 1\. The *GWAS* results to be tested

Can be a regular (space/tab/comma-separated) file or a gzipped file
(.gz) or a `data.frame`. Must contain the following columns, which can
have alternative names.  
SNP-identifier: `rs` or `rsid`, `snp`, `snpid`, `rnpid`  
Alternate allele: `a1` or `alt`, `alts`  
Reference allele: `a2` or `a0`, `ref`  
Z-statistics: `z` or `Z`, `zscore`  
If the Z-statistics is not present, it will be automatically calculated
from effect size and standard error, in which case the following columns
should be provided: Effect-size: `b` or `beta`, `beta1`  
Standard error: `se` or `std`  
If you want the prior/posterior/corrected effects to be rescaled, please
provide effect sizes and standard errors instead of (or in addition to)
Z-statistics.

#### 2\. Prior *GWASs*

Matrix files, containing Z-scores for all prior GWASs should be
downloaded separately and stored in `~/ZMatrices` or in the folder
specified when launching the analysis. These files contains the Z-scores
for all prior GWASs (before and after imputation) :  
*ZMatrix\_MR.csv.gz*: Z-scores (strong instruments only) used for
multivariate MR,  
*ZMatrix\_Full.csv.gz*: Z-scores (all SNPs) used to calculate the prior
Z-scores,  
*AvailableStudies.tsv*: A file containing information about the prior
GWASs available.

**Problem with switch drive, download link will be added soon\!**

<!---
On UNIX/MACOSX, from a terminal:    
``` bash
wget --no-check-certificate ##https://drive.switch.ch/index.php/s/t3xepllvzobVTCh/download## -O ZMatrices.tar.gz
tar xzvf ZMatrices.tar.gz
``` 
On WINDOWS:   --->

*If you want to use your own set of prior GWASs, please have a look
[here](docs/ZMatrices.md) to see how you can modify the files. * *We
focused on including prior GWASs that do not come from UKBB, assuming
that the focal phenotype is more likely to be obtained from UKBB. Sample
overlap between the focal phenotype and the Prior GWASs is not accounted
for by our method, so we did not include any UKBB results in the prior
GWASs. *

### Study Selection

Before running your analysis, you can select the prior GWASs you want to
include. You can use the function **`list_prioGWASs()`** to get some
information about the prior GWASs available.  
You should remove traits that by definition are not independent from
your trait of interest. For example, before analysing BMI results, make
sure to exclude “Height” from the prior GWASs used. You can use the
function **`select_priorGWASs()`** to automatically exclude/include some
traits or some files.  
You should also check for sample overlap, and remove prior GWASs that
come from the same consortium as your data. If there are individuals in
common between your conventional GWAS and prior GWASs, it might induce
some bias.

``` r
# Obtain the list of prior GWASs
AllStudies = list_priorGWASs()
# Select only the ones for specific traits
# select_priorGWASs will return the IDs of the files that are kept
MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body Mass Index", "Smoking"))
# Match these IDs against the ones in the list of prior GWASs 
AllStudies[AllStudies$ID %in% MyStudies, c("Name", "Trait", "Reference")]
```

    ##                                  Name           Trait
    ## 1             Body Mass Index (GIANT) Body Mass Index
    ## 23                Heart Rate (HRgene)      Heart Rate
    ## 35 Smoking - cigarettes per day (TAG)         Smoking
    ## 36        Smoking - ever smoked (TAG)         Smoking
    ## 37      Smoking - former smoker (TAG)         Smoking
    ## 38       Smoking - age of onset (TAG)         Smoking
    ##                                       Reference
    ## 1  https://www.ncbi.nlm.nih.gov/pubmed/25673413
    ## 23 https://www.ncbi.nlm.nih.gov/pubmed/23583979
    ## 35 https://www.ncbi.nlm.nih.gov/pubmed/20418890
    ## 36 https://www.ncbi.nlm.nih.gov/pubmed/20418890
    ## 37 https://www.ncbi.nlm.nih.gov/pubmed/20418890
    ## 38 https://www.ncbi.nlm.nih.gov/pubmed/20418890

### Analysis

  - **Example A**

<!-- end list -->

``` r
# Using a small GWAS (100,000 SNPs, Timmers et al data - stored as a data.frame)
# Please, not that this example is only for illustration, the method is designed
# to be used genome-wide, and using such a low number of SNPs can not yield
# interpretable results.
data("SmallGWAS_Timmers2019")
MyStudies = select_priorGWASs(include_traits=c("Blood Pressure", "Education"),  
include_files=c("cardiogram_gwas_results.txt", "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz"))
# 6 Prior GWASs used
list_priorGWASs(MyStudies)[, c("Name", "Trait", "Reference")]
```

    ##                                   Name                   Trait
    ## 1              Body Mass Index (GIANT)         Body Mass Index
    ## 2 Coronary Artery Disease (CARDIoGRAM) Coronary Artery Disease
    ## 3           Years of Schooling (SSGAC)               Education
    ## 4       Systolic Blood Pressure (ICBP)          Blood Pressure
    ## 5      Diastolic Blood Pressure (ICBP)          Blood Pressure
    ## 6           College Completion (SSGAC)               Education
    ##                                      Reference
    ## 1 https://www.ncbi.nlm.nih.gov/pubmed/25673413
    ## 2 https://www.ncbi.nlm.nih.gov/pubmed/21378990
    ## 3 https://www.ncbi.nlm.nih.gov/pubmed/27225129
    ## 4 https://www.ncbi.nlm.nih.gov/pubmed/21909115
    ## 5 https://www.ncbi.nlm.nih.gov/pubmed/21909115
    ## 6 https://www.ncbi.nlm.nih.gov/pubmed/23722424

``` r
A = bGWAS(name="Test_UsingSmallDataFrame",
          GWAS = SmallGWAS_Timmers2019,
          prior_studies = MyStudies,
          stepwise_threshold = 0.05)
# MR instruments will be selected using default parameters,
# MR will be performed using a threshold of 0.05 to select studies, and the default shrinkage threshold,
# A subset of prior GWASs (MyStudies) will be used to create the prior,
# Significant SNPs will be identified using default parameters (p<5e-8) and distance-pruned (500kb),
# No file will be saved.
```

  - **Example B**

<!-- end list -->

``` r
# Using a GWAS from our list our prior GWASs
# Using all other (37) GWASs to built the prior
MyGWAS = 3
list_priorGWASs(MyGWAS)[, c("Name", "Trait", "Reference")]
```

    ##                                   Name                   Trait
    ## 1 Coronary Artery Disease (CARDIoGRAM) Coronary Artery Disease
    ##                                      Reference
    ## 1 https://www.ncbi.nlm.nih.gov/pubmed/21378990

``` r
B = bGWAS(name = "Test_UsingGWASfromPriorGWASs",
         GWAS = MyGWAS)
# MR instruments will be selected using default parameters,
# MR will be performed using default parameters (stepwise selection / shrinkage threshold)
# All Prior GWASs except the one use as "GWAS" will be used to create the prior,
# Significant SNPs will be identified using default parameters (p<5e-8) and distance-pruned (500kb)
# No file will be saved.
```

### Results

**`bGWAS()`** returns an object of class “bGWAS” than can be handled in
`R`.

``` r
class(A)
```

    ## [1] "bGWAS"

``` r
print(A)
```

    ## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ 
    ##  
    ##   Analysis : "Test_UsingSmallDataFrame" 
    ## bGWAS performed on 291,583 SNPs 
    ##  
    ## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ 
    ##  
    ## 1 study used to build the prior : 
    ##                                              study    estimate  std_error
    ##  All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz -0.06433771 0.03200732
    ## 
    ## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ 
    ##  
    ## 5 significant SNPs identified : 
    ##  rs17486195, rs1160985, rs1333049, rs12993295, rs78801493
    ## 
    ## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

Functions to extract results from an object of class *bGWAS*:

``` r
hits = extract_results_bGWAS(A, "significant")
colnames(hits)
```

    ##  [1] "chrm_UK10K"               "pos_UK10K"               
    ##  [3] "rsid"                     "alt"                     
    ##  [5] "ref"                      "beta"                    
    ##  [7] "se"                       "z_obs"                   
    ##  [9] "mu_prior_estimate"        "mu_prior_std_error"      
    ## [11] "mu_posterior_estimate"    "mu_posterior_std_error"  
    ## [13] "mu_direct_estimate"       "mu_direct_std_error"     
    ## [15] "beta_prior_estimate"      "beta_prior_std_error"    
    ## [17] "beta_posterior_estimate"  "beta_posterior_std_error"
    ## [19] "beta_direct_estimate"     "beta_direct_std_error"   
    ## [21] "BF"                       "BF_p"

``` r
nrow(hits)
```

    ## [1] 5

``` r
all_results = extract_results_bGWAS(A, "all")
nrow(all_results)
```

    ## [1] 291583

``` r
MR_coefficients = extract_MRcoeffs_bGWAS(A)
colnames(MR_coefficients)
```

    ##  [1] "name"             "study"            "estimate"        
    ##  [4] "std_error"        "Tstat"            "P"               
    ##  [7] "chrm1_estimate"   "chrm1_std_error"  "chrm1_P"         
    ## [10] "chrm2_estimate"   "chrm2_std_error"  "chrm2_P"         
    ## [13] "chrm3_estimate"   "chrm3_std_error"  "chrm3_P"         
    ## [16] "chrm4_estimate"   "chrm4_std_error"  "chrm4_P"         
    ## [19] "chrm5_estimate"   "chrm5_std_error"  "chrm5_P"         
    ## [22] "chrm6_estimate"   "chrm6_std_error"  "chrm6_P"         
    ## [25] "chrm7_estimate"   "chrm7_std_error"  "chrm7_P"         
    ## [28] "chrm8_estimate"   "chrm8_std_error"  "chrm8_P"         
    ## [31] "chrm9_estimate"   "chrm9_std_error"  "chrm9_P"         
    ## [34] "chrm10_estimate"  "chrm10_std_error" "chrm10_P"        
    ## [37] "chrm11_estimate"  "chrm11_std_error" "chrm11_P"        
    ## [40] "chrm12_estimate"  "chrm12_std_error" "chrm12_P"        
    ## [43] "chrm13_estimate"  "chrm13_std_error" "chrm13_P"        
    ## [46] "chrm14_estimate"  "chrm14_std_error" "chrm14_P"        
    ## [49] "chrm15_estimate"  "chrm15_std_error" "chrm15_P"        
    ## [52] "chrm16_estimate"  "chrm16_std_error" "chrm16_P"        
    ## [55] "chrm17_estimate"  "chrm17_std_error" "chrm17_P"        
    ## [58] "chrm18_estimate"  "chrm18_std_error" "chrm18_P"        
    ## [61] "chrm19_estimate"  "chrm19_std_error" "chrm19_P"        
    ## [64] "chrm20_estimate"  "chrm20_std_error" "chrm20_P"        
    ## [67] "chrm21_estimate"  "chrm21_std_error" "chrm21_P"        
    ## [70] "chrm22_estimate"  "chrm22_std_error" "chrm22_P"

``` r
MR_coefficients[1, 1:8]
```

    ##                      name
    ## 1 Body Mass Index (GIANT)
    ##                                               study    estimate  std_error
    ## 1 All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz -0.06433771 0.03200732
    ##       Tstat          P chrm1_estimate chrm1_std_error
    ## 1 -2.010094 0.04900041    -0.07819836      0.03233513

Functions for graphic representations:

``` r
manhattan_plot_bGWAS(A)
```

<img src="doc/Figures/README-PlotsA-1.png" width="100%" />

``` r
coefficients_plot_bGWAS(A) 
```

<img src="doc/Figures/README-PlotsA-2.png" width="100%" />

##### Aditionnaly, if `save_files=T`, several files are created in the folder `./name/` :

  - **name.log** - log file  
  - **PriorGWASs.tsv** - contains information about all prior GWASs
    (general info + status (used/removed) + univariate/multivariate MR
    estimates)  
  - **CoefficientsByChromosome.csv** - contains the multivariate MR
    estimates when masking the focal chromosome (22 coefficients for
    each prior GWASs used for prior estimation)  
  - **PriorBFp.csv** - contains BF and p-values, prior, posterior and
    direct effects estimates for all SNPs  
  - **SignificantSNPs.csv** - contains BF and p-values, prior, posterior
    and direct effects estimates for a subset of significant SNPs
    (identified according to specified parameters)

A detailed description of these files can be found
[here](doc/OutputFiles.md).

## Runtime

Analysis using all the 38 prior GWASs available, for a conventional GWAS
containing \~7M SNPs in common with the prior studies \~ 25 minutes (see
complete Lifespan Analysis [here]()).

Analysis using 6 prior GWASs, for a conventional GWAS containing \~
400,000 SNPs in common with prior studies (see example A) \~ 2 minutes.

## Contact

<mounier.ninon@gmail.com>
