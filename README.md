
[![](https://travis-ci.org/n-mounier/bGWAS.svg?branch=master)](https://travis-ci.org/n-mounier/bGWAS)
[![](https://img.shields.io/badge/version-1.0.0-informational.svg)](https://github.com/n-mounier/bGWAS)
[![](https://img.shields.io/badge/lifecycle-maturing-9cf.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/github/last-commit/n-mounier/bGWAS.svg)](https://github.com/n-mounier/bGWAS/commits/master)
[![](https://img.shields.io/badge/license-GPL--2.0-lightgrey.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# bGWAS

<!--- <img src="doc/Figures/logo.svg" align="right" height=140/> --->

:information\_source: `bGWAS` has been updated to version 1.0.0.  
Check the [NEWS](NEWS.md) to learn more about what has been modified\!

:warning: If you downloaded the Z-Matrix files before 20/08/2019, they
are now obsolete and you will not be able to use them with the newest
version of the package.  
Note: some Prior GWASs have been removed, you can find more details
[here](docs/ZMatrices.md).

## Overview

bGWAS is an R-package to perform a Bayesian GWAS, using summary
statistics from a conventional GWAS as input. Briefly, it compares the
observed Z-scores for a focal phenotype to prior effects. These prior
effects are directly estimated from publicly available GWASs (currently,
a set of 38 studies, last update 20-08-2019 - hereinafter referred to as
“prior GWASs”). Only prior GWASs having a significant causal effect on
the focal phenotype, identified using a multivariate Mendelian
Randomization (MR) approach, are used to calculate the prior effects.
Causal effects are estimated masking the focal chromosome to ensure
independence.  
Observed and prior effects are compared using Bayes Factors.
Significance is assessed by calculating the probability of observing a
value larger than the observed BF (P-value) given the prior
distribution. This is done by decomposing the analytical form of the BFs
and using an approximation for most BFs to make the computation faster.
Prior, posterior and direct effects, alongside BFs and p-values are
returned. Note that prior, posterior and direct effects are estimated on
the Z-score scale, but are automatically rescaled to beta scale if
possible.

The principal functions available are:

  - **`bGWAS()`**  
    main function that calculates prior effects from prior GWASs,
    compares them to observed Z-scores and returns an object of class
    *bGWAS*

  - **`list_priorGWASs()`**  
    directly returns information about the prior GWASs that can be used
    to calculate prior effects

  - **`select_priorGWASs()`**  
    allows a quick selection of prior GWASs (to include/exclude specific
    studies when calculating prior effects)

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
    creates a Coefficients Plot (causal effect of each prior GWASs on
    the focal phenotype) from an object of class *bGWAS*

  - **`heatmap_bGWAS()`**  
    creates a heatmap to represent, for each significant SNP, the
    contribution of each prior GWAS to the estimated prior effect from
    an object of class *bGWAS*

All the functions available and more details about their usage can be
found in the [manual](doc/bGWAS-manual.pdf).

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
have alternative names:  

<ul>

SNP-identifier: `rs` or `rsid`, `snp`, `snpid`, `rnpid` Alternate
allele: `a1` or `alt`, `alts`  
Reference allele: `a2` or `a0`, `ref`  
Z-statistics: `z` or `Z`, `zscore`

</ul>

If the Z-statistics is not present, it will be automatically calculated
from effect size and standard error, in which case the following columns
should be provided:  

<ul>

Effect-size: `b` or `beta`, `beta1`  
Standard error: `se` or `std`

</ul>

If you want the prior/posterior/corrected effects to be rescaled, please
make sure to provide effect sizes and standard errors instead of (or in
addition to) Z-statistics.

#### 2\. Prior *GWASs* - Z-Matrix files

These files should be downloaded separately and stored in `~/ZMatrices`
or in the folder specified when launching the analysis. These files
contains the Z-scores for all prior GWASs :  

<ul>

<ul>

*ZMatrix\_MR.csv.gz*: Z-scores (strong instruments only) used for
multivariate MR,  
*ZMatrix\_Full.csv.gz*: Z-scores (all SNPs) used to calculate the prior
Z-scores,  
*AvailableStudies.tsv*: A file containing information about the prior
GWASs available.

</ul>

</ul>

<font color="red">**Problem with switch drive, download link will be
added soon\!**</font>

<!---
On UNIX/MACOSX, from a terminal:    
``` bash
wget --no-check-certificate ##https://drive.switch.ch/index.php/s/t3xepllvzobVTCh/download## -O ZMatrices.tar.gz
tar xzvf ZMatrices.tar.gz
``` 
On WINDOWS:   --->

<font color="grey"><small> If you want to use your own set of prior
GWASs, please have a look [here](doc/ZMatrices.md) to see how you can
modify the files.  
We focused on including prior GWASs that do not come from UKBB, assuming
that the focal phenotype results are more likely to be obtained from
UKBB. Sample overlap between the focal phenotype and the prior GWASs is
not accounted for by our method, so we did not include any UKBB results
in the prior GWASs.</font> </small>

### Study Selection

Before running your analysis, you can select the prior GWASs you want to
include. You can use the function **`list_prioGWASs()`** to get some
information about the prior GWASs available.

You should remove traits that by definition are not independent from
your trait of interest. For example, before analysing BMI results, make
sure to exclude “Height” from the prior GWASs used. You can use the
function **`select_priorGWASs()`** to automatically exclude/include some
traits or some files. You should also check for sample overlap, and
remove prior GWASs that come from the same consortium as your data. If
there are individuals in common between your conventional GWAS and prior
GWASs, it might induce some bias.

``` r
# Obtain the list of prior GWASs
AllStudies = list_priorGWASs()
# Select only the ones for specific traits
# select_priorGWASs will return the IDs of the files that are kept
MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body Mass Index", "Smoking"))
# Match these IDs against the ones in the list of prior GWASs 
AllStudies[AllStudies$ID %in% MyStudies, ]
```

    ## # A tibble: 6 x 10
    ##   File   Name      ID Trait  Consortium Reference Download  Remarks  N_SNPs
    ##   <chr>  <chr>  <dbl> <chr>  <chr>      <chr>     <chr>     <chr>     <dbl>
    ## 1 All_a… Body …     1 Body … GIANT      https://… https://… <NA>     6.81e6
    ## 2 META_… Heart…    23 Heart… HRgene     https://… http://w… <NA>     6.81e6
    ## 3 tag.c… Smoki…    35 Smoki… TAG        https://… https://… Other G… 6.63e6
    ## 4 tag.e… Smoki…    36 Smoki… TAG        https://… https://… Other G… 6.78e6
    ## 5 tag.f… Smoki…    37 Smoki… TAG        https://… https://… Other G… 6.78e6
    ## 6 tag.l… Smoki…    38 Smoki… TAG        https://… https://… Other G… 6.77e6
    ##   N_Instruments
    ##           <dbl>
    ## 1         10052
    ## 2          3229
    ## 3           558
    ## 4            94
    ## 5           184
    ## 6            34

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
list_priorGWASs(MyStudies)[,c("Name", "Trait", "Reference")]
```

    ## # A tibble: 6 x 3
    ##   Name                                 Trait                  
    ##   <chr>                                <chr>                  
    ## 1 Body Mass Index (GIANT)              Body Mass Index        
    ## 2 Coronary Artery Disease (CARDIoGRAM) Coronary Artery Disease
    ## 3 Years of Schooling (SSGAC)           Education              
    ## 4 Systolic Blood Pressure (ICBP)       Blood Pressure         
    ## 5 Diastolic Blood Pressure (ICBP)      Blood Pressure         
    ## 6 College Completion (SSGAC)           Education              
    ##   Reference                                   
    ##   <chr>                                       
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

    ## # A tibble: 1 x 3
    ##   Name                                 Trait                  
    ##   <chr>                                <chr>                  
    ## 1 Coronary Artery Disease (CARDIoGRAM) Coronary Artery Disease
    ##   Reference                                   
    ##   <chr>                                       
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
    ## # A tibble: 1 x 3
    ##   study                                             estimate std_error
    ##   <chr>                                                <dbl>     <dbl>
    ## 1 All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz  -0.0643    0.0320
    ## 
    ## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ 
    ##  
    ## 5 significant SNPs identified : 
    ##  rs17486195, rs1160985, rs1333049, rs12993295, rs78801493
    ## 
    ## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

``` r
print_log_bGWAS(A)
```

    ## <<< Preparation of analysis >>> 
    ## 
    ## > Checking parameters 
    ## The name of your analysis is: "Test_UsingSmallDataFrame". 
    ## The Z-Matrix files are stored in "/Users/nmounier/ZMatrices".  
    ## # Preparation of the data... 
    ## The conventional GWAS used as input the object: "GWAS".  
    ##    SNPID column, ok - ALT column, ok - REF column, ok - BETA column, ok - SE column, ok
    ## Posterior effects will be rescaled using BETA and SE.
    ## The analysis will be run in the folder: "/Users/nmounier/ZMatrices".  
    ## The p-value threshold used for selecting MR instruments is: 1e-06.  
    ## The minimum number instruments required for each trait is: 3.  
    ## The distance used for pruning MR instruments is: 500Kb.  
    ## Distance-based pruning will be used for MR instruments.  
    ## No shrinkage applied before performing MR.
    ## The p-value threshold used for stepwise selection is 0.05.  
    ## No shrinkage applied before performing calculating the prior.
    ## Significant SNPs will be identified according to p-value. The threshold used is :5e-08.  
    ## The distance used for pruning results is: 500Kb.  
    ## Distance-based pruning will be used for results.  
    ## 
    ## 
    ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ## <<< Identification of significant prior GWASs for MR >>>  
    ## 
    ## > Creating the Z-Matrix of strong instruments 
    ## # Loading the ZMatrix... 
    ## Selecting studies :
    ## 6 studies 
    ## 209,840 SNPs 
    ## # Adding data from the conventional GWAS :  "GWAS" 
    ## Done! 
    ## 9,125 SNPs in common between prior studies and the conventional GWAS 
    ## # Thresholding... 
    ## 786 SNPs left after thresholding 
    ## 6 studies left after thresholding 
    ## Pruning MR instruments... 
    ##    distance : 500Kb 
    ## 162 SNPs left after pruning 
    ## 6 studies left after thresholding+pruning 
    ## 
    ## > Performing MR  
    ## #Preparation of the MR analyses to identify significant studies... 
    ## Studies tested : Body Mass Index (GIANT) - Coronary Artery Disease (CARDIoGRAM) - Years of Schooling (SSGAC) - Systolic Blood Pressure (ICBP) - Diastolic Blood Pressure (ICBP) - College Completion (SSGAC)
    ## Conventionnal GWAS of interest : GWAS
    ## # Univariate regressions for each trait... 
    ##   Number of trait-specific instruments per univariate regression: 
    ##   . Body Mass Index (GIANT) : 60 
    ##   . Coronary Artery Disease (CARDIoGRAM) : 6 
    ##   . Years of Schooling (SSGAC) : 84 
    ##   . Systolic Blood Pressure (ICBP) : 5 
    ##   . Diastolic Blood Pressure (ICBP) : 8 
    ##   . College Completion (SSGAC) : 5 
    ## Done! 
    ## # Stepwise selection (all traits)... 
    ## Adding the first study :Body Mass Index (GIANT) 
    ## #Test if any study can be added with p<0.05 
    ## #Test if any study has p>0.05 now 
    ## It converged! 
    ## # Final regression... 
    ## The studies used are: 
    ## - Body Mass Index (GIANT)
    ## 
    ## Estimating adjusted R-squared: 
    ## - in-sample adjusted R-squared for the all-chromosomes multivariate regression is 0.0482 
    ## - out-of-sample R-squared (masking one chromosome at a time), for the multivariate regression will be estimated when calculating the prior. 
    ## 
    ## 
    ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ## <<< Estimation of the prior >>>  
    ## 
    ## > Creating the full Z-Matrix  
    ## # Loading the ZMatrix... 
    ## Selecting studies :
    ## 1 studies 
    ## 6,811,310 SNPs 
    ## # Adding data from the conventional GWAS :  "GWAS" 
    ## Done! 
    ## 291,583 SNPs in common between prior studies and the conventional GWAS 
    ## 
    ## > Computing prior  
    ## # Calculating the prior chromosome by chromosome... 
    ##    Chromosome 1
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 2
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 3
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 4
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 5
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 6
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 7
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 8
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 9
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 10
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 11
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 12
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 13
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 14
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 15
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 16
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 17
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 18
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 19
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 20
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 21
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ##    Chromosome 22
    ## Running regression, 
    ## Calculating prior estimates for SNPs on this chromosome, 
    ## Calculating prior standard errors for SNPs on this chromosome, 
    ## ## Out-of-sample R-squared for MR instruments across all chromosomes is 0.0161
    ## ## Out-of-sample squared correlation for MR instruments across all chromosome is 0.0489
    ## ## Correlation between prior and observed effects for all SNPs is 0.001
    ## ## Correlation between prior and observed effects for SNPs with GWAS p-value < 0.001 is 0.0319
    ## Done! 
    ## 
    ## 
    ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ## <<< Calculation of Bayes Factors and p-values >>>  
    ## 
    ## > Calculating them for all SNPs  
    ## # Computing observed Bayes Factor for all SNPs... 
    ## Done! 
    ## # Computing BF p-values... 
    ##    using a distribution approach: 
    ## ... getting approximated p-values using non-linear quantiles  
    ## ... checking p-values near significance threshold  
    ##     everything is ok!  
    ## Done! 
    ## > Pruning and identifying significant SNPs 
    ##    Starting with 291,583 SNPs 
    ## # Selecting significant SNPs according to p-values... 
    ## 22 SNPs left 
    ## Done! 
    ## # Pruning significant SNPs... 
    ##    distance : 500Kb 
    ## 5 SNPs left 
    ## Done! 
    ## 
    ## 
    ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ## Time of the analysis: 2 minute(s) and 2 second(s).

Functions to extract results from an object of class *bGWAS*:

``` r
hits = extract_results_bGWAS(A, "significant")
hits
```

    ## # A tibble: 5 x 22
    ##   chrm_UK10K pos_UK10K rsid  alt   ref      beta      se z_obs
    ##        <dbl>     <dbl> <chr> <chr> <chr>   <dbl>   <dbl> <dbl>
    ## 1          3  29790763 rs17… A     G     2.65e-4 2.40e-5 11.1 
    ## 2          1 220209813 rs11… T     C     7.79e-5 1.11e-5  7.04
    ## 3          6 132242660 rs13… C     G     1.04e-4 1.53e-5  6.77
    ## 4         15  89497100 rs12… T     G     4.34e-5 8.81e-6  4.92
    ## 5         10  85631344 rs78… T     C     6.16e-5 1.13e-5  5.46
    ##   mu_prior_estima… mu_prior_std_er… mu_posterior_es… mu_posterior_st…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.182              1.01             5.66            0.710
    ## 2           0.247              1.01             3.68            0.711
    ## 3          -0.0291             1.00             3.38            0.708
    ## 4           0.538              1.09             2.93            0.738
    ## 5           0.0335             1.00             2.75            0.708
    ## # … with 10 more variables: mu_direct_estimate <dbl>,
    ## #   mu_direct_std_error <dbl>, beta_prior_estimate <dbl>,
    ## #   beta_prior_std_error <dbl>, beta_posterior_estimate <dbl>,
    ## #   beta_posterior_std_error <dbl>, beta_direct_estimate <dbl>,
    ## #   beta_direct_std_error <dbl>, BF <dbl>, BF_p <dbl>

``` r
all_results = extract_results_bGWAS(A, "all")
nrow(all_results)
```

    ## [1] 291583

``` r
extract_MRcoeffs_bGWAS(A)[,1:12]
```

    ## # A tibble: 1 x 12
    ##   name       study           estimate std_error Tstat      P chrm1_estimate
    ##   <chr>      <chr>              <dbl>     <dbl> <dbl>  <dbl>          <dbl>
    ## 1 Body Mass… All_ancestries…  -0.0643    0.0320 -2.01 0.0490        -0.0782
    ##   chrm1_std_error chrm1_P chrm2_estimate chrm2_std_error chrm2_P
    ##             <dbl>   <dbl>          <dbl>           <dbl>   <dbl>
    ## 1          0.0323  0.0192        -0.0380          0.0312   0.229

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
300,000 SNPs in common with prior studies (see example A) \~ 2 minutes.

## Contact

<mounier.ninon@gmail.com>
