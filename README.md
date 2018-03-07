
# bGWAS
[//]:========================================

:information_source: package under development  
    
    
:warning: if you downloaded the Z-Matrix files before 25/02/2018, they are obsolete, you need to delete the old ones and download the new ones!   

## Overview
[//]:*******

bGWAS is an R-package to perform a Bayesian GWAS, using summary statistics as input. Briefly, it compares the observed Z-scores from a conventional GWAS to prior Z-scores. These prior Z-scores can be provided by the user or directly calculated from publicly available GWASs (currently, a set of 58 studies, last update 25-02-2018 - hereinafter referred to as "prior GWASs"). In this second scenario, only prior GWASs having a significant influence on the conventional GWAS (identified using a multivariate Mendelian Randomization (MR) approach) are used to calculate the prior Z-scores. Causal effect are estimated masking the focal chromosome to ensure independence.          
Observed and prior Z-scores are compared using Bayes Factors, and empirical p-values are calculated using a permutation approach.   

The main functions are:  

-   **`bGWAS()`**  
core function that calculates prior Z-scores from prior GWASs, compares them to observed Z-scores and returns an object of class "bGWAS"   

-   **`list_priorGWASs()`**
directly returns information about the prior GWASs that can be used to calculate prior Z-scores   

-   **`select_priorGWASs()`**   
allows a quick selection of prior GWASs (to include/exclude specific studies when calculating prior Z-scores)  

-   **`extract_results_bGWAS()`**    
returns results (prior estimate / standard-error + p-value from BF for SNPs) from an object of class "bGWAS"   

-   **`manhattan_plot_bGWAS()`**   
creates a Manhattan Plot from an object of class "bGWAS"  

-   **`extract_MRcoeffs_bGWAS()`**    
returns multivariate MR coefficients (1 estimate using all chromosomes + 22 estimates with 1 chromosome masked) from an object of class "bGWAS"   

-   **`coefficients_plot_bGWAS()`**   
creates a Coefficients Plot (causal effect of Prior GWASs) from an object of class "bGWAS"  

-   **`bGWAS_fromPrior()`**    
alternative function that compares prior Z-scores provided by the user to observed Z-scores and returns an object of class "bGWAS" # NOT IMPLEMENTED YET   



## Installation
[//]:*******

* Install R-package
``` r
# Directly install the package from github
# install.packages("devtools")
devtools::install_github("n-mounier/bGWAS")
```

If you want to use the **`bGWAS()`** function to calculate prior Z-scores, you should download Z-Matrix files. If you are planning to calculate your own prior Z-scores and analyse them using the **`bGWAS_fromPrior()`** you can skip this step.

* Download Z-Matrix files (Size ~ 2.09 GB):   
These files contains the Z-scores for all prior GWASs (before and after imputation) :   
Z-scores before imputation are used for multivariate MR,   
Z-scores after imputation are used to calculate the prior Z-scores.
``` bash
wget --no-check-certificate https://drive.switch.ch/index.php/s/BpRrDXvFPbnKCM6/download -O ZMatrices.tar.gz
tar xzvf ZMatrices.tar.gz
```
  
## Usage
[//]:*******

To run the analysis with `bGWAS` two inputs are needed:

#### 1. The *GWAS* results to be tested   
Can be a regular (space/tab/comma-separated) file or a gzipped file (.gz), must contain the following columns, which can have alternative names.  
SNP-identifier:  `rs` or `rsid`, `snp`, `snpid`, `rnpid`    
Alternate allele:  `a1` or `alt`, `alts`    
Reference allele: `a2` or `a0`, `ref`    
Z-statistics: `z` or `Z`, `zscore`      
If the Z-statistics is not present, it can be calculated from effect size and standard error, in which case the following columns should be provided:
Effect-size: `b` or `beta`, `beta1`    
Standard error:  `se` or `std`     

#### 2. Prior *GWASs*   
(see above for downloading instructions) - not needed for the function **`bGWAS_fromPrior()`**    
   
Matrix files, containing Z-scores for all prior GWASs should be downloaded separately and stored in `~/ZMatrices` or in the folder specified when launching the analysis.   

*Currently, this is not possible for users to directly add their own GWASs to this set of prior GWASs. Feel free to contact us for more information.*   
 
<!---  Format?
Can I add one more?--->


### Study Selection
[//]:-------------------------------

Before running your analysis, you can select the prior GWASs you want to include. You can use the function **`list_prioGWASs()`** to get some information about the prior GWASs available.   
You should remove traits that by definition are not independent from your trait. For example, before analysing BMI results, make sure to exclude "Height" from the prior GWASs used. You can use the function **`select_priorGWASs()`** to automatically exclude/include some traits or some files.   
*+ also check consortium (no function for selection) but if same individuals in your GWAS + prior GWASs, not ok*

``` r
AllStudies = list_priorGWASs()
MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body mass index", "Smoking"))
AllStudies[AllStudies$ID %in% MyStudies, ]
```

### Analysis
[//]:-------------------------------

``` r
## Example A
# Using a GWAS from our list our prior GWASs
# Using all other (57) GWASs to built the prior
MyGWAS = 5
list_priorGWASs(MyGWAS)
# Coronary Artery Disease GWAS (CARDIoGRAM)

A = bGWAS(name = "Test_UsingGWASfromPriorGWASs",
         GWAS = MyGWAS,
         verbose=T)
# MR instruments will be selected using default parameters,
# All Prior GWASs except the one use as "GWAS" will be used to create the prior,
# Significant SNPs will be identified using default parameters (p<5e-8) and will not be pruned,
# No file will be saved.
``` 

``` r
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
# MR instruments will be selected using default parameters,
# A subset of Prior GWASs will be used to create the prior,
# Significant SNPs will be identified using default parameters (p<5e-8) and will not be pruned,
# No file will be saved.         
         
      
print(B)
manhattan_plot_bGWAS(B)
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
# MR instruments will be selected using default parameters,
# A subset of Prior GWASs will be used to create the prior,
# Significant SNPs will be identified using default parameters (p<5e-8) and will not be pruned,
# Files will be saved.         
         
print(C)
         
         
# Note that B and C are using the same data (stored differently) and give the same results.
all.equal(B,C)
         
```

### Results
[//]:-------------------------------


**`bGWAS()`** returns an object of class "bGWAS" than can be handled in `R`.    

```r
## Results from example A

class(A)
# "bGWAS"
print(A)
# -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#
#  bGWAS performed on 7,104,649 SNPs
#
# -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#
# 7 studies used to build the prior :
#                                                      Study    Estimate    StdError
#          All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz  0.06707833  0.016167132
#                                    DIAGRAMv3.2013MAY07.zip  0.08591397  0.026578888
#                                  EDUyears_2016_sumstat.txt -0.03609801  0.016561887
#                                     jointGwasMc_LDL.txt.gz  0.12351084  0.011990610
#                                      jointGwasMc_TG.txt.gz  0.07034726  0.013960301
#  MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt...BMI  0.14637027  0.057695521
#                                             MetaSum.height  0.02236071 0.005837044
#
# -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
#
# 19 significant SNPs identified :
# rs599839, rs964184, rs56289821, rs11591147, rs944801, rs2954029, rs4953023, rs137897428, rs55730499, rs579459,
# rs7578326, rs2410621, rs62120566, rs34967177, rs668948, rs2228603, rs1169288, rs2148194, rs139794538
#
# -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
 
```
Functions to extract results from an object of class "bGWAS":   
```r
hits = extract_results_bGWAS(myObj, "significant")
print(hits)
# rs         chrm   pos   alt ref    obs         fit        se        log.bf  nullls_count     BF_p
# rs599839    1   109822166  A   G   6.304567   4.2059161 0.5362215  18.037043    0     1.407529e-10
# rs964184    11  116648917  C   G  -6.144531  -3.3559249 0.5488481  15.757918    0     1.407529e-10
# rs56289821  19  11188247   A   G  -5.023050  -3.3812354 0.4482898  11.401707    0     1.407529e-10
# rs11591147  1   55505647   T   G  -5.087970  -3.0578829 0.4256008  11.115867    0     1.407529e-10
# rs944801    9   22051670   C   G   10.054911  0.5824288 0.2967470   9.275391    1     2.815058e-10
# rs2954029   8   126490972  T   A  -4.189963  -3.5102907 0.4525705   8.493019    3     5.630116e-10
# rs4953023   2   44074000   A   G  -4.738788  -2.2908836 0.3391445   8.486565    3     5.630116e-10
# rs137897428 10  60793223   T   G   5.688300   1.2563294 0.5323997   8.401426    3     5.630116e-10
# rs55730499  6   161005610  T   C   7.384310   0.9791115 0.2818131   8.221794    3     5.630116e-10
# rs579459    9   136154168  C   T   5.299689   1.6495696 0.3000976   7.888927    7     1.126023e-09
# rs7578326   2   227020653  G   A  -4.736911  -1.5943275 0.6309521   7.519732    15    2.252047e-09
# rs2410621   8   19854682   C   T  -4.337672  -1.9841644 0.5199402   7.107929    30    4.363340e-09
# rs62120566  19  45197732   G   C  -4.126030  -2.3491911 0.3593355   7.053283    34    4.926352e-09
# rs34967177  7   45176432   A   G   12.500800  0.1682764 0.2636460   6.998161    37    5.348610e-09
# rs668948    2   21291529   A   G   3.663595   3.7404512 0.4433215   6.618781    65    9.289692e-09
# rs2228603   19  19329924   T   C  -3.916886  -2.2765867 0.4031095   6.438471    90    1.280851e-08
# rs1169288   12  121416650  C   A   4.639650   1.5821789 0.3014183   6.434916    90    1.280851e-08
# rs2148194   10  3975619    T   C   7.013320   0.7858808 0.2603974   6.401342    93    1.323077e-08
# rs139794538 11  116078910  T   C  -7.526290  -0.5417488 0.2736795   5.594159    285   4.025533e-08

all_results =extract_results_bGWAS(A, "all")
nrow(all_results)
# 7104649
head(all_results)
# rs         chrm   pos   alt ref    obs         fit        se        log.bf  nullls_count     BF_p
# rs599839    1   109822166  A   G   6.304567   4.2059161 0.5362215  18.037043    0     1.407529e-10
# rs964184    11  116648917  C   G  -6.144531  -3.3559249 0.5488481  15.757918    0     1.407529e-10
# rs56289821  19  11188247   A   G  -5.023050  -3.3812354 0.4482898  11.401707    0     1.407529e-10
# rs11591147  1   55505647   T   G  -5.087970  -3.0578829 0.4256008  11.115867    0     1.407529e-10
# rs944801    9   22051670   C   G   10.054911  0.5824288 0.2967470   9.275391    1     2.815058e-10
# rs2954029   8   126490972  T   A  -4.189963  -3.5102907 0.4525705   8.493019    3     5.630116e-10
```

```r
MR_coefficients = extract_MRcoeffs_bGWAS(A)
colnames(MR_coefficients)
# "Study"            "Estimate"        "StdError"        "T"     "P"
# "chrm1_Estimate"   "chrm1_StdError"  "chrm1_PValue"   "chrm2_Estimate"   "chrm2_StdError"   "chrm2_PValue"    
# "chrm3_Estimate"   "chrm3_StdError"  "chrm3_PValue"   "chrm4_Estimate"   "chrm4_StdError"   "chrm4_PValue"    
# "chrm5_Estimate"   "chrm5_StdError"  "chrm5_PValue"   "chrm6_Estimate"   "chrm6_StdError"   "chrm6_PValue"
# "chrm7_Estimate"   "chrm7_StdError"  "chrm7_PValue"   "chrm8_Estimate"   "chrm8_StdError"   "chrm8_PValue"
# "chrm9_Estimate"   "chrm9_StdError"  "chrm9_PValue"   "chrm10_Estimate"  "chrm10_StdError"  "chrm10_PValue"
# "chrm11_Estimate"  "chrm11_StdError" "chrm11_PValue"  "chrm12_Estimate"  "chrm12_StdError"  "chrm12_PValue"   
# "chrm13_Estimate"  "chrm13_StdError" "chrm13_PValue"  "chrm14_Estimate"  "chrm14_StdError"  "chrm14_PValue"
# "chrm15_Estimate"  "chrm15_StdError" "chrm15_PValue"  "chrm16_Estimate"  "chrm16_StdError"  "chrm16_PValue"   
# "chrm17_Estimate"  "chrm17_StdError" "chrm17_PValue"  "chrm18_Estimate"  "chrm18_StdError"  "chrm18_PValue"
# "chrm19_Estimate"  "chrm19_StdError" "chrm19_PValue"  "chrm20_Estimate"  "chrm20_StdError"  "chrm20_PValue"   
# "chrm21_Estimate"  "chrm21_StdError" "chrm21_PValue"  "chrm22_Estimate"  "chrm22_StdError"  "chrm22_PValue"

MR_coefficients[1:4, 1:8]
#                                            Study    Estimate   StdError        T       P          chrm1_Estimate  chrm1_StdError  chrm1_PValue
# All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz  0.06707833 0.01616713   4.149056 3.479263e-05     0.07218810    0.01722797  2.921263e-05
#                           DIAGRAMv3.2013MAY07.zip  0.08591397 0.02657889   3.232414 1.247588e-03     0.08567881    0.02664810  1.326470e-03
#                         EDUyears_2016_sumstat.txt -0.03609801 0.01656189  -2.179583 2.940545e-02    -0.03737184    0.01738187  3.168322e-02
#                            jointGwasMc_LDL.txt.gz  0.12351084 0.01199061   0.300630 2.821272e-24     0.12018871    0.01357595  1.978510e-18
```


Functions for graphic representations:   
```r 
manhattan_plot_bGWAS(A)
```
![](https://drive.switch.ch/index.php/apps/files_sharing/ajax/publicpreview.php?x=2506&y=584&a=true&file=ManhattanPlot.png&t=LQ8bk8p5PZFUCho&scalingup=0)  
```r 
coefficients_plot_bGWAS(A) 
```
![](https://drive.switch.ch/index.php/apps/files_sharing/ajax/publicpreview.php?x=2506&y=584&a=true&file=CoefficientsPlot.png&t=PnfY6EC7R1SoLUe&scalingup=0)  
   
   
##### Aditionnaly, if `save_files=T`, several files are created...   
... in the working directory :   
-   **`name`.log** - log file   

... in the folder `./name/` :   
-   **PriorGWASs.tsv** - contains Prior GWASs information (general info + status (used/removed) + MR coefficients)   
-   **CoefficientsByChromosome.csv** - contains the MR estimates when masking the focal chromosome (22 coefficients / study)   
-   **PriorBFp.csv** - contains BF and p-values estimated for all SNPs   
-   **SignificantSNPs.csv** - contains BF and p-values estimated for a subset of SNPs    


-> Further description of the files  


## Runtime
[//]:*******

Analysis using all the 58 prior GWASs available, for a conventional GWAS containing ~7M SNPs in common with the prior studies ~ 145 minutes.

Analysis using 9 prior GWASs, for a conventional GWAS containing 400,000 SNPs in common with prior studies (see example B) ~ 8 minutes



## Contact
<mounier.ninon@gmail.com>

    
      

       
       
<!---    
## Improvements to be implemented
[//]:-------------------------------

document results files    

bGWASfromPrior()    

selection from consortium (maybe not possible - complicated, too many cohorts in each consortium...)    
use of a subset of SNPs   
use of re-imputed studies for prior   
use of additional studies for prior   

optimize null-BF calculation

--->  




