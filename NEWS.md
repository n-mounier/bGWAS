# bGWAS 0.3.1 (2018-10-04)

## Changes
- Modification of the selection of traits in the multivariate MR model

The matrix of instruments use for each model should be specific to the traits used by the model, i.e. all the models we want to compare will have a different number of observations, and we can't use AIC to compare them anymore.    
We implemented a stepwise approach, based on p-value to decide if a trait should be included or not. We start with a model including only the most significant study (univariate), and try to add, one by one, all the remaining studies (and update the instruments used for each model). If one of the studies has a p-value < 0.05 / number of studies tested, the most significant one is added to the model (if the univariate and multivariate effect estimate are in the same direction, otherwise the study is discarded). The model is updated to include these 2 studies, and we test if any of the studies in the updated model has a p-value > 0.05 / number of studies tested. We keep adding/removing studies until we obtain a definitive set of studies (convergence).    

- Modification of the out-of-sample R2

Out-of-sample R2 is not corrected for the number of traits/observation (no adjustement needed), and is calculated once all the regressions masking one chromosome are performed. We now obtain a global out-of-sample R2, and not one per chromosome (no need to use the mean/median, more robust if some chromosome have only a small number of SNPs to be predict).   



# bGWAS 0.3.0 (2018-07-31)

## Changes
- Modification of the way BF p-values are estimated 

We modified the already implemented distribution approach to be usable even if no shrinkage has been applied.   
This approach uses a quantile approximation for all SNPs. The quantiles used are non-linearly distributed to get a better approximation (more quantiles for low/high prior values, because using only a 100 linealrly distributed quantiles (percentiles), we can't capture correctly the effect of low/high prior values on the p-values - but the weight are proportionnal to the number of prior values falling in these intervals). For SNPs near significance threshold (potentially false positive/negative), the full formula is used to check that the results are correct, and to re-estimate the p-values more precisely if needed.    
This method is fast (between 15 and 30 minutes, depending of the dataset), more accurate than the permutation approach (for small p-values in particular) and reproducible (no randomness involved).

- Modification of the way Z-Matrices are handled   

The file "AvailableStudies.tsv"" is now stored in the Z-Matrix folder (it was previously included in the package). It contains basic information about the prior GWASs that can be used to create the prior. It is mainly used by the function 'select_priorGWASs()', but it is also helpful to select the columns needed when loading the Z-Matrices (and avoiding loading the full matrices uselessly). The lines of this file must match the columns names from the Z-Matrices (i.e. correctly describe the prior GWASs) and this is automatically checked by the package.   
Taking this file out from the package allow an easy external modification of the prior GWASs used if an user needs to include some other studies to the set of prior GWASs available.   


- Modification of the distance pruning method

Minor modification to decrease runtime.


# bGWAS 0.2.0

## Changes
Modification of the way BF p-values are estimated (distribution approach with approximation instead of permutation) to make the analyses faster and the estimation of small p-values more accurate.  

This approach is a lot faster when the large majority of SNPs have a zero-prior. As the number of SNPs with a non-zero prior increase, so does the runtime. When using the permutation approach to get empirical p-values, the runtime only depended on the total number of SNPs analyzed and could be easily predicted. Here, both the number of SNPs with a non-zero prior and the distribution of these prior effects influence the runtime so it cannot be predicted. However, in all our tests, this approach remained faster than the pertutation approach previously implemented.


# bGWAS 0.1.0 (2018-03-15)

initial version of the package, as described in McDaid et al (2017).

<!--- 
## Bug fixes

## New functions

## Documentation

## Error messages

## Performance


--->  
