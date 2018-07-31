# bGWAS 0.3.0 (2018-07-31)

## Changes
- Modification of the way BF p-values are estimated 

We modified the already implemented distribution approachto be usable even if no shrinkage has been applied.   
This approach uses a percentile approximation for SNPs with low BF, the full formula in the area near significance to avoid false positives/negatives, and a better approximation (using more quantiles) for significant SNPs, to ensure that the highest a BF, the smallest the respective p-value.   
This method is fast (between 15 and 30 minutes, depending of the dataset), more accurate than the permutation approach () and flexible, since it automatically determine the number of quantiles needed to keep the BF/p-value ordering consistent.

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
