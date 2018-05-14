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
