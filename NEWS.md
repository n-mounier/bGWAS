# bGWAS 1.0.1 (2019-11-04)
    

## Changes
- Variance of the prior effect    

The variance of the prior effect was wrongly estimated (+1 unnecessarily added). The correct standard errors are now returned. This may affect BFs and p-values, but also posterior and direct effect estimates / standard errors.

- Corrected to raw ratio (CRR)    

The ratio between direct effect and observed effect (also named corrected to raw ratio or CRR) has been added to the output files.

# bGWAS 1.0.0 (2019-08-26)
    

## Changes
- Posterior, corrected effects and rescaling    

Now, the posterior effects (`mu_posterior_estimate` - `mu_posterior_std_error`) and the direct effects (`mu_direct_estimate` - `mu_direct_std_error`) are estimated simultaneously with the prior effects (`mu_prior_estimate` - `mu_prior_std_error`). The analysis is performed on the z-score scale, and these effects are automatically rescaled to the observed effect size beta scale when possible (`beta_posterior_estimate`, `beta_posterior_std_error`, `beta_direct_estimate`, `beta_direct_std_error`).   

- Stepwise selection : `stepwise_threshold` and studies tested       

By default, the stepwise selection threshold was set to : 0.05 / (the number of prior GWASs used). Now, we allow the user to provide a custom threshold (that could be less stringent, since some of the prior GWASs might be correlated).  
We also changed the subset of studies included in the stepwise selection procedure: to increase the robustness of the study identification, only the ones reaching nominal significance in univariate regression are included and tested (note that this does not affect the `stepwise_threshold` value).    
If convergence is not achieved after 50 iterations in the stepwise procedure, there is probably a loop (adding/removing the same studies iteratively). In such cases, the analysis is stopped and the user is invited to relaunch the analysis with a different set of `prior_studies` to avoid this behavior.     

- Update of the Z-matrices files    

We updated the Z-matrices files, to reduce the set of prior GWASs and exclude some files that were not meeting our inclusion criteria (see [here](doc/ZMatrices.md) for more details). We also added a column "Name" to the file `AvailableStudies.tsv` to use in the figures and increase readability.
   

- Minimum number of instruments for MR    

Initially, all prior GWASs with at least two instruments were kept to be tested in the multivariate model. This number was arbitrary, and we now allow the user to specify the minimum number of instruments that should be used  (`MR_ninstruments`).    

- Pruning of MR instruments    

When the default threshold was used to select MR instruments, the step to "select strong instruments" before pruning them was skipped, assuming that all the SNPs in the matrix were already strong instruments. That was the case only when all the prior GWASs were used, and if a subset of prior GWASs was selected, pruning was done based on all prior GWASs, not only the ones included. This is now fixed.      

## New functions
- `heatmap_bGWAS()`   

This function creates a heatmap to represent the contribution of each risk factor to the prior effect estimated (for all significant SNPs).   

- `get_RSquared_bGWAS()`   

This function returns the (out-of-sample) squared correlation between prior and observed effects (for all SNPs, the ones having a moderate effect on the trait, or MR instruments only).   

- `print_log_bGWAS()`   

This function prints the log (summary of the \code{bGWAS} analysis performed, everything that is printed using \code{verbose=TRUE}).   



## Documentation
- Description of Z-matrices files

We added a file to provide a more in depth description of the Z-matrices files [here](docs/ZMatrices.md).    

- Description of output files    

We added a file to provide a more in depth description of the outpute files [here](docs/OutputFiles.md).    


## Performance

- Windows Compatibilty 

Takes advantage of the ability to directly read `.gz` files using `data.table::fread` (no need to use `system("zcat < ...")` anymore.    


# bGWAS 0.3.2 (2019-04-17)

## Bug fixes    
- Correction of the instruments selection process        

When only one prior GWAS is used, the previous version was using z-score > threshold instead of abs(z-score) > threshold to select instruments. This could lead to a smaller number of instruments used and a potentially incorrect estimation of the causal effects, and the R2 / squared correlation values reported.        
    

## Changes   
- Clarification of the out-of-sample R2 and squared correlation estimation reported 

These values are estimated using only the instruments initially selected for the multivariate MR model (not all SNPs). The predicted values of the instruments obtained when masking the chromosomes they are located on are used to calculate the out-of-sample R2 and the out-of-sample squared correlation. This is now more explicitly indicated in the log file.     
Note that the R2 and the squared correlation can differ because the regression model does not include an intercept.       

- Information about the correlation between prior and observed effects    

In addition to the out-of-sample R2 and squared correlation, the correlation between prior and observed effects is now also reported (for all SNPs + for SNPs moderately associated with the trait, p<0.001).     


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
