# Output files
[//]:========================================

If `save_files=T`, several files are created in the folder `./name/`.

## Description
[//]:*******

The output files created are:    
-   **name.log** - report containing everything that has been printed during the bGWAS analysis (is `verbose=T`)     
    
    
-   **PriorGWASs.tsv** - contains information about all prior GWASs (general info + status (used/removed) + MR coefficients):
<ul>    File : c      
        Name : prior GWAS name     
        Trait : phenotypical category (both Systolic and Diastiolic Blood Pressure are part of "Blood Pressure" for example)    
        status : has the file been used to create the prior? (can be "USED", "Excluded by user", "Excluded for MR: removale reason")    
        uni_estimate : causal effect estimate from univariate MR     
        uni_std_error	: causal effect standard error from univariate MR    
        uni_T	: T-statistic from univariate MR    
        uni_P : p-value from univariate MR    
        uni_adj_Rsquared : adjusted R squared from univariate MR     
        uni_Rsquared : R squared from univariate MR     	
        multi_estimate : causal effect estimate from multivariate MR     
        multi_std_error	: causal effect standard error from multivariate MR     
        multi_T	: T-statistic from multivariate MR     
        multi_P: p-value estimate from multivariate MR     </ul>
-   **CoefficientsByChromosome.csv** - contains the MR estimates when masking the focal chromosome (22 coefficients / prior GWASs used for prior estimation)   
<ul>    study : everything that is printed during a bGWAS analysi    
        estimate : causal effect estimate from multivariate MR (masking one chromosome)    
        std_error : causal effect standard error from multivariate MR (masking one chromosome)    
        T : T-statistics from multivariate MR (masking one chromosome)    
        P : p-value from multivariate MR (masking one chromosome)     
        chrm : chromosome masked      </ul>
        -   **PriorBFp.csv** - contains BF and p-values, prior, posterior and direct effects estimates for all SNPs      
<ul>    chrm_UK10K : chromosome (obtained from UK10K data)     
        pos_UK10K : position (obtained from UK10K data)    
        rsid : rs number     
        alt : alternative (effect) allele     
        ref : reference allele     
        beta : observed effect size     
        se : observed standard error     
        z_obs : observed Z-score    
        mu_prior_estimate : prior effect estimate (z-score scale)     
        mu_prior_std_error : prior effect standard error (z-score scale)    
        mu_posterior_estimate : posterior effect estimate (z-score scale)    
        mu_posterior_std_error : posterior effect standard error (z-score scale)    
        mu_direct_estimate : direct effect estimate (z-score scale)   
        mu_direct_std_error : direct effect standard error (z-score scale)   
        beta_prior_estimate : prior effect estimate (beta scale)   
        beta_prior_std_error : prior effect standard error (beta scale)    
        beta_posterior_estimate : posterior effect estimate (beta scale)    
        beta_posterior_std_error : posterior standard error estimate (beta scale)    
        beta_direct_estimate : direct effect estimate (beta scale)   
        beta_direct_std_error : direct effect standard error (beta scale)   
        BF : Bayes Factor        
        BF_p/BF_fdr : Bayes Factor p-value / fdr    </ul>
        -   **SignificantSNPs.csv** - contains BF and p-values, prior, posterior and direct effects estimates for a subset of significant SNPs (subset of **PriorBFp.csv**)  

