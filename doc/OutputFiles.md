# Output files
[//]:========================================

If `save_files=T`, several files are created in the folder `./name/`.

## Description
[//]:*******

The output files created are:    
-   **name.log** - report containing everything that has been printed during the bGWAS analysis (is `verbose=T`)     
    
    
-   **PriorGWASs.tsv** - contains information about all prior GWASs (general info + status (used/removed) + MR coefficients):
<ul>    File : prior GWAS file name <br/>      
        Name : prior GWAS name  <br/>   
        Trait : phenotypical category (both Systolic and Diastiolic Blood Pressure are part of "Blood Pressure" for example)   <br/> 
        status : has the file been used to create the prior? (can be "USED", "Excluded by user", "Excluded for MR: removale reason")    <br/>
        uni_estimate : causal effect estimate from univariate MR     <br/>
        uni_std_error	: causal effect standard error from univariate MR    <br/>
        uni_T	: T-statistic from univariate MR    <br/>
        uni_P : p-value from univariate MR    <br/>
        uni_adj_Rsquared : adjusted R squared from univariate MR     <br/>
        uni_Rsquared : R squared from univariate MR     	<br/>
        multi_estimate : causal effect estimate from multivariate MR     <br/>
        multi_std_error	: causal effect standard error from multivariate MR   <br/>  
        multi_T	: T-statistic from multivariate MR     <br/>
        multi_P: p-value estimate from multivariate MR   <br/>  </ul>
-   **CoefficientsByChromosome.csv** - contains the MR estimates when masking the focal chromosome (22 coefficients / prior GWASs used for prior estimation)   
<ul>    study : everything that is printed during a bGWAS analysis (if `verbose=TRUE`) <br/>     
        estimate : causal effect estimate from multivariate MR (masking one chromosome)    <br/>
        std_error : causal effect standard error from multivariate MR (masking one chromosome)  <br/>  
        T : T-statistics from multivariate MR (masking one chromosome)   <br/> 
        P : p-value from multivariate MR (masking one chromosome)  <br/>   
        chrm : chromosome masked     <br/> </ul>
        -   **PriorBFp.csv** - contains BF and p-values, prior, posterior and direct effects estimates for all SNPs      
<ul>    chrm_UK10K : chromosome (obtained from UK10K data)     <br/>
        pos_UK10K : position (obtained from UK10K data)   <br/> 
        rsid : rs number     <br/>
        alt : alternative (effect) allele     <br/>
        ref : reference allele     <br/>
        beta : observed effect size    <br/> 
        se : observed standard error   <br/>  
        z_obs : observed Z-score    <br/>
        mu_prior_estimate : prior effect estimate (z-score scale)     <br/>
        mu_prior_std_error : prior effect standard error (z-score scale)    <br/>
        mu_posterior_estimate : posterior effect estimate (z-score scale)    <br/>
        mu_posterior_std_error : posterior effect standard error (z-score scale)   <br/> 
        mu_direct_estimate : direct effect estimate (z-score scale)   <br/>
        mu_direct_std_error : direct effect standard error (z-score scale)   <br/>
        beta_prior_estimate : prior effect estimate (beta scale)   <br/>
        beta_prior_std_error : prior effect standard error (beta scale)  <br/>  
        beta_posterior_estimate : posterior effect estimate (beta scale)   <br/> 
        beta_posterior_std_error : posterior standard error estimate (beta scale)   <br/> 
        beta_direct_estimate : direct effect estimate (beta scale)   <br/>
        beta_direct_std_error : direct effect standard error (beta scale)   <br/>
        BF : Bayes Factor        <br/>
        BF_p/BF_fdr : Bayes Factor p-value / fdr  <br/>  </ul>
        -   **SignificantSNPs.csv** - contains BF and p-values, prior, posterior and direct effects estimates for a subset of significant SNPs (subset of **PriorBFp.csv**)  

