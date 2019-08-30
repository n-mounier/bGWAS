###### Function to calculate BF and p-values ######



# #' Calculate Bayes Factor and empirical p-values from observed Z-scores and priors
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED




request_BFandP <- function(Prior, sign_thresh, use_permutations = F, 
                           sign_method, save_files=F, verbose=F) {
  Log = c()
  tmp = paste0("# Computing observed Bayes Factor for all SNPs... \n")
  Log = update_log(Log, tmp, verbose)
  
  
  
  # calculate BF
  Prior %>%
    mutate(BF = stats::dnorm(mean= .data$mu_prior_estimate, 
                             sd=sqrt(1 + .data$mu_prior_std_error**2), 
                             x = .data$z_obs) /
             stats::dnorm(mean = 0.0, 
                   sd=sqrt(1), 
                   x = .data$z_obs)) %>%  # and order by BF
    arrange(desc(.data$BF)) -> Prior
  
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  # true Z = just rs + log.bf ???
  
  
  
  tmp = "# Computing BF p-values... \n"
  Log = update_log(Log, tmp, verbose)
  
  if(use_permutations){# The true one : calculate BF
    tmp = "   using a permutation approach: \n"
    
    
    # how many null z-scores should we simulate?
    nrow(Prior) -> N.one.set.of.nulls
    
    # Function to count how many larger BF were observed in the simulation
    Rcpp::cppFunction(' IntegerVector count_larger(NumericVector real, NumericVector null) {
                      int N_real = real.size();
                      int N_null = null.size();
                      // both files are in descending order. Make use of this
                      IntegerVector result(N_real);
                      int j = 0; // the observed one is smaller than (not equal) this many nulls
                      for(int i=0; i<N_real; ++i) {
                      double real_i = real.at(i);
                      while(j < N_null && null.at(j) > real_i) {
                      ++j;
                      }
                      result.at(i) = j;
                      }
                      return result;
  } ')
    
    
    
    # counts for each SNPs : start with all 0
    running.total = rep.int(0 , nrow(Prior))
    
    tmp = "# Computing null Bayes Factors (needed to compute empirical p-values)... \n"
    Log = update_log(Log, tmp, verbose)
    
    
    NUMBER_OF_NULLRUNS = 1000
    full.number.of.nulls.for.comparison = NUMBER_OF_NULLRUNS* N.one.set.of.nulls
    
    tmp = paste0(NUMBER_OF_NULLRUNS, " random z-scores simulations for the whole set of SNPs (",
                 format(nrow(Prior), big.mark = ",", scientific = F), ") : ", format(full.number.of.nulls.for.comparison, big.mark = ",", scientific = F), " null BFs are calculated \n")
    Log = update_log(Log, tmp, verbose)
    
    for( i in 1:NUMBER_OF_NULLRUNS ) {
      # null bf
      if(i %% 20 == 1){
        tmp = paste0("Simulation ",i, "... \n")
        Log = update_log(Log, tmp, verbose)
        
      }
      
      # simulate z-score from normal distribution
      # ... but this does not account for the real "true distribution" (inflation?)
      # should we shuffle z-scores instead?
      # Prior$z <-   rnorm(nrow(Prior), 0, 1)
      Prior$z <- Prior$z_obs[sample(1:nrow(Prior))]
      
      Prior$BF.null =  stats::dnorm(mean= Prior$mu_prior_estimate, 
                                    sd=sqrt(1+Prior$mu_prior_std_error**2), 
                                    x=Prior$z) /
        stats::dnorm(mean= 0.0, 
                     sd=sqrt(1),
                     x=Prior$z)
      
      TempP = Prior[order(-Prior$BF.null),]
      nullbf = TempP$BF.nul
      
      N = length(nullbf)
      
      stopifnot(N == N.one.set.of.nulls)
      #(!is.unsorted(-nullbf)) %|%stopifnot
      
      count_larger(Prior$BF, nullbf) -> counts
      stopifnot(length(counts) == nrow(Prior))
      running.total = running.total + counts
    }
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
    
    tmp = "# Computing empirical p-values... \n"
    Log = update_log(Log, tmp, verbose)
    
    
    Prior$z = NULL
    Prior$BF.null = NULL
    
    Prior$BF_P = (running.total+1) / full.number.of.nulls.for.comparison
    
    
  } else {
    tmp = "   using a distribution approach: \n"
    Log = update_log(Log, tmp, verbose)
    
    
    ### Functions to estimate p-value from BF/prior distribution
    ## Calculate a p-value using the non-linear quantiles
    
    # Get p-value, for one BF value, using one vector of quantiles, one vector of weights, and one variance value
    get_OneP <- function(myBF, perc, sigma2, weights){
      suppressWarnings({ # warnings when p-value = NA (small BFS), deal with it below
        pPercentiles = stats::pnorm( (-perc/sigma2) +
                                sqrt(1+sigma2)/sigma2 *
                                sqrt(perc**2 + 2*sigma2 * log(myBF * sqrt(1+sigma2))), lower.tail=F) +
          stats::pnorm((-perc/sigma2) -
                  sqrt(1+sigma2)/sigma2 *
                  sqrt(perc**2 + 2*sigma2 * log(myBF * sqrt(1+sigma2))), lower.tail=T)
      })
      pPercentiles[is.na(pPercentiles)]=1
      # multiply each p-value by the according weight
      p = pPercentiles %*% weights
      return(p)
    }
    
    ## Get p-value, for one BF value, using all priors ("true value")
    get_fullP <- function(myBF, Data){ # BF : value - Data : data.frame with prior_estimate and prior_std_error columns
      suppressWarnings({ # warnings when p-value = NA (small BFS), deal with it below
        Data %>%
          mutate( pSNP = stats::pnorm( (-.data$mu_prior_estimate/(.data$mu_prior_std_error**2)) +
                                  sqrt(1+.data$mu_prior_std_error**2)/(.data$mu_prior_std_error**2) *
                                  sqrt(.data$mu_prior_estimate**2 + 2*.data$mu_prior_std_error**2 * log(myBF * sqrt(1+.data$mu_prior_std_error**2))), lower.tail=F) +
                    stats::pnorm((-.data$mu_prior_estimate/(.data$mu_prior_std_error**2)) -
                            sqrt(1+.data$mu_prior_std_error**2)/(.data$mu_prior_std_error**2) *
                            sqrt(.data$mu_prior_estimate**2 + 2*.data$mu_prior_std_error**2 * log(myBF * sqrt(1+.data$mu_prior_std_error**2))), lower.tail=T)) -> Data
      })
      Data %>%
        mutate(pSNP = tidyr::replace_na(.data$pSNP, 1)) %>%  # if NA -> replace by 1, that means that the polynome has no real solution and Pr(allBF>BF)=1
        summarize(mean(.data$pSNP)) -> p # here, we use all the prior values, so we can use an unweighted mean
      return(p)
    }
    
    # Get all p-values for a vector of BFs, return a list with 2 elements : "pValue" & "log_info"
    get_pValue<- function(myBFs, Data, sign_thr=5e-8, verbose=F){ # BF : vector of values - Data : data.frame with prior_estimate and prior_std_error columns
      Log=c()
      
      tmp = "... getting approximated p-values using non-linear quantiles  \n"
      Log = update_log(Log, tmp, verbose)
      
      d=data.frame(BF=myBFs)
      d$my_p=NA
      
      # here, we define quantiles on the 10^... scale
      # from 0 to 0.5, and then from 0.5 to 1
      quant = c(0, 10^(seq(-10, 0, 0.1))/2  , 1-rev(10^(seq(-10, 0, 0.1))/2)[-1], 1)
      
      
      #c(0, 10^(seq(-10, 0, 10/49))/2  , 1-rev(10^(seq(-10, 0, 10/49))/2)[-1])
      
      
      # for the weights, we take quant[x] - quant[x-1] to get the "size" of the intervals
      # and then, set the weight of a "quantile value" as half of the "size" of the interval before 
      #                                                 + half of the "size" of the interval after
      sizes = quant[-1]-quant[-length(quant)]
      wghts = numeric(length(quant))
      for(i in 1:length(quant)){
        if(i==1){
          wghts[i] = sizes[i]/2
        } else if(i==length(quant)){
          wghts[i] = sizes[i-1]/2
        } else {
          wghts[i] = sizes[i]/2 + sizes[i-1]/2
        }
      }
      
      # get the quantiles from the data
      #quantUsed = quantile(Data$prior_estimate, quant[-c(length(quant))])
      quantiles = stats::quantile(Data$mu_prior_estimate, quant)
      
      # get the mean prior variance from the data
      sigma2 = mean(Data$mu_prior_std_error**2)
      
      # apply the function to each BF, using the parameters estimated just before
      d$my_p = apply(d, 1, function(x) get_OneP(x[1], quantiles, sigma2, wghts))
      
      tmp = "... checking p-values near significance threshold  \n"
      Log = update_log(Log, tmp, verbose)
      
      # check if we have some false positive / false negative
      count_corrected=0
      # Identify last non-significant SNP
      snp = min(which(d$my_p>sign_thr))
      p_true = get_fullP(d$BF[snp], Data)
      NeedToCorrect = T
      while(NeedToCorrect){
        if(p_true<sign_thr){
          d$my_p[snp] = p_true
          count_corrected = count_corrected+1
          snp=snp+1
          p_true = get_fullP(d$BF[snp], Data)
        } else {
          NeedToCorrect=F
        }
      }
      
      
      # Identify first significant SNP
      snp = max(which(d$my_p<sign_thr))
      p_true = get_fullP(d$BF[snp], Data)
      NeedToCorrect = T
      while(NeedToCorrect){
        if(p_true>sign_thr){
          d$my_p[snp] = p_true
          count_corrected = count_corrected+1
          snp=snp-1
          p_true = get_fullP(d$BF[snp], Data)
        } else {
          NeedToCorrect=F
        }
      }
      
      if(count_corrected==0) tmp = "    everything is ok!  \n"
      if(count_corrected==1)  tmp = paste0("   ", count_corrected," p-value has been re-estimated using the exact formula.  \n")
      if(count_corrected>1)  tmp = paste0("   ", count_corrected, " p-values have been re-estimated using the exact formula.  \n")
      
      Log = update_log(Log, tmp, verbose)
      
      res = list()
      res$pValue = d$my_p
      res$log_info = Log
      
      return(res)
      
      
    }
    
    
    PVal = get_pValue(Prior$BF, Prior, sign_thresh, verbose)
    Prior %>%
      mutate(BF_p = PVal$pValue) -> Prior
    Log=c(Log, PVal$log_info)
    
  }
  
  
  if(sign_method == "fdr"){
    tmp = "# Estimating FDR (Benjamini-Hochberg procedure)... \n"
    Log = update_log(Log, tmp, verbose)
    
    Prior %>%
      mutate(BF_fdr = stats::p.adjust(.data$BF_p)) -> Prior
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  
  ### add p-values for posterior / direct effects
  tmp = "# Estimating p-values for posterior effects... \n"
  Log = update_log(Log, tmp, verbose)
  
  Prior %>%
    mutate(p_posterior = 2 * stats::pnorm(-abs(.data$z_posterior))) -> Prior
  
  
  if(sign_method == "fdr"){
    tmp = "# Estimating FDR (Benjamini-Hochberg procedure) for posterior effects... \n"
    Log = update_log(Log, tmp, verbose)
    
    Prior %>%
      mutate(fdr_posterior = stats::p.adjust(.data$p_posterior)) -> Prior
    
    
  }
  
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  tmp = "# Estimating p-values for direct effects... \n"
  Log = update_log(Log, tmp, verbose)
  
  Prior %>%
    mutate(p_direct = 2 * stats::pnorm(-abs(.data$z_direct))) -> Prior
  
  
  if(sign_method == "fdr"){
    tmp = "# Estimating FDR (Benjamini-Hochberg procedure) for direct effects... \n"
    Log = update_log(Log, tmp, verbose)
    
    Prior %>%
      mutate(fdr_direct = stats::p.adjust(.data$p_direct)) -> Prior
    
    
  }
  
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  
  
  if(save_files){
    system("rm Prior.csv")
    readr::write_csv(Prior, "PriorBFp.csv")
    tmp = "The file Prior.csv has been updated into Prior_BFp.csv \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  
  
  res=list(log_info = Log,
           SNPs = Prior)
  return(res)
}
