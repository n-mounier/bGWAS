###### Function to calculate BF and p-values ######



# #' Calculate Bayes Factor and empirical p-values from observed effects and priors
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED




request_BFandP <- function(Prior, parent0, n_permutations, save_files=F, verbose=F) {
  Log = c()
  tmp = paste0("# Computing observed Bayes Factor for all SNPs... \n")
  Log = update_log(Log, tmp, verbose)
  
  
  
  # calculate BF
  Prior %>%
    mutate(BF = stats::dnorm(mean= .data$mu_prior_estimate, 
                             sd= sqrt(.data$mu_prior_std_error^2+.data$std_se^2), 
                             x = .data$std_beta) /
             stats::dnorm(mean = 0.0, 
                   sd=.data$std_se, 
                   x = .data$std_beta)) %>%  # and order by BF
    arrange(desc(.data$BF)) -> Prior
  
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  
  tmp = "# Computing BF p-values... \n"
  Log = update_log(Log, tmp, verbose)
  
  tmp = "   using a permutation approach: \n"
  
  
  # how many null z-scores should we simulate?
  nrow(Prior) -> N.one.set.of.nulls
  
  # Function to count how many larger BF were observed in the simulation
  count_larger <- Rcpp::cppFunction(' IntegerVector count_larger(NumericVector real, NumericVector null) {
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
  } ', env=NULL)
  
  
  
  # counts for each SNPs : start with all 0
  running.total = rep.int(0 , nrow(Prior))
  
  tmp = "# Computing null Bayes Factors (needed to obtain empirical p-values)... \n"
  Log = update_log(Log, tmp, verbose)
  
  full.number.of.nulls.for.comparison = n_permutations* N.one.set.of.nulls
    
  tmp = paste0(n_permutations, " random standardised effects simulations for the whole set of SNPs (",
               format(nrow(Prior), big.mark = ",", scientific = F), ") : ", format(full.number.of.nulls.for.comparison, big.mark = ",", scientific = F), " null BFs are calculated \n")
  Log = update_log(Log, tmp, verbose)
  

  for( i in 1:n_permutations ) {
    # null bf
    if(i %% 20 == 1){
      tmp = paste0("Simulation ",i, "... \n")
      Log = update_log(Log, tmp, verbose)
      
    }
    
    # simulate z-score from normal distribution
    # ... but this does not account for the real "true distribution" (inflation?)
    # should we shuffle z-scores instead?
    # Prior %>%
    #   mutate(z = stats::rnorm(N.one.set.of.nulls, 0, 1),
    #          BF.null =  stats::dnorm(mean= .data$mu_prior_estimate, 
    #                               sd=sqrt(1+.data$mu_prior_std_error**2), 
    #                               x=.data$z) /
    #            stats::dnorm(mean= 0.0, 
    #                         sd=sqrt(1),
    #                         x=.data$z)) %>%
    #   arrange(desc(.data$BF.null))-> Prior_BF
    mySE = mean(Prior$std_se)
    Prior %>%
      mutate(std_beta = stats::rnorm(N.one.set.of.nulls, 0, 1/sqrt(.data$N)),
             BF.null =  stats::dnorm(mean= .data$mu_prior_estimate,
                                     sd=sqrt(.data$std_se^2+.data$mu_prior_std_error**2),
                                     x=.data$std_beta) /
               stats::dnorm(mean= 0.0,
                            sd=.data$std_se,
                            x=.data$std_beta)) -> Prior_BF

    Prior_BF %>%
      arrange(desc(.data$BF.null)) %>%
      pull(.data$BF.null) -> nullbf
    
    
    count_larger(Prior$BF, nullbf) -> counts
    stopifnot(length(counts) == nrow(Prior))
    running.total = running.total + counts
  }
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  
  Prior$z = NULL
  Prior$BF.null = NULL
  
  Prior$BF_p = (running.total+1) / full.number.of.nulls.for.comparison
  
  
  
  if(save_files){
    if(parent0=="mat0"){
      system("rm Prior_mat0.csv.gz")
      data.table::fwrite(Prior, "PriorBFp_mat0.csv.gz", sep = ",", compress="gzip")
      tmp = paste0("The file Prior_mat0.csv.gz has been updated into Prior_BFp_mat0.csv.gz. \n")
    } else {
      system("rm Prior_pat0.csv.gz")
      data.table::fwrite(Prior, "PriorBFp_pat0.csv.gz", sep = ",", compress="gzip")
      tmp = paste0("The file Prior_pat0.csv.gz has been updated into Prior_BFp_pat0.csv.gz. \n")
    }
    Log = update_log(Log, tmp, verbose)
    
  }
  
  
  
  res=list(log_info = Log,
           SNPs = Prior)
  return(res)
}
