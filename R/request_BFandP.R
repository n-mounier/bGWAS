###### Function to calculate BF and p-values ######



# #' Calculate Bayes Factor and empirical p-values from observed Z-scores and priors
# #'
# #' From a pruned Z-matrix of strong MR instruments, performs multivariate MR
# #' @inheritParams bGWAS
# #' @param Prior from the function makePrior()
#' @importFrom magrittr "%>%"
# #'
# #' @return Log + data.table containing rs-chr-pos-alt-ref-obs-fit-se-...
# #'
# Function not exported, no need of extended documentation?



request_BFandP <- function(Prior, sign_thresh, use_perm = F, save_files=F, verbose=F) {
  Log = c()
  tmp = paste0("# Computing observed Bayes Factor for all SNPs... \n")
  Log = update_log(Log, tmp, verbose)



  # calculate BF
  Prior$BF =  dnorm(mean= Prior$prior_estimate    , sd=sqrt(
    1+Prior$prior_std_error**2
  ) , x=Prior$observed_Z) /
    dnorm(mean= 0.0       , sd=sqrt(
      1
    ) , x=Prior$observed_Z)
  # order by BF
  Prior = Prior[order(Prior$BF),]



  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)

  # true Z = just rs + log.bf



  tmp = "# Computing BF p-values... \n"
  Log = update_log(Log, tmp, verbose)

  if(use_perm){# The true one : calculate BF
    tmp = "   using a permutation approach: \n"
    Log = update_log(Log, tmp, verbose)



    Prior = Prior[order(-Prior$BF),]


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

      Prior$z <-   rnorm(nrow(Prior), 0, 1)

      Prior$BF.null =  dnorm(mean= Prior$prior_estimate    , sd=sqrt(
        1+Prior$prior_std_error**2
      ) , x=Prior$z) /
        dnorm(mean= 0.0       , sd=sqrt(
          1
        ) , x=Prior$z)

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
    ## Calculate a p-value using the percentiles or any quantiles
    get_PercentileP<- function(BF, perc, sigma2){ # BF : value -  perc : vector of quantiles - sigma2 : (variance) value
      # here we need to suppress warnings because if BF small, delta<0, x1/x2 NA -> pb with p value
      suppressWarnings({
        a = pnorm( (-perc/sigma2) +
                     sqrt(1+sigma2)/sigma2 *
                     sqrt(perc**2 + 2*sigma2 * log(BF * sqrt(1+sigma2))), lower.tail=F) +
          pnorm((-perc/sigma2) -
                  sqrt(1+sigma2)/sigma2 *
                  sqrt(perc**2 + 2*sigma2 * log(BF * sqrt(1+sigma2))), lower.tail=T)
      })
      a[is.na(a)]=1 # if NA -> replace by 1, that means that the polynome has no real solution and Pr(allBF>BF)=1
      p=mean(a) # here, we assume that the quantiles are linearly distributed, so we can use an unweighted mean
      return(p)
    }

    # Get all p-values using the percentiles or any quantiles (wrapper function for "apply")
    get_PercentileAll<- function(BF, perc, sigma2){ # BF : vector of values value -  perc : vector of quantiles - sigma2 : (variance) value
      res=data.frame(BF=BF)
      res$my_p=NA

      res$my_p = apply(res, 1, function(x) get_PercentileP(x[1], perc, sigma2))

      return(res$my_p)
    }


    ## Calculate a p-value using all priors
    get_fullP <- function(BF, Data){ # BF : value - Data : data.frame with prior_estimate and prior_std_error columns
      suppressWarnings({
        a = pnorm( (-Data$prior_estimate/(Data$prior_std_error**2)) +
                     sqrt(1+Data$prior_std_error**2)/(Data$prior_std_error**2) *
                     sqrt(Data$prior_estimate**2 + 2*Data$prior_std_error**2 * log(BF * sqrt(1+Data$prior_std_error**2))), lower.tail=F) +
          pnorm((-Data$prior_estimate/(Data$prior_std_error**2)) -
                  sqrt(1+Data$prior_std_error**2)/(Data$prior_std_error**2) *
                  sqrt(Data$prior_estimate**2 + 2*Data$prior_std_error**2 * log(BF * sqrt(1+Data$prior_std_error**2))), lower.tail=T)
      })
      a[is.na(a)]=1 # if NA -> replace by 1, that means that the polynome has no real solution and Pr(allBF>BF)=1
      p=mean(a) # here, we use all the prior values, so we can use an unweighted mean
      return(p)
    }

    # Get all p-values using all priors (wrapper function for "apply")
    get_fullAll<- function(BF, Data){ # BF : vector of values - Data : data.frame with prior_estimate and prior_std_error columns
      res=data.frame(BF=BF)
      res$my_p=NA

      res$my_p = apply(res, 1, function(x) get_fullP(x[1], Data))

      return(res$my_p)
    }


    # Get all p-values using an hybrid approach
    get_pValue <-  function(BF, Data,
                            size_window=100, sign_thr=5e-8){

      n_window = ceiling(length(BF)/size_window)

      if(is.unsorted(rev(BF))) stop("BF should be ordered decreasingly!")

      res=data.frame(BF=BF)
      res$new_p=NA

      # get the "percentile approximation"
      perc = quantile(Data$prior_estimate, probs = seq(0,1,0.01))
      sigma2 = mean(Data$prior_std_error**2)

      res$new_p =  get_PercentileAll(res$BF, perc, sigma2)

      ### Different steps, for different versions ####
      # use "full p near significance" (only for the delicate window)
      # + more quantiles for smaller p (identify the number of quantiles needed to keep the right order)
      start_quantile = 0
      for(wind in 1:n_window){
        print(wind)
        min_w = (wind-1)*size_window+1
        max_w = ifelse(wind!=n_window,  wind*size_window, nrow(res) )

        # each time, test the first one and the last one
        pmin = get_fullP(res$BF[min_w], Data)
        if(wind ==1 && pmin>sign_thr){
          break # no SNP significant
        }
        # need to check if the limit was between two windows
        if(wind != 1 && pmin>sign_thr){ # all the significant SNPs were in the previous windows
          start_quantile = min_w-1
          # + check the previous window, if perc approach was ok or if we need to reestimate it
          wind=wind-1
          min_w = (wind-1)*size_window+1
          max_w = ifelse(wind!=n_window,  wind*size_window, nrow(res) )
          if(res$new_p[min_w]>sign_thr){
            start_quantile = min_w-1
            res$new_p[min_w:max_w] = get_fullAll(res[min_w:max_w,"BF"], Data)
          }
          break
        }
        pmax =  get_fullP(res$BF[max_w], Data)
        # 1 if all sign / if all non sign, 0 if discordant
        trueP = sign((pmin-sign_thr) * (pmax-sign_thr))
        if(trueP<1){
          start_quantile=min_w-1
          res$new_p[min_w:max_w] = get_fullAll(res[min_w:max_w,"BF"], Data)
          break
        }
      }
      if(start_quantile>0){
        try_more_quantiles = T
        p_larger = res$new_p[start_quantile+1]
        nquant = 5000 # by 1000
        while(try_more_quantiles == T){
          nquant = nquant + 1000
          print(nquant)
          quant = quantile(Data$prior_estimate, probs = seq(0,1,1/nquant))
          p_try = get_PercentileP(res$BF[start_quantile], quant, sigma2)
          if(p_larger > p_try) try_more_quantiles = F
        }
        res$new_p[1:start_quantile] = get_PercentileAll(res[1:start_quantile,"BF"], quant, sigma2)
      }
      return(res$new_p)
    }

    Prior = Prior[order(-Prior$BF),]
    Prior$BF_P = get_pValue(Prior$BF, Prior,
                            size_window = 100, sign_thr = sign_thresh)

  }

  if(save_files){
    system("rm Prior.csv")
    readr::write_csv(Prior, "PriorBFp.csv")
    tmp = "The file Prior.csv has been updated into Prior_BFp.csv \n"
    Log = update_log(Log, tmp, verbose)

  }
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)



  res=list()
  res$log_info = Log
  res$SNPs = Prior
  return(res)
}
