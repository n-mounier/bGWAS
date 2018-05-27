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

    Prior = Prior[order(-Prior$BF),]




    # true Z = just rs + log.bf

    #trueZ[, !is.unsorted(-log.bf)] %|%stopifnot # should be sorted, why test again we just did it !!!

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


    Prior[, z := NULL]
    Prior[, BF.null := NULL]
   # Prior[, nullls_count := running.total]

    #stopifnot(max(running.total) <= full.number.of.nulls.for.comparison)
    Prior[, BF_P := (running.total+1) / full.number.of.nulls.for.comparison ]


  } else {

  ### Functions to estimate p-value from BF/prior distribution

  ## for zero prior (works with a vector of se2 -> returns a vector or p)
  # function to calculate "Chi2" value, and return the p-value
  # take variance : sigma**2 as input
  zero_snp <- function(se2, BF){
    p = pchisq(2 * (1+se2)/se2 * log(BF*sqrt(1+se2)), df=1, lower.tail=F)
    return(p)
  }

  ## for non-zero prior (works with vectors of mu/sigma2 -> returns a vector or p)
  # function to calculate x1 and x2, and return the p-value
  # for SNP with non-zero prior (mu, sigma**2) - observed BF
  nonzero_snp <- function(mu, sigma2, BF){
    p = numeric(length(mu))

    a = -1/(2*(sigma2+1)) + 1/2
    b = mu/(sigma2+1)
    c = -mu**2/(2*(sigma2+1)) - log(sqrt((sigma2+1))*BF)
    delta = b**2 - 4*a*c

    p[delta<0]=1

    x1 = (-b[delta>0] - sqrt(delta[delta>0]))/(2*a[delta>0])
    x2 = (-b [delta>0]+ sqrt(delta[delta>0]))/(2*a[delta>0])
    p1 = pnorm(x1, lower.tail = T)
    p2 = pnorm(x2, lower.tail = F)
    p[delta>0] = p1+p2

    return(p)
  }

  # calculate p-value for an observed BF, using formula from part 1 + part 2
  # BF can be a vector of BFs
  get_pValue <- function(BF, se_bychrm, zero_bychrm, non_zero){
    res=data.frame(BF=BF)
    res$my_p=NA
    ntotal = sum(zero_bychrm) + nrow(non_zero)

    # deal with all zero SNPs / deal with non-zero SNPs
    res$my_p = apply(res, 1, function(x) ((sum(zero_snp(se_bychrm**2, x[1]) * zero_bychrm)
                                           + sum(nonzero_snp(non_zero$prior_estimate, non_zero$prior_std_error**2, x[1])))) / ntotal)

    return(res$my_p)
  }


  # calculate p-value for an observed BF, using formula from part 1 only (zero prior SNPs)
  # BF can be a vector of BFs
  get_pValue_part1 <- function(BF, se_bychrm, zero_bychrm){
    res=data.frame(BF=BF)
    res$my_p=NA
    ntotal = sum(zero_bychrm)

    # deal with all zero SNPs
    res$my_p = apply(res, 1, function(x) sum(zero_snp(se_bychrm**2, x[1]) * zero_bychrm) / ntotal)

    return(res$my_p)
  }


  # calculate p-value for an observed BF, using formula from part 2 only (non-zero prior SNPs)
  # BF can be a vector of BFs
  get_pValue_part2 <- function(BF, non_zero){
    res=data.frame(BF=BF)
    res$my_p=NA
    ntotal =  nrow(non_zero)

    # deal with all zero SNPs / deal with non-zero SNPs
    res$my_p =  apply(res, 1, function(x) sum(nonzero_snp(non_zero$prior_estimate, non_zero$prior_std_error**2, x[1])) / ntotal)

    return(res$my_p)
  }


  ## Percentiles

  # funtion to approximate the second part of the formula, using percentiles
  # sigma_mu + the number of non-zero SNPs
  get_app_part2 <- function(BF, perc, sigma2, n){
    # here we need to suppress warnings because if BF small, delta<0, x1/x2 NA -> pb with p value
    suppressWarnings({
      a = pnorm( (-perc/sigma2) +
                   sqrt(1+sigma2)/sigma2 *
                   sqrt(perc**2 + 2*sigma2 * log(BF * sqrt(1+sigma2))), lower.tail=F)
      + pnorm((-perc/sigma2) -
                sqrt(1+sigma2)/sigma2 *
                sqrt(perc**2 + 2*sigma2 * log(BF * sqrt(1+sigma2))), lower.tail=T)
    })
    a[is.na(a)]=1
    p=n*mean(a)
    return(p)
  }


  get_pValue_percentiles <- function(BF, se_bychrm, zero_bychrm, non_zero){
    res=data.frame(BF=BF)
    res$my_p=NA
    ntotal = sum(zero_bychrm) + nrow(non_zero)
    n2 = nrow(non_zero)
    perc_mu = quantile(non_zero$prior_estimate, probs = seq(0,1,0.01))
    sigma2_mu = mean(non_zero$prior_std_error)

    # deal with all zero SNPs / deal with non-zero SNPs
    res$my_p = apply(res, 1, function(x) (sum(zero_snp(se_bychrm**2, x[1]) * zero_bychrm)
                                          + get_app_part2(x[1], perc_mu, sigma2_mu, n2)) / ntotal)

    return(res$my_p)
  }


  # calculate p-value for an observed BF, using an hybrid approach
  # BF can be a vector of BFs
  # for each window, use only part 1, then compare for the last SNP : part 1 vs full p value
  # if "high" difference, compare for the last SNP : percentile vs full p value
  # if ok, use the percentile approach for this window
  # if still "high" difference, re-calculate all p values from this window using part 1 + part 2
  # if any p of the window "close" to significance, also re-calculate
  get_pValue_hybrid_final <-  function(BF, se_bychrm, zero_bychrm, non_zero,
                                       size_window=1000, sign_thr=5e-8, near_sign=1, tolerance=0.05){

    n_window = ceiling(length(BF)/size_window)
    near_sign = -log10(near_sign * sign_thr)

    if(!is.unsorted(-BF)) stop("BF should be ordered decreasingly!")

    res=data.frame(BF=BF)
    res$new_p=NA

    # get the "zero-prior SNPs approximation"
    res$new_p =  get_pValue_part1(res$BF, se_bychrm, zero_bychrm)

    for(wind in 1:n_window){
      min_w = (wind-1)*size_window+1
      max_w = ifelse(wind!=n_window,  wind*size_window, nrow(res) )


      # are we close to our significance threshold ?
      # if true, we should use the full formula to avoid false positive / false negative
      if(any(-log10(abs(res$new_p[min_w:max_w]-sign_thr))>near_sign)){
        res$new_p[min_w:max_w] = get_pValue(res$BF[min_w:max_w], se_bychrm, zero_bychrm, non_zero)
      } else { # sanity check : is our approximation "too" different from the true value ?
        # true p-value, standard to compare with
        p2 = get_pValue(res$BF[max_w], se_bychrm, zero_bychrm, non_zero)
        # tolerance for this window
        my_tolerance = -log10(tolerance*p2)
        if(-log10(abs(res$new_p[max_w]-p2))<my_tolerance){ # if true
          # compare the "percentile approximation
          p_perc = get_pValue_percentiles(res$BF[max_w], se_bychrm, zero_bychrm, non_zero)
          # is our approximation "too" different from the true value ?
          if(-log10(abs(p_perc-p2))<my_tolerance){ # if true
            res$new_p[min_w:max_w] = get_pValue(res$BF[min_w:max_w], se_bychrm, zero_bychrm, non_zero)
          } else { # if false
            # use the "percentile approximation" to get the p-values for this window
            res$new_p[min_w:max_w] = get_pValue_percentiles(res$BF[min_w:max_w], se_bychrm, zero_bychrm, non_zero)
          }
        }
      }
    }
    return(res$new_p)
  }


  non_zero = Prior[Prior$prior_estimate!=0,]
  zero = Prior[Prior$prior_estimate==0,]
  zero_bychrm = table(zero$chrm)
  se_bychrm = numeric(length = 22)
  names(se_bychrm) = c(1:22)
  se_bychrm = zero$prior_std_error[match(names(se_bychrm), zero$chrm)]



  Prior$BF_P = get_pValue_hybrid_final(Prior$BF,
                                       se_bychrm = se_bychrm, zero_bychrm = zero_bychrm, non_zero = non_zero,
                                       size_window = 1000, sign_thr = sign_thresh, tolerance=0.05, near_sign = 1)

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
