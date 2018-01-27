###### Function to calculate BF and p-values ######



#' Calculate Bayes Factor and empirical p-values from observed Z-scores and priors
#'
#' From a pruned Z-matrix of strong MR instruments, performs multivariate MR
#' @inheritParams bGWAS
#' @param Prior from the function makePrior()
#
#' @importFrom magrittr "%>%"
#'
#' @return Log + data.table containing rs-chr-pos-alt-ref-obs-fit-se-...
#'
# Function not exported, no need of extended documentation?



request_BFandP <- function(Prior, saveFiles=F, verbose=F) {
  Log = c()
  tmp = paste0("# Computing observed Bayes Factor for all SNPs... \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  # The true one : calculate BF
  Prior$log.bf =  dnorm(log=T, mean= Prior$fit    , sd=sqrt(
    1+Prior$se**2
  ) , x=Prior$obs) -
    dnorm(log=T, mean= 0.0       , sd=sqrt(
      1
    ) , x=Prior$obs)
  # order by BF
  Prior = Prior[order(-Prior$log.bf),]


  tmp = "Done! \n"
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
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
  	Log = c(Log, tmp)
  	if(verbose) cat(tmp)

  	NUMBER_OF_NULLRUNS = 100
  	full.number.of.nulls.for.comparison = NUMBER_OF_NULLRUNS* N.one.set.of.nulls

  	tmp = paste0(NUMBER_OF_NULLRUNS, " random z-scores simulations for the whole set of SNPs (",
  	             format(nrow(Prior), big.mark = ",", scientific = F), ") : ", format(full.number.of.nulls.for.comparison, big.mark = ",", scientific = F), " null BFs are calculated \n")
  	Log = c(Log, tmp)
  	if(verbose) cat(tmp)
  	for( i in 1:NUMBER_OF_NULLRUNS ) {
  	  # null bf
  	  if(i %% 20 == 1){
  	    tmp = paste0("Simulation ",i, "... \n")
  	    Log = c(Log, tmp)
  	    if(verbose) cat(tmp)
  	  }

  	  Prior$z <-   rnorm(nrow(Prior), 0, 1)

  	  Prior$log.bf.null =  dnorm(log=T, mean= Prior$fit    , sd=sqrt(
  	    1+Prior$se**2
  	  ) , x=Prior$z) -
  	    dnorm(log=T, mean= 0.0       , sd=sqrt(
  	      1
  	    ) , x=Prior$z)
  #	  Prior = Prior[!is.na(Prior$log.bf.null),]
  #	  Prior = Prior[order(-Prior$log.bf.null),]
      TempP = Prior[order(-Prior$log.bf.null),]
  	  nullbf = TempP$log.bf.nul

  	  N = length(nullbf)

  	  stopifnot(N == N.one.set.of.nulls)
  	  #(!is.unsorted(-nullbf)) %|%stopifnot

  	  count_larger(Prior$log.bf, nullbf) -> counts
  	  stopifnot(length(counts) == nrow(Prior))
  	  running.total = running.total + counts
  	}

  	tmp = "Done! \n"
  	Log = c(Log, tmp)
  	if(verbose) cat(tmp)

  	tmp = "# Computing empirical p-values... \n"
  	Log = c(Log, tmp)
  	if(verbose) cat(tmp)

  	Prior[, z := NULL]
  	Prior[, log.bf.null := NULL]
  	Prior[, nullls_count := running.total]

  	#stopifnot(max(running.total) <= full.number.of.nulls.for.comparison)
  	Prior[, BF_p := (running.total+1) / full.number.of.nulls.for.comparison ]


  	if(saveFiles){
      system("rm Prior.csv")
  	  readr::write_csv(Prior, "PriorBFp.csv")
  	  tmp = "The file Prior.csv has been updated into Prior_BFp.csv \n"
  	  Log = c(Log, tmp)
  	  if(verbose) cat(tmp)
  	}

  	tmp = "Done! \n"
  	Log = c(Log, tmp)
  	if(verbose) cat(tmp)


  	res=list()
  	res$Log = Log
  	res$SNPs = Prior
  	return(res)
}
