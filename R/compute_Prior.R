###### Function to compute prior for each SNP######



# #' Compute Prior
# #' 
# #' From a list of significant studies (from \code{identify_StudiesMR()}), the pruned Z-matrix
# #' of MR instruments (from \code{makeMR_ZMatrix()}) and the full Z-matrix of all SNPs for these
# #' significant studies (from \code{makeFull_ZMatrix()}), compute the prior
# #' 
# #' @inheritParams bGWAS
# #' @param selected_studies list of significant prior GWASs file names (character)
# #' @param MR_ZMatrix, used to get per-chromosome estimate (data.frame)
# #' @param All_ZMatrix, used to get prior effects (data.frame)
# #' @param GWASData, used to add prior/posterior/direct to results (data.frame)
# #' @param rescaling, should we get the effects on beta scale? (logical)
# NOT EXPORTED


compute_prior <- function(selected_studies, MR_ZMatrix, All_ZMatrix, GWASData, rescaling, MR_shrinkage, prior_shrinkage, Z_matrices, save_files=FALSE, verbose=FALSE){
  
  Log = c()
 
  
  # Function to automatically generate the formula for linear model
  generate_Formula <- function(outcome, study_names, with_intercept=F ) {
    formula = ifelse(with_intercept,
                     paste(outcome, ' ~  1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))),
                     paste(outcome, ' ~ -1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))))
    return(formula)
  }
  
  MR_ZMatrix %>%
    names() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    slice(ncol(MR_ZMatrix)) %>%
    pull() -> my_outcome
  
  generate_Formula(my_outcome, selected_studies) -> myFormula
  
  all.priors =  tibble()
  all.coefs = tibble()
  
  tmp = "# Calculating the prior chromosome by chromosome... \n"
  Log = update_log(Log, tmp, verbose)
  
  OutOfSample = tibble(Obs=numeric(), Pred=numeric(), 
                           Res=numeric())
  
  for(chr in 1:22) {
    tmp = paste0("   Chromosome ", chr, "\n")
    Log = update_log(Log, tmp, verbose)
    
    if(is.na(table(All_ZMatrix$chrm)[chr])){
      tmp = "No SNP on this chromosome \n"
      Log = update_log(Log, tmp, verbose)
      
      # ALSO TEST IF NSNP < NSTUDIES?
      next
    }
    
    # create the dataset without this chromosome - for MR instruments
    MR_ZMatrix %>%
      filter(.data$chrm != chr) -> MR_ZMatrix_masked
    
    tmp = "Running regression, \n"
    Log = update_log(Log, tmp, verbose)
    
    
    stats::lm(data=MR_ZMatrix_masked, formula = myFormula
    ) -> fit_masked # fit, without one chromosome
    fit_masked %>%  
      summary %>%
      stats::coef() %>%
      as.data.frame() %>% 
      tibble::rownames_to_column("nm") %>%
      mutate(chrm=chr) -> coefs
    
    # stopifnot(length(selected_studies) ==  nrow(coefs))
    ## this can happen, if all non-0 Z-score of a study are one the same chromosome
    ## but this should not be a problem
    # just add some warning message here :
    if(length(selected_studies) !=  nrow(coefs)){
      missingSt = selected_studies[!selected_studies %in% rownames(coefs)]
      tmp = paste0("No effect estimate when masking this chromosome for : ",
                   paste0(missingSt, collapse=" - "),
                   " (all SNPs having non-0 Z-scores are on this chromosome) \n")
      Log = update_log(Log, tmp, verbose)
      
    }
    
    # Do we need it?
    coefs %>%
      arrange(desc(.data$`Pr(>|t|)`)) -> coefs
    
    stopifnot(nrow(coefs) >= 1) # why would this happen?
    # if only one study selected, and all instruments on the same chromosome?
    # better handle this!
    
    all.coefs %>% 
      bind_rows(coefs) -> all.coefs
    
    # check the predictions using data from the training set
    #in_train = data.frame( CRG=d_masked[,..outcome], prd=predict.lm(fit_masked) )
    #residual.variance.in.training=var(in_train[,1]-in_train$prd)
    
    tmp = "Calculating prior estimates for SNPs on this chromosome, \n"
    Log = update_log(Log, tmp, verbose)
    
    # all SNPs we need to predict (the ones on this chromosome)
    All_ZMatrix %>%
      filter(.data$chrm==chr) -> d_test
    
    suppressWarnings({ #In predict.lm(fit_masked, d_test, se.fit = T) : prediction from a rank-deficient fit may be misleading
      stats::predict(fit_masked, d_test, se.fit=T) -> preds_test
    })
    
    tmp = "Calculating prior standard errors for SNPs on this chromosome, \n"
    Log = update_log(Log, tmp, verbose)
    
    
    #A few lines to compute the extra.variance if we allow that the z's are random too
    cv <- stats::vcov(fit_masked)
    extra.variance.1 <- sum(diag(cv)) * length(selected_studies)
    
    summary(fit_masked) %>%
      stats::coef() %>%
      as.data.frame() %>%
      pull(.data$Estimate) -> Estimates
    
    t(Estimates) %*% Estimates -> extra.variance.2
    
    extra.variance <- extra.variance.1 + extra.variance.2
    
    d_test %>%
      select(c(1:5, ncol(d_test))) %>%
      mutate(fit=preds_test$fit,
             se = sqrt(preds_test$se.fit^2 + c(extra.variance))) -> d_test

    
    d_test %>%
      bind_rows(all.priors, .data) -> all.priors
    
  
    # subset - MR instruments only
    
    d_test %>%
      filter(.data$rs %in% MR_ZMatrix$rs) %>%
      select(Obs=6,
             Pred=7) %>%
      mutate(Res=.data$Obs-.data$Pred) %>%
      bind_rows(OutOfSample, .data) -> OutOfSample
    
  }
  
  # Out of sample R2
  SS.total     <- sum((OutOfSample$Obs-mean(OutOfSample$Obs))^2)
  SS.residuals <- sum((OutOfSample$Res)^2)
  
  R2 <- 1 - SS.residuals/SS.total
  
  tmp = paste0("## Out-of-sample R-squared for MR instruments across all chromosomes is ", round(R2, 4), "\n")
  Log = update_log(Log, tmp, verbose)
  tmp = paste0("## Out-of-sample squared correlation for MR instruments across all chromosome is ", round(stats::cor(OutOfSample$Obs, OutOfSample$Pred)^2, 4), "\n")
  Log = update_log(Log, tmp, verbose)
  
  
  tmp = paste0("## Correlation between prior and observed effects for all SNPs is ", round(stats::cor(all.priors[,6:7])[1,2], 4), "\n")
  Log = update_log(Log, tmp, verbose)
  
  Zlimit <- stats::qnorm(0.001/2, lower.tail = F)
  all.priors %>%
    filter_at(vars(6), any_vars(abs(.data$.)>Zlimit)) -> SNPs_moderateeffect
  tmp = paste0("## Correlation between prior and observed effects for SNPs with GWAS p-value < 0.001 is ", round(stats::cor(SNPs_moderateeffect[,6:7])[1,2], 4), "\n")
  Log = update_log(Log, tmp, verbose)
  
  
  # Merge with GWASData, to keep all columns, nicely aligned
  # keep the SNPs in our Z matrix and order them correctly
  GWASData %>%
    slice(match(all.priors$rs, .data$rsid)) -> GWASData
  
  GWASData %>%
    mutate(mu_prior_estimate=case_when(
      .data$z_obs == all.priors[,my_outcome] ~ all.priors$fit, # aligned
      TRUE ~ - all.priors$fit),  #swapped
      mu_prior_std_error= all.priors$se) -> Results
  

  
  ## also add the posterior / direct
  Results %>%
    mutate(mu_posterior_estimate = (.data$mu_prior_std_error**2/
                                      (.data$mu_prior_std_error**2+1)) *
             ((.data$mu_prior_estimate/.data$mu_prior_std_error**2) + .data$z_obs),
           mu_posterior_std_error = sqrt(.data$mu_prior_std_error**2/
                                           (.data$mu_prior_std_error**2+1)),
           z_posterior = .data$mu_posterior_estimate/.data$mu_posterior_std_error,
           mu_direct_estimate = (.data$z_obs - .data$mu_prior_estimate),
           mu_direct_std_error = sqrt(1 + .data$mu_prior_std_error**2),
           z_direct = .data$mu_direct_estimate/.data$mu_direct_std_error) -> Results

  
  # if rescaling
  if(rescaling){
    # prior
    Results %>%
      mutate(beta_prior_estimate = .data$se * .data$mu_prior_estimate,
             beta_prior_std_error = .data$se * .data$mu_prior_std_error) -> Results
    
    # posterior
      Results %>%
        mutate(beta_posterior_estimate = .data$se * .data$mu_posterior_estimate,
               beta_posterior_std_error = .data$se * .data$mu_posterior_std_error) -> Results
    # direct
    Results %>%
      mutate(beta_direct_estimate = .data$se * .data$mu_direct_estimate,
             beta_direct_std_error = .data$se * .data$mu_direct_std_error) -> Results
      
    # CRR
    Results %>%
      mutate(CRR = .data$beta_direct_estimate / .data$beta) -> Results
    
  } else {
    Results %>%
      mutate(CRR = .data$mu_direct_estimate / .data$z_obs) -> Results
    
  }
  
  # add UK10K chr/pos data.frame for SNPs that are kept
  All_ZMatrix %>%
    slice(match(Results$rsid, .data$rs)) %>%
    transmute(chrm_UK10K=.data$chrm, 
           pos_UK10K=.data$pos) -> ChrPos
  bind_cols(ChrPos, Results) -> Results
  
  
  
  all.coefs %>%
    arrange(.data$nm) %>%
    set_names(c("study", "estimate", "std_error", "T", "P", "chrm")) -> all.coefs
  
  if(save_files){
    readr::write_csv(path="CoefficientsByChromosome.csv", x=all.coefs)
    tmp = paste0("The file ", "CoefficientsByChromosome.csv has been successfully written. \n")
    Log = update_log(Log, tmp, verbose)
    
  }
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  
  
  
  if(save_files){
    readr::write_csv(path="Prior.csv", x=Results)
    tmp = paste0("The file ", "Prior.csv has been successfully written. \n")
    Log = update_log(Log, tmp, verbose)
    
  }
  
  res=list(log_info = Log,
           prior = Results,
           all_coeffs = all.coefs) 
  return(res)
}
 
