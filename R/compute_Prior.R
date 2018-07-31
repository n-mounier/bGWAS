###### Function to compute prior for each SNP######



# #' Compute Prior
# #'
# #' From a list of significant studies (from \code{identify_StudiesMR()}), the pruned Z-matrix
# #' of MR instruments (from \code{makeMR_ZMatrix()}) and the full Z-matrix of all SNPs for these
# #' significant studies (from \code{makeFull_ZMatrix()}), compute the prior
# #'
# #' @inheritParams bGWAS
# #' @param selected_studies data table
# #' @param MR_ZMatrix data table
# #' @param All_ZMatrix data table
# #'
#  Function not exported, no need of extended documentation?


compute_prior <- function(selected_studies, MR_ZMatrix, All_ZMatrix, MR_shrinkage, prior_shrinkage, Z_matrices, save_files=FALSE, verbose=FALSE){

  Log = c()


  tmp = paste0("# Preparation of the data... \n")
  Log = update_log(Log, tmp, verbose)


  push.extreme.zs.back.a.little.towards.zero <- function(d) { # Some z-scores are just too far from zero
    max.allowed.z = abs(stats::qnorm(1e-300 / 2)) # p=1e-300 is the max allowed now, truncate z-scores accordingly
    studies.here = names(d)
    studies.here = studies.here[!studies.here %in% c("rs","chrm","pos","alt","ref")]
    for(n in studies.here) {
      d[[n]] [ d[[n]] < -max.allowed.z ] <- -max.allowed.z
      d[[n]] [ d[[n]] >  max.allowed.z ] <-  max.allowed.z
    }
    d
  }

  MR_ZMatrix= push.extreme.zs.back.a.little.towards.zero(MR_ZMatrix)
  All_ZMatrix = push.extreme.zs.back.a.little.towards.zero(All_ZMatrix)
  # Set the z-scores to 0 for the regression if too low
  if(MR_shrinkage < 1.0) {
    for(column_of_zs in selected_studies) {
      print(column_of_zs)
      threshold = abs(stats::qnorm(MR_shrinkage/2))
      MR_ZMatrix[c(abs(MR_ZMatrix[,column_of_zs]) < threshold) , column_of_zs] <- 0
    }
    # check how many chromosomes have non-zero zs
    # is that really interesting??
   # data.table::rbindlist(lapply(selected_studies, function(column_of_zs) {
   #   MR_ZMatrix[ , .(study_name = column_of_zs, nz=sum(.SD[[column_of_zs]] != 0)) , by=chrm]
    #})) -> counts.by.chrm
   # counts.by.chrm[,.(chrms.with.nz=sum(nz!=0)), by=study_name] [order(-chrms.with.nz)] -> chrms.with.nz
#    if(save_files){
#      readr::write_csv(chrms.with.nz, 'chrms.with.nz.csv')
#    }
  }

  if(prior_shrinkage < 1.0) {
    for(column_of_zs in selected_studies) {
      threshold = abs(stats::qnorm(prior_shrinkage/2))
      All_ZMatrix[abs(All_ZMatrix[,column_of_zs]) < threshold , column_of_zs] <- 0
    }
    tmp = paste0("Applying shrinkage (threshold = ", prior_shrinkage, ") before calculating the prior. \n")
    Log = update_log(Log, tmp, verbose)
  }

  if(prior_shrinkage<1){
    NonZero_ZMatrix = All_ZMatrix[apply(All_ZMatrix[,6:(ncol(All_ZMatrix)-1)], 1, function(x) any(unlist(abs(x)>threshold))),1:(ncol(All_ZMatrix)-1)]
  } else {
    NonZero_ZMatrix = All_ZMatrix
  }


  generate.formula <- function(outcome, study_names ) {
    paste(paste0(outcome,' ~ -1 + '),paste(collapse=' + ', paste(sep='','`',study_names,'`')))
  }

  outcome = colnames(MR_ZMatrix)[ncol(MR_ZMatrix)]
  generate.formula(outcome, selected_studies) -> form

  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)



  all.priors =  data.table::data.table()
  all.coefs = data.table::data.table()

  tmp = "# Calculating the prior chromosome by chromosome... \n"
  Log = update_log(Log, tmp, verbose)

  R2=numeric(22)
  dynamic.study.names = unlist(selected_studies)

  for(chrm in 1:22) {
    tmp = paste0("   Chromosome ", chrm, "\n")
    Log = update_log(Log, tmp, verbose)


    if(is.na(table(All_ZMatrix$chrm)[chrm])){
      tmp = "No SNP on this chromosome \n"
      Log = update_log(Log, tmp, verbose)

      # ALSO TEST IF NSNP < NSTUDIES

      next
    }

    # create the dataset without this chromosome - for MR instruments
    train    =   (MR_ZMatrix$chrm != chrm)
    d_masked = MR_ZMatrix[train,]

    tmp = "Running regression, \n"
    Log = update_log(Log, tmp, verbose)


    lm(data=d_masked, formula = generate.formula(outcome, dynamic.study.names)
    ) -> fit_masked # fit, without one chromosome
    coefs <- data.frame(coef(summary(fit_masked)))

   # stopifnot(length(dynamic.study.names) ==  nrow(coefs))
     ## this can happen, if all non-0 Z-score of a study are one the same chromosome
     ## but this should not be a problem
    # just add some warning message here :
    if(length(dynamic.study.names) !=  nrow(coefs)){
      missingSt = dynamic.study.names[!dynamic.study.names %in% rownames(coefs)]
      tmp = paste0("No effect estimate when masking this chromosome for : ",
                  paste0(missingSt, collapse=" - "),
                  " (all SNPs having non-0 Z-scores are on this chromosome) \n")
      Log = update_log(Log, tmp, verbose)

    }
    coefs = coefs[order(-coefs$Pr...t..),,drop=F]
  # stopifnot(length(dynamic.study.names) ==  nrow(coefs))

    stopifnot(nrow(coefs) >= 1)

    coefs.DT =  data.table::data.table(chrm, data.table::data.table(coefs, keep.rownames=T))
    data.table::setnames(coefs.DT, 'rn', 'study_name')
    coefs.DT[ , study_name := gsub('`','',study_name) ]

    all.coefs = rbind(all.coefs, coefs.DT)

    # check the predictions using data from the training set
    #in_train = data.frame( CRG=d_masked[,..outcome], prd=predict.lm(fit_masked) )
    #residual.variance.in.training=var(in_train[,1]-in_train$prd)


    chrm_ = chrm

    tmp = "Calculating prior estimates for SNPs on this chromosome, \n"
    Log = update_log(Log, tmp, verbose)


    # all SNPs we need to predict (the ones on this chromosome)
    d_test             = subset(All_ZMatrix, chrm==chrm_)

    suppressWarnings({ #In predict.lm(fit_masked, d_test, se.fit = T) : prediction from a rank-deficient fit may be misleading
      predict(fit_masked, d_test, se.fit=T) -> preds
    })

    tmp = "Calculating prior standard errors for SNPs on this chromosome, \n"
    Log = update_log(Log, tmp, verbose)


    #A few lines to compute the extra.variance if we allow that the z's are random too
    {
      cv = vcov(fit_masked)
      sum(diag(cv))            -> extra.variance.1

      Estimate = coef(summary(fit_masked))[,"Estimate"]
      names(Estimate) %>% {. == rownames(cv)} %>% stopifnot
      Estimate = t(t(Estimate))
      t(Estimate) %*% Estimate -> extra.variance.2

      extra.variance = c(extra.variance.1 + extra.variance.2)

      rownames(cv) %>% {gsub('`','',.)} -> rns

    }


    nice.table = d_test[,c(
      "rs"
      ,"chrm"
      ,"pos"
      ,"alt"
      ,"ref"
      ,outcome)]
    colnames(nice.table)[6] = "obs"
    nice.table$fit = preds$fit
    nice.table$se = preds$se.fit
    # add extra variance
    nice.table$se = sqrt( (nice.table$se^2) + extra.variance )


    all.priors = rbind(all.priors, nice.table)

    tmp = "Calculating out of sample prediction for SNPs on this chromosome, \n"
    Log = update_log(Log, tmp, verbose)

    outcome = colnames(MR_ZMatrix)[ncol(MR_ZMatrix)]


    # all SNPs we need to predict (the ones on this chromosome)
    d_test             = MR_ZMatrix[chrm==chrm_          ,,drop=T]

    suppressWarnings({ #In predict.lm(fit_masked, d_test, se.fit = T) : prediction from a rank-deficient fit may be misleading
      predict(fit_masked, d_test, se.fit=T) -> preds
    })


    test.outcome    <- d_test[,outcome]

    SS.total      <- sum((test.outcome - mean(test.outcome))^2)
    SS.regression <- sum((preds$fit - mean(test.outcome))^2)


    # fraction of variability explained by the model : SS.regression/SS.total
    R2_chrm =  1-(1-SS.regression/SS.total)*(nrow(d_masked)-1)/(nrow(d_masked)-length(dynamic.study.names)-1)
    tmp = paste0("Adjusted R-squared : ", round(R2_chrm, 4), "\n")
    Log = update_log(Log, tmp, verbose)

    R2[chrm] =R2_chrm

  }
  colnames(all.priors)[6:8] = c("observed_Z", "prior_estimate", "prior_std_error")

  tmp = paste0("## Mean out-of-sample adjusted R-squared across all chromosomes is ", round(mean(R2), 4), "\n")
  Log = update_log(Log, tmp, verbose)
  tmp = paste0("## Median out-of-sample adjusted R-squared across all chromosomes is ", round(median(R2), 4), "\n")
  Log = update_log(Log, tmp, verbose)

  ## we need to add one to the prior variance to account for the fact that we are
  ## predicting a noisy variable : observedZ ~ N(trueZ, 1)
  all.priors$prior_std_error = sqrt(all.priors$prior_std_error**2 + 1)


  ## also add the posterior
  all.priors$posterior_estimate = (all.priors$prior_std_error**2/(all.priors$prior_std_error**2+1)) *
                           ((all.priors$prior_estimate/all.priors$prior_std_error**2)+all.priors$observed_Z)
  all.priors$posterior_std_error = sqrt(all.priors$prior_std_error**2/(all.priors$prior_std_error**2+1))



  data.table::setkey(all.coefs, study_name, chrm)

  colnames(all.coefs) = c("chrm", "study", "estimate", "std_error", "T", "P")

  if(save_files){
    readr::write_csv(path="CoefficientsByChromosome.csv", x=all.coefs)
    tmp = paste0("The file ", "CoefficientsByChromosome.csv has been successfully written. \n")
    Log = update_log(Log, tmp, verbose)

  }

  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)




  if(save_files){
    readr::write_csv(path="Prior.csv", x=all.priors)
    tmp = paste0("The file ", "Prior.csv has been successfully written. \n")
    Log = update_log(Log, tmp, verbose)

  }

  res=list()
  res$log_info = Log
  res$prior = as.data.frame(all.priors)
  res$all_coeffs = as.data.frame(all.coefs)
  res$non_zero = NonZero_ZMatrix
  return(res)
}
