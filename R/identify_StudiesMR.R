###### Function to identifiy studies significantly affecting trait of interest - multivariate MR ######



# #' Identify studies for prior using multivariate MR
# #'
# #' From a pruned Z-matrix of strong MR instruments, performs multivariate MR
# #'
# #' @inheritParams bGWAS
# #' @param ZMatrix The pruned Z-Matrix of strong instrument or the path of the file
# #'        containing it (data table or character)
# #'
# #' @return Log file + list of studies significant + create files (list them) if save_files=T
# #'
 #' @importFrom magrittr "%>%"
 #' @importFrom magrittr "%<>%"
# #'



identify_studiesMR <- function(ZMatrix, MR_shrinkage, save_files=FALSE, verbose=FALSE, power=FALSE){

  Log = c()

  if(save_files) Files_Info = data.table::fread("PriorGWASs.tsv")


  # if p-value studies > MR_shrinkage , Z-Score is set to zero

  # In the pipeline, should always be a data.table,
  # but for other cases, allow the use of a file
  if(data.table::is.data.table(ZMatrix)){ # ok if it comes from the pipeline, but still need to check that the format is ok
    if(!all(c("rs", "chrm", "pos", "alt", "ref") %in% colnames(ZMatrix))) stop("There is something wrong with colnames of the matrix")
  } else if(is.character(ZMatrix)){
    # TO BE DONE

  } else {
    stop("ZMatrix provided not correct")
  }

  # function to truncate extreme z-scores
  push.extreme.zs.back.a.little.towards.zero <- function(d) { # Some z-scores are just too far from zero
    max.allowed.z = abs(qnorm(1e-300 / 2)) # p=1e-300 is the max allowed now, truncate z-scores accordingly
    studies.here = names(d)
    studies.here = studies.here[!studies.here %in% c("rs","chrm","pos","alt","ref")]
    for(n in studies.here) {
      #pp(minmax(d[[n]]),n)
      d[[n]] [ d[[n]] < -max.allowed.z ] <- -max.allowed.z
      d[[n]] [ d[[n]] >  max.allowed.z ] <-  max.allowed.z
      #pp(minmax(d[[n]]),n)
    }
    d
  }

  # function to ...
  "%..%" <- function(d,form){
    p = function(..){}
    body(p,envir=parent.frame()) <- substitute(form)
    do.call(p,list(substitute(d)), envir=parent.frame())
  }

  tmp = paste0("#Preparation of the MR analyses to identify significant studies... \n")
  Log = update_log(Log, tmp, verbose)

  # Truncate Z-scores
  ZMatrix= push.extreme.zs.back.a.little.towards.zero(ZMatrix)



  All_study_names    = names(ZMatrix) %..% ..[!.. %in% c('rs','chrm','pos','alt','ref')]

  Prior_study_names = All_study_names[1:length(All_study_names)-1]
  tmp = paste0("Studies tested : ", paste(Prior_study_names, collapse = " - "), "\n")
  Log = update_log(Log, tmp, verbose)

  tmp = paste0("Conventionnal GWAS of interest : ", All_study_names[length(All_study_names)], "\n")
  Log = update_log(Log, tmp, verbose)


  # Set the z-scores to 0 for the regression if too low
  if(MR_shrinkage < 1.0) {
    for(column_of_zs in Prior_study_names) {
      threshold = abs(qnorm(MR_shrinkage/2))
      ZMatrix[c(abs(ZMatrix[,..column_of_zs]) < threshold) , column_of_zs] <- 0
    }
    tmp = paste0("Applying shrinkage (threshold = ", MR_shrinkage, ") before performing MR. \n")
    Log = update_log(Log, tmp, verbose)
  }


  # Function to automatically generate the formula for linear model
  generate.formula <- function(outcome, study_names, with.intercept=F ) {
    if(with.intercept) {
      paste(outcome, ' ~  1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`')))
    } else {
      paste(outcome, ' ~ -1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`')))
    }
  }



  # Compute and save the univariate regressions (used to check directionnality in multivariate regression):
  tmp = paste0("# Univariate regressions for each trait... \n")
  Log = update_log(Log, tmp, verbose)


  uni.coefs.collection = data.table::data.table()
  for(one.study in Prior_study_names) {
    paste(All_study_names[length(All_study_names)], ' ~ -1 + `',one.study,'`', sep='', collapse="") -> uni.form
    lm(data=ZMatrix, formula = uni.form) -> uni.fit
    uni.coefs <- data.frame(coef(summary(uni.fit)))
    uni.coefs=cbind(nm=rownames(uni.coefs) %>% gsub("`","",.), uni.coefs)
    rownames(uni.coefs) <- NULL
    uni.coefs = uni.coefs[order(uni.coefs$Pr...t..),,drop=F]
    uni.coefs$adj.r.squared <- summary(uni.fit)$adj.r.squared
    uni.coefs$r.squared     <- summary(uni.fit)$r.squared
    uni.coefs.collection = rbind(uni.coefs.collection, uni.coefs)
  }
  data.table::setkey(uni.coefs.collection, nm)


  # if some of the r2 are too high, this is probably the same trait and you should exclude
  # it from the analysis
  if(any(uni.coefs.collection$r.squared>0.7)){
    too.high=which(uni.coefs.collection$r.squared>0.7)
    Names = uni.coefs.collection[too.high,"nm"]
    r2 = uni.coefs.collection[too.high,"r.squared"]
    S = list_priorGWASs()
    IDs = S$ID[match(Names$nm, S$File)]
    Ref = paste0(Names$nm, " (ID = ", IDs, " - r2 = ", r2$r.squared, ")")

    tmp = paste0("The study ", Ref,
                " has a really high r squared, it probably corresponds to the",
                " same trait than your conventionnal GWAS - please remove it",
                " before re-rerunning the analysis. \n")
    Log = update_log(Log, tmp, verbose)

    res=list()
    res$log_info = Log
    res$stop = T
    return(res)
  }

  colnames(uni.coefs.collection) = c("study", "estimate", "std_error", "T", "P",
                                 "adj_Rsquared", "Rsquared")

  if(save_files){ # add univariate coeffs
   Files_Info$uni_estimate = numeric()
   Files_Info$uni_estimate[match(uni.coefs.collection$study, Files_Info$File)] = uni.coefs.collection$estimate
   Files_Info$uni_std_error = numeric()
   Files_Info$uni_std_error[match(uni.coefs.collection$study, Files_Info$File)] = uni.coefs.collection$std_error
   Files_Info$uni_T = numeric()
   Files_Info$uni_T[match(uni.coefs.collection$study, Files_Info$File)] = unlist(uni.coefs.collection[,"T"])
   Files_Info$uni_P = numeric()
   Files_Info$uni_P[match(uni.coefs.collection$study, Files_Info$File)] = uni.coefs.collection$P
   Files_Info$uni_adj_Rsquared = numeric()
   Files_Info$uni_adj_Rsquared[match(uni.coefs.collection$study, Files_Info$File)] = uni.coefs.collection$adj_Rsquared
   Files_Info$uni_Rsquared = numeric()
   Files_Info$uni_Rsquared[match(uni.coefs.collection$study, Files_Info$File)] = uni.coefs.collection$Rsquared
  }


  tmp = paste0("Done! \n")
  Log = update_log(Log, tmp, verbose)



  tmp = paste0("# Stepwise AIC multivariate regression... \n")
  Log = update_log(Log, tmp, verbose)


  significant.studies = Prior_study_names # might be shrunk in the coming lines
  stopifnot(length(significant.studies) == length(Prior_study_names))

  initial.formula = paste0(All_study_names[length(All_study_names)],' ~ -1')
  k=2 # for AIC
  # step() looks in global environment for the data it needs,
  # trick : use the assign function !
  assign("ZMatrix", ZMatrix, envir = .GlobalEnv)

  lm(data=ZMatrix, formula = initial.formula) -> fit
  step(fit, scope = list(
    lower= initial.formula
    ,upper= generate.formula(All_study_names[length(All_study_names)], Prior_study_names)
  )
  ,trace = F
  , k = k
  ) -> fit.stepped


  coefs = data.frame( coef(summary(fit.stepped)) )
  coefs = coefs[order(coefs$Pr...t..),,drop=F]
  significant.studies = gsub("`","",rownames(coefs))



  post_AIC_studies = significant.studies

  if(save_files){ # add Status : AIC exclusion
    StudiesNotSelected =   Prior_study_names[!Prior_study_names %in%   post_AIC_studies ]
    Files_Info$status[Files_Info$File %in% StudiesNotSelected] = "Excluded during multivariate MR (AIC stepwise selection)"
  }

  if(length(post_AIC_studies)==0){
    tmp = "No study significant - Analysis failed"
    Log = update_log(Log, tmp, verbose)

    if(save_files)                write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )

    res=list()
    res$log_info = Log
    res$stop = T
    return(res)
  }


  if(length(post_AIC_studies)>1){
    tmp = paste0(length(post_AIC_studies), " studies are significant : ", paste0(post_AIC_studies, collapse = " - "), "\n")
    Log = update_log(Log, tmp, verbose)
  } else {
    tmp = paste0(length(post_AIC_studies), " study is significant : ", paste0(post_AIC_studies, collapse = " - "), "\n")
    Log = update_log(Log, tmp, verbose)
  }

  "%|%" <- function(x,y){
    do.call(y,list(substitute(x)),envir=parent.frame()) # just right
  }

  tmp = paste0("Done! \n")
  Log = update_log(Log, tmp, verbose)


  # We now have a set of studies from AIC,but we may want to prune them further
  tmp = paste0("#Further pruning of these studies... \n")
  Log = update_log(Log, tmp, verbose)

  # model
  lm(data=ZMatrix, formula = generate.formula(All_study_names[length(All_study_names)], significant.studies)) %>%summary %>%coef %>% data.table::data.table(keep.rownames=T) -> coefs
  data.table::setnames(coefs,"rn","nm")
  coefs[, nm := gsub("`","",nm)]
  na.fail(coefs)
  studies_to_remove = NULL
  reasons = NULL

  # Remove traits with multivariate p-value larger that 0.05
  coefs[order(-`Pr(>|t|)`)] -> largest.MR.p.value
  while(largest.MR.p.value$`Pr(>|t|)`[1] > 0.05) {
    studies_to_remove <- c(studies_to_remove, largest.MR.p.value$nm[1])
    reasons          <- c(reasons, paste0('MR pvalue too high: ', format(digits=3, largest.MR.p.value$`Pr(>|t|)`[1])))
    largest.MR.p.value = largest.MR.p.value[-1,]
    if(nrow(largest.MR.p.value)==0){
      Files_Info$status[Files_Info$File %in% studies_to_remove] =
        "Excluded during multivariate MR (p-value > 0.05)"

      tmp = "No study significant in multivariate analysis with p<0.05 - Analysis failed"
      Log = update_log(Log, tmp, verbose)

      if(save_files)                write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )


      res=list()
      res$log_info = Log
      res$stop = T
      return(res)
    }
  }

  if(!is.null(studies_to_remove)){
    if(save_files){ # add Status : AIC exclusion
      Files_Info$status[Files_Info$File %in% studies_to_remove] =
        "Excluded during multivariate MR (p-value > 0.05)"
    }
    tmp = paste0(paste0(studies_to_remove, collapse=" - "), " : removed because of MR p-value > 0.05 \n")
    Log = update_log(Log, tmp, verbose)
    significant.studies=significant.studies[!significant.studies %in% studies_to_remove]
  }

  # Check directionality
  uni.coefs.collection[ significant.studies ] -> corresponding.uvar.coefs
  mycoefs = coefs[coefs$nm %in% significant.studies]
  mycoefs$nm = as.character(mycoefs$nm)
  stopifnot(mycoefs$nm == corresponding.uvar.coefs$nm)
  (mycoefs$Estimate / corresponding.uvar.coefs$Estimate) -> look.for.magnitude.of.diference
  studies_to_removeD = NULL
  while( min(look.for.magnitude.of.diference) < 0) {
    which.min(look.for.magnitude.of.diference) -> study.with.biggest.sign.difference
    studies_to_remove <- c(studies_to_remove, coefs$nm[study.with.biggest.sign.difference])
    studies_to_removeD <- c(studies_to_removeD, coefs$nm[study.with.biggest.sign.difference])
    reasons          <- c(reasons, 'different sign of effect estimate between univariate and multivariate regression')
    look.for.magnitude.of.diference = look.for.magnitude.of.diference[-study.with.biggest.sign.difference]
    if(length(look.for.magnitude.of.diference)==0){
      if(save_files){ # add Status : Direction
        Files_Info$status[Files_Info$File %in% studies_to_removeD] =
          "Excluded during multivariate MR (unconsistent direction)"
      }

      tmp = "No study significant after checking for direction - Analysis failed"
      Log = update_log(Log, tmp, verbose)

      if(save_files)                write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )


      res=list()
      res$log_info = Log
      res$stop = T
      return(res)
      }
  }

  if(!is.null(studies_to_removeD)){
    if(save_files){ # add Status : Direction
      Files_Info$status[Files_Info$File %in% studies_to_removeD] =
        "Excluded during multivariate MR (unconsistent direction)"
    }
    tmp = paste0(paste0(studies_to_removeD, collapse=" - "), " : removed because of different directions (univariate vs multivariate) \n")
    Log = update_log(Log, tmp, verbose)
    significant.studies=significant.studies[!significant.studies %in% studies_to_removeD]
  }


  tmp = paste0("Done! \n")
  Log = update_log(Log, tmp, verbose)

  final_set_of_study_names  = significant.studies[!significant.studies %in% studies_to_remove]


  # compute, and then save, the 'full-22-in-one-go' regression results
  #catn('regression on entire dataset:')
  tmp = paste0("# Final regression... \n")
  Log = update_log(Log, tmp, verbose)

  coefs <- data.frame(coef(summary(lm(data=ZMatrix, formula = generate.formula(All_study_names[length(All_study_names)], final_set_of_study_names)))))
  coefs=cbind(nm=rownames(coefs), coefs)
  rownames(coefs) <- NULL
  coefs = coefs[order(coefs$Pr...t..),,drop=F]



  colnames(coefs) = c("study", "estimate", "std_error", "T", "P")

  if(save_files){
    Files_Info$multi_estimate = numeric()
    Files_Info$multi_estimate[match(coefs$study, Files_Info$File)] = coefs$estimate
    Files_Info$multi_std_error = numeric()
    Files_Info$multi_std_error[match(coefs$study, Files_Info$File)] = coefs$std_error
    Files_Info$multi_T = numeric()
    Files_Info$multi_T[match(coefs$study, Files_Info$File)] = unlist(coefs[,"T"])
    Files_Info$multi_P = numeric()
    Files_Info$multi_P[match(coefs$study, Files_Info$File)] = coefs$P
  }
  if(power){
    Files_Info$R2_adj = summary(lm(data=ZMatrix, formula = generate.formula(All_study_names[length(All_study_names)], final_set_of_study_names)))$adj.r.squared
    Files_Info$OoS_R2_adj_mean = numeric()
    Files_Info$OoS_R2_adj_median = numeric()

    tmp = "Getting out-of-sample R2 masking one chromosome \n"
    Log = update_log(Log, tmp, verbose)

    outcome = colnames(ZMatrix)[ncol(ZMatrix)]

    generate.formula <- function(outcome, study_names ) {
      paste(paste0(outcome,' ~ -1 + '),paste(collapse=' + ', paste(sep='','`',study_names,'`')))
    }

    for(chrm in 1:22) {
      tmp = paste0("   Chromosome ", chrm, "\n")
      Log = update_log(Log, tmp, verbose)
      # create the dataset without this chromosome
      train    =   ZMatrix$chrm != chrm
      d_masked = ZMatrix[train,,drop=F]

      tmp = "Running regression, \n"
      Log = update_log(Log, tmp, verbose)


      lm(data=d_masked, formula = generate.formula(outcome, final_set_of_study_names)
      ) -> fit_masked # fit, without one chromosome
      coefs <- data.frame(coef(summary(fit_masked)))

      # stopifnot(length(dynamic.study.names) ==  nrow(coefs))
      ## this can happen, if all non-0 Z-score of a study are one the same chromosome
      ## but this should not be a problem
      # just add some warning message here :
      if(length(final_set_of_study_names) !=  nrow(coefs)){
        missingSt = final_set_of_study_names[!final_set_of_study_names %in% rownames(coefs)]
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

      # all.coefs = rbind(all.coefs, coefs.DT)

      # check the predictions using data from the training set
      #in_train = data.frame( CRG=d_masked[,..outcome], prd=predict.lm(fit_masked) )
      #residual.variance.in.training=var(in_train[,1]-in_train$prd)


      chrm_ = chrm

      tmp = "Calculating out of sample prediction for SNPs on this chromosome, \n"
      Log = update_log(Log, tmp, verbose)


      # all SNPs we need to predict (the ones on this chromosome)
      d_test             = ZMatrix[chrm==chrm_          ,,drop=T]


      predict(fit_masked, d_test, se.fit=T) -> preds



      test.outcome    <- d_test$lifegen_phase2_bothpl_alldr_2017_09_18.tsv_withZ.gz

      SS.total      <- sum((test.outcome - mean(test.outcome))^2)
      SS.regression <- sum((preds$fit - mean(test.outcome))^2)


      # fraction of variability explained by the model : SS.regression/SS.total
      Files_Info[,paste0("0oS_R2_adj_",chrm) := 1-(1-SS.regression/SS.total)*(nrow(d_masked)-1)/(nrow(d_masked)-length(final_set_of_study_names)-1)]

      tmp = paste0("Done! \n")
      Log = update_log(Log, tmp, verbose)
    }
    Files_Info$OoS_R2_adj_mean = mean(unlist(Files_Info[1,(ncol(Files_Info)-21):ncol(Files_Info)]))
    Files_Info$OoS_R2_adj_median = median(unlist(Files_Info[1,(ncol(Files_Info)-21):ncol(Files_Info)]))
  }




  if(save_files)                write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )

  res=list()
  res$log_info = Log
  res$studies = data.table::data.table(study_selected=final_set_of_study_names)
  res$coeffs = coefs
  if(power){
    res$R2_adj = summary(lm(data=ZMatrix, formula = generate.formula(All_study_names[length(All_study_names)], final_set_of_study_names)))$adj.r.squared
  }


  return(res)
}
