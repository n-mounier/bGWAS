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



identify_studiesMR <- function(ZMatrix, MR_shrinkage, MR_threshold, Z_Matrices, save_files=FALSE, verbose=FALSE){

  Log = c()

  if(save_files) Files_Info = data.table::fread("PriorGWASs.tsv")


  # if p-value studies > MR_shrinkage , Z-Score is set to zero


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
      ZMatrix[abs(ZMatrix[,column_of_zs]) < threshold , column_of_zs] <- 0
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
  tmp= "  Number of trait-specific instruments per univariate regression: \n"
  Log = update_log(Log, tmp, verbose)
  
  Zlimit = qnorm(MR_threshold/2, lower.tail = F)
  
  
  for(one.study in Prior_study_names) {
    # create a ZMatrix_uni containing instruments only for one.study
    ZMatrix_uni = ZMatrix[abs(ZMatrix[,one.study])>Zlimit,]
    tmp= paste0("  . ",one.study, " : ", nrow(ZMatrix_uni), " \n")
    Log = update_log(Log, tmp, verbose)
    
    paste(All_study_names[length(All_study_names)], ' ~ -1 + `',one.study,'`', sep='', collapse="") -> uni.form
    lm(data=ZMatrix_uni, formula = uni.form) -> uni.fit
    uni.coefs <- data.frame(coef(summary(uni.fit)))
    uni.coefs=cbind(nm=rownames(uni.coefs) %>% gsub("`","",.), uni.coefs)
    rownames(uni.coefs) <- NULL
    uni.coefs = uni.coefs[order(uni.coefs$Pr...t..),,drop=F]
    uni.coefs$adj.r.squared <- summary(uni.fit)$adj.r.squared
    uni.coefs$r.squared     <- summary(uni.fit)$r.squared
    uni.coefs.collection = rbind(uni.coefs.collection, uni.coefs)
  }
  data.table::setkey(uni.coefs.collection, nm)
  
  
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
  
  
  
  tmp = paste0("# Stepwise selection (all traits)... \n")
  Log = update_log(Log, tmp, verbose)
  
  
  backward_threshold = 0.05/length(Prior_study_names)
  
  
  uni.coefs.collection=as.data.frame(uni.coefs.collection)
  uni.coefs.collection$study = as.character(uni.coefs.collection$study)
  
  if(all(uni.coefs.collection$P>backward_threshold)){
    tmp = "No study significant - Analysis failed"
    Log = update_log(Log, tmp, verbose)
    
    if(save_files)                write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )
    
    res=list()
    res$log_info = Log
    res$stop = T
    return(res)
  }
  
  
  significant.studies = uni.coefs.collection$study[which.min(uni.coefs.collection$P)]
  
  
  studies.to.test = unlist(subset(uni.coefs.collection, P<0.05, study))
  studies.to.test = studies.to.test[!studies.to.test %in% significant.studies]
  
  non.significant.studies = c(unlist(c(subset(uni.coefs.collection, !P<0.05, study), subset(uni.coefs.collection, is.na(P), study))),studies.to.test)
  
  if(save_files){ # add Status : AIC exclusion
    Files_Info$status[Files_Info$File %in% non.significant.studies] = "Excluded during stepwise selection"
  }
  
  steps = data.frame(Direction=character(), Study=character())
  
  Convergence=F
  
  it = 0
  while(!Convergence){
    if(nrow(steps)>7){ # are we in a loop?
      if(steps$Study[nrow(steps)-6]==steps$Study[nrow(steps)] & steps$Study[nrow(steps)-7]==steps$Study[nrow(steps)-1]){
        tokeep = steps$Study[nrow(steps)]
        toexclude = steps$Study[nrow(steps)-1]
        # add the last ones to our list of study
        non.significant.studies = non.significant.studies[!non.significant.studies %in% tokeep]
        significant.studies = c(significant.studies, tokeep)
        # remove the second last one
        significant.studies = significant.studies[!significant.studies %in% toremove]
        studies.to.test = studies.to.test[!studies.to.test %in% toremove]
        non.significant.studies = c(non.significant.studies, toremove)
        
        tmp=paste0("We are in a loop: ", toremove, " will be excluded from our set of studies and ", tokeep, " will be added back in the model \n")
        update_log(Log, tmp, verbose)
        
      }
      
    }
    
    it=it+1
    tmp=paste0("  iteration ", it, ": ", length(significant.studies), " studies \n")
    update_log(Log, tmp, verbose)
    
    ### UPDATE Z-MATRIX
    my.formula = generate.formula(All_study_names[length(All_study_names)], significant.studies)
    
    # keep only our studies
    ZMatrix_subset = ZMatrix
    #print(non.significant.studies)
    ZMatrix_subset[,non.significant.studies] <- NULL
    
    # remove instruments non specific to our studies
    Zlimit = qnorm(MR_threshold/2, lower.tail = F)
    if(length(significant.studies)>1){
      SNPsToKeep = apply(ZMatrix_subset[,-c(1:5,as.numeric(ncol(ZMatrix_subset)))], 1, function(x) any(abs(x)>Zlimit))
      ZMatrix_subset=ZMatrix_subset[SNPsToKeep,]
    } else {
      SNPsToKeep = ZMatrix_subset[,-c(1:5,as.numeric(ncol(ZMatrix_subset)))]>Zlimit
      ZMatrix_subset=ZMatrix_subset[SNPsToKeep,]
    }
    
    ## RUN MODEL
    tmp="#Run model \n"
    update_log(Log, tmp, verbose)
    
    model = lm(data=ZMatrix_subset, formula = my.formula) 
    lm(data=ZMatrix_subset, formula = my.formula) %>%summary %>%coef %>% data.table::data.table(keep.rownames=T) -> coefs
    data.table::setnames(coefs,"rn","nm")
    coefs[, nm := gsub("`","",nm)]
    
    no_change = T
    
    tmp = paste0("#Test if any study can be added with p<", round(backward_threshold, 4), " \n")
    Log = update_log(Log, tmp, verbose)
    PValues = numeric(length(studies.to.test))
    PValues = data.frame(study=studies.to.test, P=numeric(length(studies.to.test)), Estimate=numeric(length(studies.to.test)), stringsAsFactors = F)
    for(one.to.add in studies.to.test){
      ### UPDATE Z-MATRIX
      test_formula = generate.formula(All_study_names[length(All_study_names)], c(significant.studies, one.to.add))
      # keep only our studies
      ZMatrix_test = ZMatrix
      studies_not_used = non.significant.studies[!non.significant.studies %in% one.to.add]
      ZMatrix_test[,studies_not_used] <- NULL
      
      # remove instruments non specific to our studies
      Zlimit = qnorm(MR_threshold/2, lower.tail = F)
      # DO NOT USE THE LAST COLUMN!!
      SNPsToKeep = apply(ZMatrix_test[,-c(1:5,as.numeric(ncol(ZMatrix_test)))], 1, function(x) any(abs(x)>Zlimit))
      ZMatrix_test=ZMatrix_test[SNPsToKeep,]
      
      
      ## RUN MODEL
      model = lm(data=ZMatrix_test, formula = test_formula) 
      
      PValues[one.to.add] = coef(summary(model))[one.to.add,"Pr(>|t|)"]
      PValues[PValues$study==one.to.add,c(2,3)] = c(coef(summary(model))[one.to.add,"Pr(>|t|)"], coef(summary(model))[one.to.add,"Estimate"])
    }
    
    
    
    
    
    if(any(PValues$P<backward_threshold)){
      no_change = F
      study_to_add = PValues$study[which.min(PValues$P)]
      
      
      sign = uni.coefs.collection$estimate[uni.coefs.collection$study==study_to_add]/PValues$Estimate[which.min(PValues$P)]
      
      if(sign<0){
        tmp = paste0("Study :", study_to_add, " cannot be added, direction is not consistent between univariate and multivariate model.\n")
        Log = update_log(Log, tmp, verbose)
        
        studies.to.test = studies.to.test[-which(studies.to.test %in% study_to_add)]
        if(save_files){ # add Status : AIC exclusion
          Files_Info$status[Files_Info$File %in% study_to_add] = "Excluded during forward selection (unconsistent direction)"
        }
        
        steps=rbind(steps, data.frame(Direction="+", Study=""))
        
        
        
      } else {
        
        tmp = paste0("Adding one study :", study_to_add, " \n")
        Log = update_log(Log, tmp, verbose)
        
        steps=rbind(steps, data.frame(Direction="+", Study=study_to_add))
        
        
        significant.studies=c(significant.studies, study_to_add)
        studies.to.test = studies.to.test[-which(studies.to.test %in% study_to_add)]
        non.significant.studies = non.significant.studies[-which(non.significant.studies %in% study_to_add)]
        if(save_files){ # add Status : USE
          Files_Info$status[Files_Info$File %in% study_to_add] = "USED"
        }
        
        
        tmp = paste0("Done! \n")
        Log = update_log(Log, tmp, verbose)
      }
      
    } else {
      steps=rbind(steps, data.frame(Direction="+", Study=""))
      
    }
    
    
    # Add
    
    my.formula = generate.formula(All_study_names[length(All_study_names)], significant.studies)
    
    # keep only our studies
    ZMatrix_subset = ZMatrix
    #print(non.significant.studies)
    ZMatrix_subset[,non.significant.studies] <- NULL
    
    # remove instruments non specific to our studies
    Zlimit = qnorm(MR_threshold/2, lower.tail = F)
    # DO NOT USE THE LAST COLUMN!!
    SNPsToKeep = apply(ZMatrix_subset[,-c(1:5,as.numeric(ncol(ZMatrix_subset)))], 1, function(x) any(abs(x)>Zlimit))
    ZMatrix_subset=ZMatrix_subset[SNPsToKeep,]
    
    ## RUN MODEL
    tmp="#Update model \n"
    update_log(Log, tmp, verbose)
    
    model = lm(data=ZMatrix_subset, formula = my.formula) 
    lm(data=ZMatrix_subset, formula = my.formula) %>%summary %>%coef %>% data.table::data.table(keep.rownames=T) -> coefs
    data.table::setnames(coefs,"rn","nm")
    coefs[, nm := gsub("`","",nm)]
    
    
    tmp = paste0("#Test if any study has p>", round(backward_threshold, 4), " now \n")
    Log = update_log(Log, tmp, verbose)
    
    
    # Remove traits with multivariate p-value larger that 0.05
    coefs[order(-`Pr(>|t|)`)] -> largest.MR.p.value
    study_to_remove <- character()
    if(any(largest.MR.p.value$`Pr(>|t|)` > backward_threshold)){
      no_change = F
      
      study_to_remove <-  largest.MR.p.value$nm[1]
      
      steps=rbind(steps, data.frame(Direction="-", Study=study_to_remove))
      
      tmp = paste0("Excluding one study :", study_to_remove, " \n")
      Log = update_log(Log, tmp, verbose)
      
      significant.studies=significant.studies[!significant.studies %in% study_to_remove]
      non.significant.studies = c(non.significant.studies,study_to_remove)
      
      if(save_files){ # add Status : AIC exclusion
        Files_Info$status[Files_Info$File %in% study_to_remove] = "Excluded during stepwise selection"
      }
      
      
      tmp = paste0("Done! \n")
      Log = update_log(Log, tmp, verbose)
      
      
    } else {
      steps=rbind(steps, data.frame(Direction="-", Study=""))
      
    } 
    
    
    if(no_change){
      tmp = paste0("It converged! \n")
      Log = update_log(Log, tmp, verbose)
      Convergence=T
    }
    
    
    
  }
  
  
  

  final_set_of_study_names  = significant.studies
  

  # compute, and then save, the 'full-22-in-one-go' regression results
  #catn('regression on entire dataset:')
  tmp = paste0("# Final regression... \n")
  Log = update_log(Log, tmp, verbose)

  coefs <- data.frame(coef(summary(lm(data=ZMatrix_subset, formula = generate.formula(All_study_names[length(All_study_names)], final_set_of_study_names)))))
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
  outcome = colnames(ZMatrix_subset)[ncol(ZMatrix_subset)]
  
  tmp = paste0("Estimating adjusted R-squared: \n")
  Log = update_log(Log, tmp, verbose)
  R2_Multi = summary(lm(data=ZMatrix_subset, formula = generate.formula(outcome, final_set_of_study_names)))$adj.r.squared
  tmp = paste0("- in-sample adjusted R-squared for the all-chromosomes multivariate regression is ", round(R2_Multi,4), " \n")
  Log = update_log(Log, tmp, verbose)
  tmp = paste0("- out-of-sample R-squared (masking one chromosome at a time), for the multivariate regression will be estimated when calculating the prior. \n")
  Log = update_log(Log, tmp, verbose)


  if(save_files)                write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )

  res=list()
  res$log_info = Log
  res$studies = final_set_of_study_names
  res$coeffs = coefs
  res$ZMat = ZMatrix_subset
  

  return(res)
}
