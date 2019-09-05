###### Function to identifiy studies significantly affecting trait of interest - multivariate MR ######



# #' Identify studies for prior using multivariate MR
# #'
# #' From a pruned Z-matrix of strong MR instruments, performs multivariate MR
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED




identify_studiesMR <- function(ZMatrix, MR_shrinkage, MR_threshold, stepwise_threshold, Z_matrices, save_files=FALSE, verbose=FALSE){
  
  Log = c()
  
  if(save_files) Files_Info = readr::read_tsv("PriorGWASs.tsv", progress = FALSE, col_types = readr::cols())
  
  
  tmp = paste0("#Preparation of the MR analyses to identify significant studies... \n")
  Log = update_log(Log, tmp, verbose)
  
  names(ZMatrix) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>%
    filter(!.data$. %in% c('rs','chrm','pos','alt','ref')) %>%
    pull() ->  All_study_names

  All_study_names %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    slice(-length(.data$.)) %>%
    pull() ->  Prior_study_names 

   
  tmp = paste0("Conventionnal GWAS of interest : ", All_study_names[length(All_study_names)], "\n")
  Log = update_log(Log, tmp, verbose)
  
  
  # Function to automatically generate the formula for linear model
  generate_Formula <- function(outcome, study_names, with_intercept=F ) {
    formula = ifelse(with_intercept,
                     paste(outcome, ' ~  1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))),
                     paste(outcome, ' ~ -1 + ',paste(collapse=' + ', paste(sep='','`',study_names,'`'))))
    return(formula)
  }
  
  # Function to automatically get subsetted matrix of instruments
  get_Instruments <- function(Zmat, Zlim){
    # get exposures
    Zmat %>%
      names() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      slice(6:(ncol(Zmat)-1)) %>%
      pull() -> exposures

    # subset : any row, with abs(exposure)>Zlimit
    Zmat %>%
      filter_at(vars(exposures), any_vars(abs(.data$.)>Zlim)) -> Zmat
    return(Zmat)
  }
  
  # Compute and save the univariate regressions (used to check directionnality in multivariate regression):
  tmp = paste0("# Univariate regressions for each trait... \n")
  Log = update_log(Log, tmp, verbose)
  
  tmp= "  Number of trait-specific instruments per univariate regression: \n"
  Log = update_log(Log, tmp, verbose)
  
  Zlimit = stats::qnorm(MR_threshold/2, lower.tail = F)
  
  all_uni_coefs = data.frame()
 
   All_study_names %>%
     as.data.frame(stringsAsFactors=F) %>%
     slice(length(All_study_names)) %>%
     pull() -> my_outcome
 
  for(one_study in Prior_study_names) {
    # create a ZMatrix_uni containing instruments only for one_study
    ZMatrix %>%
      filter(abs(eval(parse(text=one_study))) > Zlimit) -> ZMatrix_uni
    
    tmp= paste0("  . ", get_names(one_study, Z_matrices), " : ", nrow(ZMatrix_uni), " \n")
    Log = update_log(Log, tmp, verbose)
    
    paste(my_outcome, ' ~ -1 + `',one_study,'`', sep='', collapse="") -> uni_form
    stats::lm(data=ZMatrix_uni, formula = uni_form) -> uni_fit
    
    uni_coefs <- data.frame(stats::coef(summary(uni_fit)))
    uni_coefs %>%
      tibble::rownames_to_column("nm") -> uni_coefs
    
    uni_coefs %>%
      mutate(adj.r.squared = summary(uni_fit)$adj.r.squared,
             r.squared = summary(uni_fit)$r.squared) -> uni_coefs
    all_uni_coefs %>%
      bind_rows(uni_coefs) -> all_uni_coefs
  }
  
  all_uni_coefs %>%
    set_names(c("study", "estimate", "std_error", "Tstat", "P",
                "adj_Rsquared", "Rsquared")) -> all_uni_coefs
  
  if(save_files){ # add univariate coeffs
    order_inFiles = match(all_uni_coefs$study, Files_Info$File)
    Files_Info$uni_estimate = NA
    Files_Info$uni_estimate[order_inFiles] = pull(all_uni_coefs, .data$estimate)
    Files_Info$uni_std_error = NA
    Files_Info$uni_std_error[order_inFiles] = pull(all_uni_coefs, .data$std_error)
    Files_Info$uni_T = NA
    Files_Info$uni_T[order_inFiles] = pull(all_uni_coefs, .data$Tstat)
    Files_Info$uni_P = NA
    Files_Info$uni_P[order_inFiles] = pull(all_uni_coefs, .data$P)
    Files_Info$uni_adj_Rsquared = NA
    Files_Info$uni_adj_Rsquared[order_inFiles] = pull(all_uni_coefs, .data$adj_Rsquared)
    Files_Info$uni_Rsquared = NA
    Files_Info$uni_Rsquared[order_inFiles] = pull(all_uni_coefs, .data$Rsquared)
  }
  
  tmp = paste0("Done! \n")
  Log = update_log(Log, tmp, verbose)
  
  tmp = paste0("# Stepwise selection (all traits)... \n")
  Log = update_log(Log, tmp, verbose)
  
  if(is.null(stepwise_threshold)){
    stepwise_threshold = 0.05/length(Prior_study_names)
    tmp = paste0("The p-value threshold used for stepwise selection is ", format(round(stepwise_threshold, 4), scientific = F), " (", length(Prior_study_names), " Prior GWASs tested).  \n")
    Log = update_log(Log, tmp, verbose)
  } 

  if(all(all_uni_coefs$P>stepwise_threshold)){
    tmp = "No study significant - Analysis failed"
    Log = update_log(Log, tmp, verbose)
    
    if(save_files) utils::write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )
    
    res=list(log_info = Log,
             stop = T)
    return(res)
  }
  
  all_uni_coefs %>%
    filter(.data$P==min(.data$P)) %>%
    pull(.data$study) -> significant_studies
  
  all_uni_coefs %>%
    filter(!.data$study %in% significant_studies) %>%
    pull(.data$study) -> non_significant_studies
  
  #studies_to_test = non_significant_studies # same now, but some studies (inconsistent direction / loop : might be removed from the test set, and not from the non.significant set)
  
  all_uni_coefs %>%
    filter(!.data$study %in% significant_studies) %>%
    filter(P<0.05) %>%
    pull(.data$study) ->  studies_to_test
  
  tmp = paste0("Studies tested (reaching p<0.05 in univariate models) : \n ", paste(get_names(studies_to_test, Z_matrices), collapse = " \n "), "\n")
  Log = update_log(Log, tmp, verbose)
  
  
  if(save_files){ # add Status : stepwise exclusion, will be updated if the studies are added
    Files_Info$status[Files_Info$File %in% non_significant_studies] = "Excluded during stepwise selection"
  }
  
  steps = data.frame(Direction=character(), Study=character())
  
  Convergence=F
  
  tmp = paste0("Adding the first study :", get_names(significant_studies, Z_matrices) , " \n")
  Log = update_log(Log, tmp, verbose)
  
  it = 0
  while(!Convergence){
    if(it>50){ # are we in a loop?
      tmp=c("We are probably in a loop...\n ", "Analysis failed \n",
            "Please have a look at the stepwise procedure results ", 
            "and update \"prior_studies\" before relaunching the analysis ",
            "to avoid looping again.")
      Log = update_log(Log, tmp, verbose)
      
      if(save_files) utils::write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )
      
      res=list(log_info = Log,
               stop = T)
      return(res)
    }
    
    
    it=it+1
    tmp=paste0("  iteration ", it, ": ", length(significant_studies), " studies \n")
    update_log(Log, tmp, verbose)
    
    ### UPDATE Z-MATRIX
    myFormula = generate_Formula(my_outcome, 
                                 significant_studies)
    
    # keep only our studies
    ZMatrix %>%
      select(-non_significant_studies) %>%
      get_Instruments(Zlim=Zlimit) -> ZMatrix_subset
    
    ## RUN MODEL - with studies already included
    tmp="#Run model \n"
    update_log(Log, tmp, verbose)
    
    model = stats::lm(data=ZMatrix_subset, formula = myFormula) 
    stats::lm(data=ZMatrix_subset, formula = myFormula) %>%
      summary %>%
      stats::coef() %>%
      as.data.frame() %>% # 1st create a data.frame (to get rownames)
      tibble::rownames_to_column("nm") %>% # then convert to tibble
      as_tibble()  -> coefs

    
    
    ## TRY TO ADD
    no_change = T
    
    tmp = paste0("#Test if any study can be added with p<", format(round(stepwise_threshold, 4), scientific = F), " \n")
    Log = update_log(Log, tmp, verbose)
    PValues = data.frame(study=studies_to_test, P=NA, Estimate=NA, stringsAsFactors = F)
    for(one_to_add in studies_to_test){
      ### UPDATE Z-MATRIX
      test_formula = generate_Formula(my_outcome, 
                                      c(significant_studies, one_to_add))
      # keep only our studies
      non_significant_studies[!non_significant_studies %in% one_to_add] -> studies_not_used 
      ZMatrix %>%
        select(-studies_not_used) %>%
        get_Instruments(Zlim=Zlimit) -> ZMatrix_test
      
      ## RUN MODEL
      model_test = stats::lm(data=ZMatrix_test, formula = test_formula) 
      
      PValues[PValues$study==one_to_add,c(2,3)] = c(stats::coef(summary(model_test))[one_to_add,"Pr(>|t|)"], stats::coef(summary(model_test))[one_to_add,"Estimate"])
    }
    
    if(any(PValues$P<stepwise_threshold)){
      no_change = F
      PValues %>%
        filter(.data$P==min(.data$P)) %>%
        pull(.data$study) -> study_to_add 
      
      # test if direction is consistent with uni before adding
      all_uni_coefs %>%
        filter(.data$study==study_to_add) %>%
        pull(.data$estimate) /
        PValues %>% 
        filter(.data$study==study_to_add) %>%
        pull(.data$Estimate) -> ratio
      
      if(ratio<0){
        tmp = paste0("Study :", get_names(study_to_add, Z_matrices), " cannot be added, direction is not consistent between univariate and multivariate model.\n")
        Log = update_log(Log, tmp, verbose)
        
        studies_to_test[!studies_to_test == study_to_add] -> studies_to_test
        
        if(save_files){ # add Status : AIC exclusion
          Files_Info$status[Files_Info$File %in% study_to_add] = "Excluded during forward selection (unconsistent direction)"
        }
        
        steps=rbind(steps, data.frame(Direction="+", Study="none: direction issue"))
      } else {
        
        tmp = paste0("Adding one study :", get_names(study_to_add, Z_matrices), " \n")
        Log = update_log(Log, tmp, verbose)
        
        steps=rbind(steps, data.frame(Direction="+", Study=study_to_add))
        
        
        significant_studies=c(significant_studies, study_to_add)
        studies_to_test = studies_to_test[-which(studies_to_test %in% study_to_add)]
        non_significant_studies = non_significant_studies[-which(non_significant_studies %in% study_to_add)]
        if(save_files){ # add Status : USE
          Files_Info$status[Files_Info$File %in% study_to_add] = "USED"
        }
        
        
        tmp = paste0("Done! \n")
        Log = update_log(Log, tmp, verbose)
      }
      
    } else {
      steps=rbind(steps, data.frame(Direction="+", Study="none"))
      
    }
    
    
    # Add it if needed
    
    myFormula = generate_Formula(my_outcome, 
                                 significant_studies)
    
    # keep only our studies
    ZMatrix %>% 
      select(-non_significant_studies) %>%
      get_Instruments(Zlim=Zlimit) -> ZMatrix_subset
    
    
    ## RUN MODEL
    tmp="#Update model \n"
    update_log(Log, tmp, verbose)
    
    model = stats::lm(data=ZMatrix_subset, formula = myFormula) 
    stats::lm(data=ZMatrix_subset, formula = myFormula) %>%
      summary %>%
      stats::coef() %>%
      as.data.frame() %>% # 1st create a data.frame (to get rownames)
      tibble::rownames_to_column("nm") %>% # then convert to tibble
      as_tibble()  -> coefs
    
    tmp = paste0("#Test if any study has p>", format(round(stepwise_threshold, 4), scientific = F), " now \n")
    Log = update_log(Log, tmp, verbose)
    
    
    # Remove traits with multivariate p-value larger that 0.05
    coefs %>%
      arrange(desc(.data$`Pr(>|t|)`)) -> coefs
    
    if(any(coefs$`Pr(>|t|)` > stepwise_threshold)){
      no_change = F
      
      coefs %>%
        pull(.data$nm) -> st
      st[1] -> study_to_remove 
      
      steps=rbind(steps, data.frame(Direction="-", Study=study_to_remove))
      
      tmp = paste0("Excluding one study :", get_names(study_to_remove, Z_matrices), " \n")
      Log = update_log(Log, tmp, verbose)
      
      significant_studies[!significant_studies %in% study_to_remove] -> significant_studies
      non_significant_studies = c(non_significant_studies,study_to_remove)
      studies_to_test = c(studies_to_test, study_to_remove)
      
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
  
  
  
  
  # compute, and then save, the 'full-22-in-one-go' regression results
  tmp = paste0("# Final regression... \n")
  Log = update_log(Log, tmp, verbose)
  
  tmp = paste0("The studies used are: \n")
  Log = update_log(Log, tmp, verbose)
  
  tmp = c(paste0(paste0("- ", get_names(significant_studies, Z_matrices)), collapse="\n"), "\n")
  Log = update_log(Log, tmp, verbose)
  
  ZMatrix %>% 
    select(1:5, significant_studies, ncol(ZMatrix)) %>%
    get_Instruments(Zlim=Zlimit) -> ZMatrix_final
  myFormula_final = generate_Formula(my_outcome, significant_studies)
  
  stats::lm(data=ZMatrix_final, formula = myFormula_final)  %>%
    summary %>%
    stats::coef() %>%
    as.data.frame() %>% # 1st create a data.frame (to get rownames)
    tibble::rownames_to_column("nm") %>% # then convert to tibble
    as_tibble()  -> coefs
  
  coefs %>%
    set_names(c("study", "estimate", "std_error", "Tstat", "P")) -> coefs
  
  if(save_files){
    order_inFiles = match(coefs$study, Files_Info$File)
    Files_Info$multi_estimate = NA
    Files_Info$multi_estimate[order_inFiles] = pull(coefs, .data$estimate)
    Files_Info$multi_std_error = NA
    Files_Info$multi_std_error[order_inFiles] = pull(coefs, .data$std_error)
    Files_Info$multi_T = NA
    Files_Info$multi_T[order_inFiles] = pull(coefs, .data$Tstat)
    Files_Info$multi_P = NA
    Files_Info$multi_P[order_inFiles] = pull(coefs, .data$P)
  }
  
  tmp = paste0("Estimating adjusted R-squared: \n")
  Log = update_log(Log, tmp, verbose)
  
  stats::lm(data=ZMatrix_final, formula = myFormula_final)  %>%
    summary -> summary_lm
  summary_lm["adj.r.squared"] %>% 
    as.numeric() -> R2_Multi 
  
  tmp = paste0("- in-sample adjusted R-squared for the all-chromosomes multivariate regression is ", round(R2_Multi,4), " \n")
  Log = update_log(Log, tmp, verbose)
  tmp = paste0("- out-of-sample R-squared (masking one chromosome at a time), for the multivariate regression will be estimated when calculating the prior. \n")
  Log = update_log(Log, tmp, verbose)
  
  
  if(save_files) utils::write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )
  
  res=list(log_info = Log,
           studies = significant_studies,
           coeffs = coefs,
           ZMat = ZMatrix_subset)
  return(res)
}


