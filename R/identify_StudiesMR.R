###### Function to identifiy studies significantly affecting trait of interest - multivariate MR ######



#' Identify studies for prior using multivariate MR
#'
#' From a pruned Z-matrix of strong MR instruments, performs multivariate MR
#'
#' @inheritParams bGWAS
#' @param ZMatrix The pruned Z-Matrix of strong instrument or the path of the file
#'        containing it (data table or character)
#'
#' @return Log file + list of studies significant + create files (list them) if saveFiles=T
#'
#' @importFrom magrittr "%>%"
#'
#' @export


identify_StudiesMR <- function(ZMatrix, saveFiles=FALSE, verbose=FALSE){

  Log = c()

  scheme_INpXXtozero = 1e-5 # if p-value studies > 1e-5, Z-Score is set to zero

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
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  # Truncate Z-scores
  ZMatrix= push.extreme.zs.back.a.little.towards.zero(ZMatrix)



  All_study_names    = names(ZMatrix) %..% ..[!.. %in% c('rs','chrm','pos','alt','ref')]

  Prior_study_names = All_study_names[1:length(All_study_names)-1]
  tmp = paste0("Studies tested : ", paste(Prior_study_names, collapse = " - "), "\n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  tmp = paste0("Conventionnal GWAS of interest : ", All_study_names[length(All_study_names)], "\n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  # Set the z-scores to 0 for the regression if too low
  if(scheme_INpXXtozero < 1.0) {
    for(column_of_zs in Prior_study_names) {
      threshold = abs(qnorm(scheme_INpXXtozero/2))
      ZMatrix[c(abs(ZMatrix[,..column_of_zs]) < threshold) , column_of_zs] <- 0
    }
    #tmp = "Z-scores were set to 0 if p-value > 1e-5 \n"
    #Log = c(Log, tmp)
    #if(verbose) cat(tmp)
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
  tmp = paste0("Univariate regressions for each trait... \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


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
    S = availableStudies()
    IDs = S$ID[match(Names$nm, S$File)]
    Ref = paste0(Names$nm, " (ID = ", IDs, " - r2 = ", r2$r.squared, ")")

    stop(paste0("The study ", Ref,
                " has a really high r squared, it probably corresponds to the",
                " same trait than your conventionnal GWAS - please remove it",
                " before re-rerunning the analysis. \n"))
  }


  if(saveFiles){
    readr::write_csv("Univariate_coefs.csv", x=uni.coefs.collection)
    tmp = "The file \"Univariate_coefs.csv\" has been successfully writed. \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }


  tmp = paste0("Done! \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  tmp = paste0("# Stepwise AIC multivariate regression...")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  significant.studies = Prior_study_names # might be shrunk in the coming lines
  stopifnot(length(significant.studies) == length(Prior_study_names))

  initial.formula = paste0(All_study_names[length(All_study_names)],' ~ -1')
  k=2 # for AIC
  # step looks in global environment for the data it needs, trick : using the assign function
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
  #cat(coefs)
  significant.studies = gsub("`","",rownames(coefs))



  post_AIC_studies = significant.studies
  if(length(post_AIC_studies)==0) stop("no study significant")

  if(length(post_AIC_studies)>1){
    tmp = paste0(length(post_AIC_studies), " studies are significant : ", paste0(post_AIC_studies, collapse = " - "), "\n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  } else {
    tmp = paste0(length(post_AIC_studies), " study is significant : ", paste0(post_AIC_studies, collapse = " - "), "\n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  "%|%" <- function(x,y){
    do.call(y,list(substitute(x)),envir=parent.frame()) # just right
  }

  tmp = paste0("Done! \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  # We now have a set of studies from AIC,but we may want to prune them further
  tmp = paste0("#Further pruning of these studies... \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  removal_reasons=data.table::data.table()

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
    if(nrow(largest.MR.p.value)==0) stop("no study significant in multivariate analysis with p<0.05")
  }

  if(!is.null(studies_to_remove)){
    tmp = paste0(paste0(studies_to_remove, collapse=" - "), " : removed because of MR p-value > 0.05 \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
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
    if(length(look.for.magnitude.of.diference)==0) stop("no study significant after checking for direction")
  }

  if(!is.null(studies_to_removeD)){
    tmp = paste0(paste0(studies_to_removeD, collapse=" - "), " : removed because of different directions (univariate vs multivariate) \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }


  if(!is.null(studies_to_remove)){
    removal_reasons %<>% rbind(.,data.table::data.table(studies_to_remove, reasons))

    if(saveFiles){
      readr::write_csv(removal_reasons, "Removal_reasons.csv")
      tmp = paste0("The file ", "Removal_reasons.csv has been successfully writed. \n")
      Log = c(Log, tmp)
      if(verbose) cat(tmp)
    }
  }

  tmp = paste0("Done! \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  final_set_of_study_names  = significant.studies[!significant.studies %in% studies_to_remove]


  # compute, and then save, the 'full-22-in-one-go' regression results
  #catn('regression on entire dataset:')
  tmp = paste0("# Final regression... \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  coefs <- data.frame(coef(summary(lm(data=ZMatrix, formula = generate.formula(All_study_names[length(All_study_names)], final_set_of_study_names)))))
  coefs=cbind(nm=rownames(coefs), coefs)
  rownames(coefs) <- NULL
  coefs = coefs[order(coefs$Pr...t..),,drop=F]

  if(saveFiles){
    readr::write_csv(path="Multivariate_coefs.csv", x=coefs)
    tmp = paste0("The file ", "Multivariate_coefs.csv has been successfully writed. \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  tmp = paste0("Done! \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  ## Do we really need it ?
#  coefs <- data.frame(coef(summary(lm(data=ZMatrix, formula = generate.formula(All_study_names[length(All_study_names)], final_set_of_study_names, with.intercept=T)))))
#  coefs=cbind(nm=rownames(coefs), coefs)
#  rownames(coefs) <- NULL
#  coefs = coefs[order(coefs$Pr...t..),,drop=F]
  #cat(coefs)
  #readr::write_csv(path="univariate.coefs.withInter.csv", x=coefs)
#  if(saveFiles){
#    readr::write_csv(data.table::data.table(study_selected=final_set_of_study_names), "MR_Studies.csv")
#    tmp = paste0("The file ", "MR_Studies.csv has been successfully writed. \n")
#    Log = c(Log, tmp)
#    if(verbose) cat(tmp)
#  }

  res=list()
  res$Log = Log
  res$Studies = data.table::data.table(study_selected=final_set_of_study_names)
  return(res)
}
