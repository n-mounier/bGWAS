###### Function to create full Z-matrix to compute prior ######



# #' Create full Z-matrix for prior computation
# #'
# #' From a list of significant studies, create the full Z-matrix that can be used to compute
# #' the prior.
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED





makeFull_ZMatrix <- function(studies=NULL, GWASData, GName,  Z_matrices="~/Z_matrices", 
                             prior_shrinkage, save_files=F, verbose=F) {

  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = update_log(Log, tmp, verbose)
  
  
  tmp = paste0("Selecting studies :\n")
  Log = update_log(Log, tmp, verbose)
  
  ZMatrix = as_tibble(data.table::fread(file.path(Z_matrices, "ZMatrix_Full.csv.gz"), select=c(1:5, studies+5 ), showProgress = FALSE, data.table = F))

  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = update_log(Log, tmp, verbose)
  
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)
  
  
  
  
  # Add conventional GWAS column, at the end (make sure alleles are aligned)
  tmp = paste0("# Adding data from the conventional GWAS : \n \"", GName,
               "\" \n")
  Log = update_log(Log, tmp, verbose)
  # keep the SNPs in our Z matrix and order them correctly + check allele alignement
  ZMatrix %>%
    filter(.data$rs %in% GWASData$rsid) -> ZMatrix
  
  GWASData %>%
    slice(match(ZMatrix$rs, .data$rsid)) %>%
    mutate( ZMat_alt = ZMatrix$alt,
            Zmat_ref = ZMatrix$ref,
            aligned_Z = case_when(
              (.data$alt == .data$ZMat_alt &
                 .data$ref == .data$Zmat_ref) ~ .data$z_obs,
              (.data$ref == .data$ZMat_alt &
                 .data$alt == .data$Zmat_ref) ~ -.data$z_obs,
              TRUE ~ NA_real_))-> GWASData
  
  ZMatrix %>%
    mutate({{GName}} := GWASData$aligned_Z) -> ZMatrix
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  # set NA to zero for prediction
  ZMatrix %>%
    mutate_all(~replace(., is.na(.), 0)) -> ZMatrix
  
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS \n")
  Log = update_log(Log, tmp, verbose)
  
  
  
  push_extreme_zs_back_a_little_towards_zero <- function(d) { # Some z-scores are just too far from zero
    maxAllowed_z = abs(stats::qnorm(1e-300 / 2)) # p=1e-300 is the max allowed now, truncate z-scores accordingly
    names(d) %>% 
      .[!. %in% c("rs","chrm","pos","alt","ref")] -> studies_here
    for(n in studies_here) {
      d %>%
        mutate(!!n := case_when(
          abs(eval(parse(text=n))) > maxAllowed_z  ~  maxAllowed_z,
          abs(eval(parse(text=n))) < -maxAllowed_z ~ -maxAllowed_z,
          TRUE ~ eval(parse(text=n)))) -> d
    }
    return(d)
  }
  
  # Truncate Z-scores
  ZMatrix %>%
    push_extreme_zs_back_a_little_towards_zero() -> ZMatrix
  
  
  
  
  # Set the z-scores to 0 for the regression if shrinkage
  if(prior_shrinkage < 1.0) {  
    names(ZMatrix) %>% 
      .[!. %in% c('rs','chrm','pos','alt','ref', GName)] -> Prior_study_names
    threshold = abs(stats::qnorm(prior_shrinkage/2))
    for(column_of_zs in Prior_study_names) { 
      ZMatrix %>%
        mutate(!!column_of_zs := case_when(
          abs(eval(parse(text=column_of_zs))) < threshold ~ 0,
          TRUE ~ eval(parse(text=column_of_zs)))) -> ZMatrix
    }
    tmp = paste0("Applying shrinkage (threshold = ", prior_shrinkage, ") before calculating the prior. \n")
    Log = update_log(Log, tmp, verbose)
  }
  
  
  res=list(
    log_info = Log,
    mat = ZMatrix)
  return(res)
}

