###### Function to compute prior for each SNP######



# #' Compute Prior
# #' 
# NOT EXPORTED

compute_prior <- function(Data, parent0, verbose){ # Data is output of tidy_inputGWAS
  Log = c()
  #rs,chrm,pos,alt,ref,observed_Z,prior_estimate,prior_std_error,posterior_estimate,posterior_std_error
  if(parent0 == "mat0"){
    my_prior = -2*Data$GWAS_add$std_beta
  } else if(parent0 == "pat0"){
    my_prior = 2*Data$GWAS_add$std_beta
  } else {
    stop("wrong model")
  }
  
  
  Data$GWAS_diff %>%
    mutate(mu_prior_estimate = my_prior,
           mu_prior_std_error = Data$GWAS_add$std_se) -> Results
  
  
  
  if(parent0=="mat0"){
    data.table::fwrite(Results, "Prior_mat0.csv.gz", sep = ",", compress="gzip")
        tmp = paste0("The file Prior_mat0.csv.gz had been successfully written. \n")
  } else {
    data.table::fwrite(Results, "Prior_pat0.csv.gz", sep = ",", compress="gzip")
    tmp = paste0("The file Prior_pat0.csv.gz had been successfully written. \n")
  }
  Log = bGWAS:::update_log(Log, tmp, verbose)
  
  
  res=list(log_info = Log,
           prior = Results) 
  return(res)
}
