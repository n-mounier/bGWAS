###### Function to create a nice data.frame for any type of input GWAS######



# #' Tidy input GWAS
# #'
# #' From the GWAS argument of the main bGWAS function, create a nice/tidy data.frame
# #' that can be used by all other functions.
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED



tidy_inputGWAS <- function(GWAS_PofO, GWAS_additive, verbose){
  Log = c()
  
  
  tmp = paste0("# Preparation of the data... \n")
  Log = bGWAS:::update_log(Log, tmp, verbose)
  
 
  # from PofO get "diff"
  my_colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", "se", "p", "N")
  
  GWAS_PofO %>%
    select(c(1:3,5,6), dplyr::ends_with("COMP"), 21) %>%
    stats::setNames(my_colNames) %>%
    mutate(p=NULL,
           z = beta/se,
           std_beta = z/sqrt(N),
           std_se = 1/sqrt(N)) -> GWAS_diff
  
  # align GWAS_additive + subset to SNPs in common
  GWAS_add = data.table::fread(GWAS_additive, showProgress = FALSE, data.table=F) %>%
    transmute(rsid=SNP, allele=ALLELE, beta=BETA, se=SE, N, z=beta/se, std_beta = z/sqrt(N),
              std_se = 1/sqrt(N))
  data <- dplyr::inner_join(GWAS_diff, GWAS_add, by=c("rsid"), suffix=c(".PofO", ".additive"))
  
  data %>%
    dplyr::mutate(std_beta.additive = dplyr::case_when(
      allele == alt ~ std_beta.additive, # aligned (we only have one allele to compare to)
      allele == ref ~ -std_beta.additive, # swapped (we only have one allele to compare to)
      TRUE ~ NA_real_ # otherwise exclude SNP
    ) ) %>%
    dplyr::filter(!is.na(std_beta.additive)) -> data
  # nrow(data) 5,349,118
  
  
  # final data
  myfinal_colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", "se", "N", "z", "std_beta", "std_se")
  
  data %>% 
    select(1:5, ends_with(".PofO")) %>%
    setNames(myfinal_colNames) -> GWAS_diff_clean
  
  data %>% 
    select(1:5, ends_with(".additive")) %>%
    setNames(myfinal_colNames) -> GWAS_add_clean
  
  res=list(log_info = Log,
           GWAS_add = GWAS_add_clean,
           GWAS_diff = GWAS_diff_clean)
  return(res)
}
