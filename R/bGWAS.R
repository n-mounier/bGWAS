###### Main Function ######

#' bGWAS - main function
#'
#' Performs a Bayesian GWAS from Summary Statistics ... work in progress to apply the
#' approach to Parent of Origin (PofO) results
#'
#'
#' @param name The name of the analysis (character)
#' @param GWAS_PofO
#' @param GWAS_Additive
#' @param noise (to be removed)
#' @param sign_method The method used to identify significant SNPs, should be \code{"p"} for
#'        p-value or \code{"fdr"} for false discovery rate, \code{default="p"} (character)
#' @param sign_thresh The threshold used to identify significant SNPs, \code{default="5e-8"}
#'        (numeric)
#' @param use_permutations  A logical indicating if BF p-values should be estimated using the permutation
#'        approach,  \code{default=FALSE}
#' @param res_pruning_dist The distance used for pruning results (in Kb), should be between 10 and 1000,
#'        (if set to NULL, no pruning is done), \code{default=500} (numeric)
#' @param res_pruning_LD The LD threshold used for pruning results, should be between 0 and 1
#'        (if 0, distance-based pruning is used), \code{default=0} (numeric)
#' @param save_files A logical indicating if the results should be saved as files,
#'        \code{default=FALSE}
#' @param verbose  A logical indicating if information on progress should be reported,
#'        \code{default=TRUE}
#'
#'          
#' @import dplyr
#' @import magrittr
#' @importFrom rlang .data
#' 
#' @export



bGWAS <- function(name,
                  GWAS_PofO,
                  GWAS_additive,
                  noise=0,
                  n_permutations = 1000,
                  save_files = TRUE,
                  verbose = TRUE) {
  
  
  # Path where the analysis has been launched
  InitPath = getwd()
  on.exit(setwd(InitPath))
  
  StartTime =  proc.time()
  
  ## Be chatty ?
  if(!is.logical(verbose)) stop("verbose : should be logical", call. = FALSE)
  
  
  # initialization of log_info file
  log_info = c()
  
  tmp = "<<< Preparation of analysis >>> \n"
  log_info = update_log(log_info, tmp, verbose)
  
  ### check the parameters ###
  tmp = c("", "> Checking parameters \n")
  log_info = update_log(log_info, tmp, verbose)
  
  ## Name of analysis
  if(!is.character(name)) stop("name : non-character argument", call. = FALSE) # should be a string
  
  tmp = paste0("The name of your analysis is: \"", name, "\". \n")
  log_info = update_log(log_info, tmp, verbose)
  

  #  if the directory already exists : error (only is save_files==T)
  if(!is.logical(save_files)) stop("save_files : should be logical", call. = FALSE)
  if(save_files){
    ### the directory to store the results ###
    Dir = file.path(InitPath, name)
    
    if(dir.exists(Dir)) stop("You already run an analysis with the same name in that directory,
                             please specify another name or choose another directory to run the analysis", call. = FALSE)
  }
  
  # Tidy input GWAS
  Data = tidy_inputGWAS(GWAS_PofO, GWAS_additive, verbose)
  
  tmp = paste0("The analysis will be run in the folder: \"", InitPath, "\".  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  ## save_files
  if(save_files){
    tmp = paste0("Files will be saved in: \"", file.path(InitPath,  name), "\".  \n")
    log_info = update_log(log_info, tmp, verbose)
  }
  
  # Go into the analysis' directory
  if(save_files){
    dir.create(Dir)
    setwd(Dir)
  }
  
  # 1 : Compute Prior
  log_info = update_log(log_info, c("", ""), F)
  tmp = c("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n", 
          "<<< Estimation of the prior >>>  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  tmp = c("", "> Computing prior for null paternal effects  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  Prior_pat0 = compute_prior_PofO(Data, "pat0", noise, verbose)
  log_info = c(log_info, Prior_pat0$log_info)
  
  tmp = c("", "> Computing prior for null maternal effects  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  Prior_mat0 = compute_prior_PofO(Data, "mat0", noise, verbose)
  log_info = c(log_info, Prior_mat0$log_info)
  
  ##### COMPUTE THE BAYES FACTOR AND THE P-VALUE #####
  log_info = update_log(log_info, c("", ""), F)
  tmp = c("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n",
          "<<< Calculation of Bayes Factors and p-values >>>  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  
  # This script create a file containing all SNPs in common between prior file / imputed files
  tmp = c("", "> Calculating them for null maternal effects  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  PriorWithBF_mat0 = request_BFandP_PofO(Prior_mat0$prior, pa, n_permutations, save_files, verbose)
  log_info = update_log(log_info, PriorWithBF_mat0$log_info, F)
  
  tmp = c("", "> Calculating them for null paternal effects  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  PriorWithBF_pat0 = request_BFandP_PofO(Prior_pat0$prior, "pat0", n_permutations, save_files, verbose)
  log_info = update_log(log_info, PriorWithBF_pat0$log_info, F)
  
  
  ##### IDENTIFY SIGNIFICANT SNPS + PRUNING #####
  # tmp = paste0("", "> Pruning and identifying significant SNPs \n")
  # log_info = update_log(log_info, tmp, verbose)
  # 
  # 1) we need to get correlation between the two sets of BFs to adjust the sign_thresh
  # 
  # 2a) get significant SNPs for null maternal effects
  # 2b) get significant SNPs for null paternal effects
  #
  # 3) bind them into a single data.frame, with a column for the model
  #
  # Results = get_significantSNPs(PriorWithBF$SNPs, sign_method, sign_thresh, res_pruning_dist, res_pruning_LD, 
  #                               res_MR$studies, matrix_all$mat, save_files, verbose)
  # log_info = update_log(log_info, Results$log_info, F)
 
  
  ### go back to inital folder ###
  log_info = update_log(log_info, c("", ""), F)
  tmp = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"
  log_info = update_log(log_info, tmp, verbose)
  
  
  ### write log_info File ###
  Time = as.integer((proc.time()-StartTime)[3])
  minutes <- as.integer(trunc(Time/60))
  seconds <- Time - minutes * 60
  tmp = paste0("Time of the analysis: ", minutes, " minute(s) and ", seconds, " second(s).  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  log_info = apply(as.array(log_info), 1,function(x) gsub("\n", "", x, fixed=T))
  
  if(save_files) write(log_info, paste0(name,".log"))

  
  results=list()
  results$log_info = log_info
  results$all_BFs_mat0 = PriorWithBF_mat0$SNPs
  results$all_BFs_pat0 = PriorWithBF_pat0$SNPs

  
  class(results) = "bGWAS"
  return(results)
  }
