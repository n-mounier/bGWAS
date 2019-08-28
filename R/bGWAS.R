###### Main Function ######

#' bGWAS - main function
#'
#' Performs a Bayesian GWAS from Summary Statistics, using publicly available results
#' to calculate the prior effects of the SNPs and compare it to observed z-scores
#'
#'
#' @param name The name of the analysis (character)
#' @param GWAS The path to the conventional GWAS of interest, the ID of the GWAS from the
#'        list of studies available (prior GWASs), or a \code{data.frame} (character, numeric or \code{data.frame})
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @param prior_studies The IDs of prior GWASs to use for the analysis, \code{default=NULL},
#'        will include all the prior GWASs available (numeric vector)
#' @param MR_threshold The threshold used to select strong instruments for MR, should be lower
#'        than 1e-5, \code{default=1e-6} (numeric)
#' @param MR_ninstruments The minimum number of strong instruments needed to use a prior GWAS, 
#'        should be between 2 and 8, \code{default=3} (numeric)
#' @param MR_pruning_dist The distance used for pruning MR instruments (in Kb), should be between 10 and 1000,
#'        \code{default=500} (numeric)
#' @param MR_pruning_LD The LD threshold used for pruning MR instruments, should be between 0 and 1
#'        (if 0, distance-based pruning is used), \code{default=0} (numeric)
#' @param prior_shrinkage The p-value threshold used for shrinkage before calculating the prior,
#'        should be between \code{MR_threshold} and 1, \code{default=1e-5} (numeric)
#' @param MR_shrinkage The p-value threshold used for shrinkage before performing MR, should be between
#'        \code{MR_threshold} and 1 (no shrinkage), \code{default=1} (numeric)
#'        #' @param prior_shrinkage The p-value threshold used for shrinkage before calculating the prior,
#'        should be between \code{MR_threshold} and 1, \code{default=1e-5} (numeric)
#' @param stepwise_threshold The p-value threshold used for inclusion/exclusion of Prior GWASs during the
#'        stepwise selection approach, should be between 0.05 and 0.0005, \code{default=NULL} will use 0.05 
#'        divided by the number of Prior GWASs tested (numeric)
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
#' @details
#' \code{Name} and \code{GWAS} are required arguments.
#' If \code{GWAS} is a path to a file (regular or .gz) or a \code{data.frame}, it should contain the following
#' columns : \cr
#' SNPID (rs numbers) should be : \code{rs}, \code{rsid}, \code{snp}, \code{snpid}, \code{rnpid} \cr
#' A1 should be : \code{a1}, \code{alt}, \code{alts} \cr
#' A2 should be : \code{a2}, \code{a0}, \code{ref} \cr
#' Z should be : \code{z}, \code{Z}, \code{zscore} \cr
#' If Z is not present, it can be calculated from BETA and SE. \cr
#' BETA should be : \code{b}, \code{beta}, \code{beta1} \cr
#' SE should be : \code{se}, \code{std} \cr
#' Note: in order to get rescaled (prior/posterior/corrected) effects, BETA and SE should be provided.
#' 
#' Z-Matrix files, containing Z-scores for all prior GWASs should be downloaded separately
#' and stored in \code{"~/ZMatrices"} or in the folder specified with the argument
#' \code{Z_matrices}. \cr
#' See [here](https://github.com/n-mounier/bGWAS) for more informations.
#'
#' Use \code{\link{list_priorGWASs}()} to see all the prior GWASs available.
#' Using one of them as your conventionnal GWAS (argument \code{GWAS} = numeric ID) will automatically
#' remove it from the list of prior GWASs used to build the prior.
#'
#' Use \link{select_priorGWASs}() to automatically select the prior GWASs to
#' be included/excluded when building the prior (argument \code{prior_studies}).
#'
#'
#' @return \code{bGWAS}() returns an object of class "bGWAS". \cr
#' Additionnaly, if \code{save_files=T}, several files are created in the folder \code{./name/} :
#' \itemize{
#' \item "PriorGWASs.tsv" - contains information about all prior GWASs
#' (general info + status (used/excluded) + MR coefficients)
#' \item "CoefficientsByChromosome.csv" - contains the MR estimates when masking the focal
#' chromosome (22 coefficients / prior GWASs used for prior estimation)
#' \item "PriorBFp.csv" - contains BF and p-values, prior, posterior and direct effects estimates for all SNPs
#' \item "SignificantSNPs.csv" - contains BF and p-values, prior, posterior and direct effects estimates for 
#' a subset of significant SNPs
#' }
#'
#' @examples
#' # Permorm bGWAS, using a small conventional GWAS included in the package (data.frame) 
#' # and selecting a subset of studies for the prior
#'\dontrun{top
#' data("SmallGWAS_Timmers2019")
#' MyStudies = select_priorGWASs(include_traits=c("Blood Pressure", "Education"),  
#'                               include_files=c("cardiogram_gwas_results.txt", 
#'                                              "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz"))
#' # 6 Prior GWASs used
#' list_priorGWASs(MyStudies) 
#'
#'  A = bGWAS(name="Test_UsingSmallDataFrame",
#'           GWAS = SmallGWAS_Timmers2019,
#'           prior_studies=MyStudies,
#'           stepwise_threshold=0.05,
#'           save_files=T)
#'           }
#'           
#' # Permorm bGWAS, using a conventional GWAS from the list of prior GWASs
#' \dontrun{MyGWAS = 3
#' list_priorGWASs(MyGWAS)
#' # Coronary Artery Disease GWAS (CARDIoGRAM)
#' B = bGWAS(name = "Test_UsingGWASfromPriorGWASs",
#'          GWAS = MyGWAS)
#'          }
#'          
#'          
#' @import dplyr
#' @import magrittr
#' @importFrom rlang .data
#' 
#' @export



bGWAS <- function(name,
                  GWAS,
                  Z_matrices = "~/ZMatrices/",
                  prior_studies = NULL,
                  MR_threshold = 1e-6,
                  MR_ninstruments = 3,
                  MR_pruning_dist = 500,
                  MR_pruning_LD = 0,
                  MR_shrinkage = 1,
                  stepwise_threshold = NULL,
                  prior_shrinkage = 1,
                  sign_method = "p",
                  sign_thresh = 5e-8,
                  use_permutations= FALSE,
                  res_pruning_dist = 500,
                  res_pruning_LD = 0,
                  save_files = FALSE,
                  verbose = TRUE) {
  
  
  # Path where the analysis has been launched
  InitPath = getwd()
  on.exit(setwd(InitPath))
  
  StartTime =  proc.time()
  
  # platform identification 
  # if platform is windows, an path should be provided for Z-matrices
  platform = .Platform$OS.type
  if(platform=="windows" && Z_matrices=="~/ZMatrices/") stop("Windows operating system, please provide a file path to the \"Z_matrices\" argument", call. = FALSE)
  
  
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
  if(save_files){
    ### the directory to store the results ###
    Dir = file.path(InitPath, name)
    
    if(dir.exists(Dir)) stop("You already run an analysis with the same name in that directory,
                             please specify another name or choose another directory to run the analysis", call. = FALSE)
  }
  
  
  ## ZMatrices
  if (is.character(Z_matrices)){
    # get absolute path
    Z_matrices = normalizePath(Z_matrices)
    
    if(!file.exists(file.path(Z_matrices, "ZMatrix_Full.csv.gz"))) stop("No \"ZMatrix_Full.csv.gz\" file in specified Z_matrices folder", call. = FALSE)
    if(!file.exists(file.path(Z_matrices, "ZMatrix_MR.csv.gz"))) stop("No \"ZMatrix_MR.csv.gz\" file in specified Z_matrices folder", call. = FALSE)
    if(!file.exists(file.path(Z_matrices, "AvailableStudies.tsv"))) stop("No \"AvailableStudies.tsv\" file in specified Z_matrices folder", call. = FALSE)
    

  } else stop("Z_matrices : wrong format, should be character", call. = FALSE)
  
  # use this to throw an error if not using the most recent Zmat files
  Prior_GWASs = list_priorGWASs(Z_matrices=Z_matrices)
  
  # fread deals with .gz files now
  MR_Files =  colnames(data.table::fread(file.path(Z_matrices, "ZMatrix_MR.csv.gz"), nrows=0, showProgress = FALSE))[-(1:5)]
  Prior_Files =  colnames(data.table::fread(file.path(Z_matrices, "ZMatrix_Full.csv.gz"), nrows=0, showProgress = FALSE))[-(1:5)]
  List_Files = list_files(Z_matrices=Z_matrices)
  if(!all(MR_Files==Prior_Files)) stop("Z_matrices : columns unconsistent between the two files")
  if(!all(MR_Files==List_Files))  stop("Z_matrices : columns unconsistent with \"AvailableStudies.tsv\"")
  
  tmp = paste0("The Z-Matrix files are stored in \"", Z_matrices, "\".  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  
  # Tidy input GWAS
  Data = tidy_inputGWAS(GWAS, Z_matrices, verbose)
  log_info = update_log(log_info, Data$log_info, F)
  
  
  
  
  tmp = paste0("The analysis will be run in the folder: \"", InitPath, "\".  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  ## save_files
  if(!is.logical(save_files)) stop("save_files should be logical", call. = FALSE)
  if(save_files){
    tmp = paste0("Files will be saved in: \"", file.path(InitPath,  name), "\".  \n")
    log_info = update_log(log_info, tmp, verbose)
  }
  
  
  ## prior_studies
  # check that all the files required exist in our list of studies
  # should be specified as "File ID"
  if(is.null(prior_studies)) prior_studies = c(1:length(list_files(Z_matrices = Z_matrices)))
  if(!all(prior_studies %in% c(1:length(list_files(Z_matrices = Z_matrices))))) stop("prior_studies : all the IDs provided should belong to the ones available", call. = FALSE)
  # if GWAS from data, make sure to remove it + the studies using the same trait!!!
  if(is.numeric(GWAS) && GWAS %in% prior_studies){
    # remove GWAS of interest
    prior_studies = prior_studies[prior_studies != GWAS]
    tmp = paste0("The study ", list_files(IDs=GWAS, Z_matrices = Z_matrices), " (ID=", GWAS, ") has been removed from the prior GWASs used to build the prior since it is used as conventionnal GWAS. \n")
    log_info = update_log(log_info, tmp, verbose)
  }
  if(is.numeric(GWAS)){
    # remove GWAS(s) for the same trait if needed
    Info = list_priorGWASs(Z_matrices = Z_matrices)
    Info %>%
      filter(.data$ID==GWAS) %>%
      pull(.data$Trait) -> trait
    Info %>%
      filter(.data$Trait==trait) %>%
    if(GWAS %in% ToExclude) ToExclude=ToExclude[ToExclude!=GWAS]
    if(any(ToExclude %in% prior_studies)){
      ToExclude = ToExclude[ToExclude %in% prior_studies]
      nbToExclude = length(ToExclude)
      prior_studies = prior_studies[!prior_studies %in% ToExclude]
      removed = ToExclude
      if(nbToExclude==1){
        tmp = paste0("The study ", list_files(IDs=ToExclude, Z_matrices = Z_matrices), " has been removed from the prior GWASs used to build the prior since it's a GWAS for the same trait as the study used as conventionnal GWAS. \n")
      }else{
        tmp = paste0("The studies : ", paste0(list_files(IDs=ToExclude, Z_matrices = Z_matrices), collapse=" - "), " have been removed from the prior GWASs used to build the prior since they are GWASs for the same trait as the study used as conventionnal GWAS. \n")
      }
      log_info = update_log(log_info, tmp, verbose)
    }
    
    # check that the user did not use exactly the same study for GWAS and for Prior,
    # in this case the prior study is removed here... and we don't have any study left
    if(length(prior_studies)==0) stop("You did not select any prior GWAS different from your conventionnal GWAS")
    
  }
  
  
  ## MR_threshold -> should not be larger than 10-5, can only be more stringent
  if(!is.numeric(MR_threshold)) stop("MR_threshold : non-numeric argument", call. = FALSE)
  if(MR_threshold>10^5) stop("MR_threshold : superior to the threshold limit (10^-5)", call. = FALSE)
  
  tmp = paste0("The p-value threshold used for selecting MR instruments is: ", format(MR_threshold, scientific = T), ".  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  ## MR_ninstruments -> should not be larger than 10-5, can only be more stringent
  if(!is.numeric(MR_ninstruments)) stop("MR_ninstruments : non-numeric argument", call. = FALSE)
  if(MR_ninstruments>8) stop("MR_ninstruments : too large, should be between 2 and 8", call. = FALSE)
  if(MR_ninstruments<2) stop("MR_ninstruments : too small, should be between 2 and 8", call. = FALSE)
  
  tmp = paste0("The minimum number instruments required for each trait is: ", MR_ninstruments, ".  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  ## MR_pruning_dist
  if(!is.numeric(MR_pruning_dist)) stop("MR_pruning_dist : non-numeric argument", call. = FALSE)
  if(MR_pruning_dist<10) stop("MR_pruning_dist : should be higher than 10Kb", call. = FALSE)
  if(MR_pruning_dist>1000) stop("MR_pruning_dist : should be lower than 1Mb", call. = FALSE)
  
  
  tmp = paste0("The distance used for pruning MR instruments is: ", MR_pruning_dist, "Kb.  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  
  
  ## MR_pruning_LD
  if(!is.numeric(MR_pruning_LD)) stop("MR_pruning_LD : non-numeric argument", call. = FALSE)
  if(MR_pruning_LD<0) stop("MR_pruning_LD : should be positive", call. = FALSE)
  
  if(MR_pruning_LD>0){
    tmp = paste0("The LD threshold used for pruning MR instruments is: ",MR_pruning_LD, ".  \n")
    log_info = update_log(log_info, tmp, verbose)
  } else {
    tmp = "Distance-based pruning will be used for MR instruments.  \n"
    log_info = update_log(log_info, tmp, verbose)
  }
  
  ## MR_shrinkage, should be between MR_threshold and 1
  if(!is.numeric(MR_shrinkage)) stop("MR_shrinkage : non-numeric argument", call. = FALSE)
  if(MR_shrinkage<0) stop("MR_shrinkage : should be positive", call. = FALSE)
  if(MR_shrinkage<MR_threshold) stop("MR_shrinkage : should be higher than the threshold used to select MR instruments", call. = FALSE)
  if(MR_shrinkage>1) stop("MR_shrinkage : should not be higher than 1", call. = FALSE)
  
  if(MR_shrinkage==1){
    tmp = "No shrinkage applied before performing MR."
    log_info = update_log(log_info, tmp, verbose)
  } else {
    tmp = paste0("All p-values lower than ",format(MR_shrinkage, scientific = T), " will be shrunk to 0 before performing MR.  \n")
    log_info = update_log(log_info, tmp, verbose)
  }
  
  ## stepwise_threshold, should be between 0.05 and 0.0005 or NULL
  if(is.null(stepwise_threshold)){
    tmp = "The p-value threshold used for stepwise selection will be derived according to the number of Prior GWASs used.  \n"
    log_info = update_log(log_info, tmp, verbose)
  } else if(is.numeric(stepwise_threshold)){
    if(stepwise_threshold<0.0005) stop("stepwise_threshold : should not be lower than 0.0005", call. = FALSE)
    if(stepwise_threshold>0.05) stop("stepwise_threshold : should not be higher 0.05", call. = FALSE)
    tmp = paste0("The p-value threshold used for stepwise selection is ", format(stepwise_threshold, scientific = F), ".  \n")
    log_info = update_log(log_info, tmp, verbose)
  } else{
    stop("stepwise_threshold : should not numeric or NULL", call. = FALSE)
  }
  
  
  ## prior_shrinkage, should be between MR_threshold and 1
  if(!is.numeric(prior_shrinkage)) stop("prior_shrinkage : non-numeric argument", call. = FALSE)
  if(prior_shrinkage<0) stop("prior_shrinkage : should be positive", call. = FALSE)
  if(prior_shrinkage<MR_threshold) stop("prior_shrinkage : should be higher than the threshold used to select MR instruments", call. = FALSE)
  if(prior_shrinkage>1) stop("prior_shrinkage : should not be higher than 1", call. = FALSE)
  
  if(prior_shrinkage==1){
    tmp = "No shrinkage applied before performing calculating the prior."
    log_info = update_log(log_info, tmp, verbose)
  } else {
    tmp = paste0("All p-values lower than ", format(prior_shrinkage, scientific = T), " will be shrunk to 0 before calculating the prior.  \n")
    log_info = update_log(log_info, tmp, verbose)
  }
  
  
  ## sign_method -> should not be "p" or "fdr"
  if(!sign_method %in% c("p", "fdr")) stop("sign_method : method not accepted, should be p or fdr", call. = FALSE)
  
  ## sign_thresh -> should be numeric and lower (or equal) to 1
  if(!is.numeric(sign_thresh)) stop("sign_thresh : non numeric threshold", call. = FALSE)
  if(sign_thresh>1) stop("sign_thresh : a threshold higher than 1 does not make sense", call. = FALSE)
  
  if(sign_method=="p"){
    tmp = paste0("Significant SNPs will be identified according to p-value. The threshold used is :",
                 format(sign_thresh, scientific = T), ".  \n")
  } else if(sign_method=="fdr"){
    tmp = paste0("Significant SNPs will be identified according to FDR. The threshold used is :",
                 format(sign_thresh, scientific = T), ".  \n")
  }
  log_info = update_log(log_info, tmp, verbose)
  
  ## res_pruning_dist
  if(!is.numeric(res_pruning_dist)) stop("res_pruning_dist : non-numeric argument", call. = FALSE)
  
  if(is.null(res_pruning_dist)){
    tmp = "Results will not be pruned.  \n"
    log_info = update_log(log_info, tmp, verbose)
  } else {
    if(res_pruning_dist<10) stop("res_pruning_dist : should be higher than 10Kb", call. = FALSE)
    if(res_pruning_dist>1000) stop("res_pruning_dist : should be lower than 1Mb", call. = FALSE)
    
    
    tmp = paste0("The distance used for pruning results is: ", res_pruning_dist, "Kb.  \n")
    log_info = update_log(log_info, tmp, verbose)
    
    ## res_pruning_LD
    if(!is.numeric(res_pruning_LD)) stop("res_pruning_LD : non-numeric argument", call. = FALSE)
    if(res_pruning_LD<0) stop("res_pruning_LD : should be positive", call. = FALSE)
    
    if(res_pruning_LD>0){
      tmp = paste0("The LD threshold used for pruning results is: ", res_pruning_LD, ".  \n")
      log_info = update_log(log_info, tmp, verbose)
    } else {
      tmp = "Distance-based pruning will be used for results.  \n"
      log_info = update_log(log_info, tmp, verbose)
    }
  }
  
  
  # Go into the analysis' directory
  if(save_files){
    dir.create(Dir)
    setwd(Dir)
  }
  
  ### 1 : create "Summarize_file" ###
  
  if(save_files){
    tmp = paste0("# Initializing the summary information file\n")
    log_info = update_log(log_info, tmp, verbose)
    
    Files_Info = list_priorGWASs(Z_matrices = Z_matrices)
    # keep only interesting columns + add "Status" column
    Files_Info %>%
      select("File", "Name", "Trait") %>%
      mutate(status= "Exluded by user") -> Files_Info
    Files_Info[prior_studies, "status"] = "USED"
    if(is.numeric(GWAS))  Files_Info[GWAS, "status"] = "Conventionnal GWAS"
    
    utils::write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )
    
    tmp = paste0("List of files : ", Dir, "/PriorGWASs.csv has been successfully created.  \n")
    log_info = update_log(log_info, tmp, verbose)
  }
  
  # 2 : Z-Matrix for MR
  log_info = update_log(log_info, c("", ""), F)
  tmp = c("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n",
          "<<< Identification of significant prior GWASs for MR >>>  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  # We should keep the Z-Matrix creation outside of the study identification function
  # so that we can quickly re-run the second part using a file containing the Z-Matrix
  tmp = c("", "> Creating the Z-Matrix of strong instruments \n")
  log_info = update_log(log_info, tmp, verbose)
  
  # the "global z_matrix of strong instruments" for all GWASs already exists, just select the studies kept for the prior + prune + add the GWAS of interest
  # makeMR_ZMatrix() create a ZMatrix file and returns the log_info
  matrix_MR = makeMR_ZMatrix(prior_studies, Data$GWAS, Data$GWAS_Name, MR_threshold, MR_ninstruments, MR_pruning_dist, MR_pruning_LD, Z_matrices, MR_shrinkage, save_files, verbose)
  log_info = update_log(log_info, matrix_MR$log_info, F)
  # here, add potential stop() in function(s) and check for it?
  
  
  tmp = c("", "> Performing MR  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  res_MR = identify_studiesMR(matrix_MR$mat, MR_shrinkage, MR_threshold, stepwise_threshold, Z_matrices, save_files, verbose)
  log_info = c(log_info, res_MR$log_info)
  # if error/problem in identify_studiesMR (i.e no study significant)
  if(isTRUE(res_MR$stop)){
    
    log_info = apply(as.array(log_info), 1,function(x) gsub("\n", "", x, fixed=T))
    
    if(save_files) write(log_info, paste0(name,".log"))
    
    print("Analysis stopped : see log for more informations.")
    
    results=list()
    results$log_info = log_info
    class(results) = "bGWAS"
    return(results)
  }
  
  
  # 3 : make MR ZMat
  # 4 : Compute Prior
  log_info = update_log(log_info, c("", ""), F)
  tmp = c("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n", 
          "<<< Estimation of the prior >>>  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  tmp = c("", "> Creating the full Z-Matrix  \n")
  log_info = update_log(log_info, tmp, verbose)
  Studies = select_priorGWASs(include_files=res_MR$studies, Z_matrices = Z_matrices)
  matrix_all = makeFull_ZMatrix(Studies, Data$GWAS,  Data$GWAS_Name, Z_matrices, prior_shrinkage, save_files, verbose)
  log_info = update_log(log_info, matrix_all$log_info, F)
  
  tmp = c("", "> Computing prior  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  Prior = compute_prior(res_MR$studies, res_MR$ZMat, matrix_all$mat, Data$GWAS, Data$rescaling, MR_shrinkage, prior_shrinkage, Z_matrices, save_files, verbose)
  log_info = c(log_info, Prior$log_info)
  
  
  
  ##### COMPUTE THE BAYES FACTOR AND THE P-VALUE #####
  log_info = update_log(log_info, c("", ""), F)
  tmp = c("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n",
          "<<< Calculation of Bayes Factors and p-values >>>  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  
  # This script create a file containing all SNPs in common between prior file / imputed files
  tmp = c("", "> Calculating them for all SNPs  \n")
  log_info = update_log(log_info, tmp, verbose)
  
  PriorWithBF = request_BFandP(Prior$prior, sign_thresh, use_permutations, sign_method, save_files, verbose)
  log_info = update_log(log_info, PriorWithBF$log_info, F)
  
  
  
  ##### IDENTIFY SIGNIFICANT SNPS + PRUNING #####
  tmp = paste0("", "> Pruning and identifying significant SNPs \n")
  log_info = update_log(log_info, tmp, verbose)
  
  
  Results = get_significantSNPs(PriorWithBF$SNPs, sign_method, sign_thresh, res_pruning_dist, res_pruning_LD, 
                                res_MR$studies, matrix_all$mat, save_files, verbose)
  log_info = update_log(log_info, Results$log_info, F)
  
  
  
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
  results$significant_SNPs = Results$SNPs
  results$all_BFs = PriorWithBF$SNPs
  results$significant_studies = res_MR$coeffs
  results$all_MRcoeffs = Prior$all_coeffs
  results$matrix_heatmap = Results$mat
  
  class(results) = "bGWAS"
  return(results)
  }
