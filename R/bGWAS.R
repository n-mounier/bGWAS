###### Main Function ######



#' bGWAS
#'
#' Performs a bayesian GWAS from Summary Statistics, using publicly available results
#' to calculate the prior effects of the SNPs and compare it to observed z-scores
#'
#'
#' @param name The name of the analysis (character)
#' @param GWAS The path to the conventionnal GWAS of interest, the ID of the GWAS from the
#'        list of studies availables (prior GWASs), or a data.frame (character, numeric or data.frame)
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @param prior_studies The IDs of prior GWASs to use for the analysis, \code{default=NULL},
#'        will include all the prior GWASs available (numeric vector)
#' @param SNPs_list, The path to a file containing the rsids of the SNPs to use ,
#'        \code{default=NULL}, will use all the SNPs in common between prior GWASs and the
#'        conventionnal GWAS (character)
#'              # NOT IMPLEMENTED YET
#' @param MR_threshold The threshold used to select strong instruments for MR, should be lower
#'        than 1e-5, \code{default=1e-5} (numeric)
#' @param sign_method The method used to identify significant SNPs, should be \code{"p"} for
#'        p-value or \code{"fdr"} for false discovery rate, \code{default="p"} (character)
#' @param sign_thresh The threshold used to identify significant SNPs, \code{default="5.10e-8"}
#'        (numeric)
#' @param prune_res A logical indicating if the results should be pruned (by distance, 500kb),
#'        \code{default=FALSE}
#' @param save_files A logical indicating if the results should be saved as files,
#'        \code{default=FALSE}
# #' @param OutPath character, path to the outputs, needed if saveFiles is TRUE, by default,
# #' current working dictory
#' @param verbose  A logical indicating if information on progress should be reported,
#'        \code{default=FALSE}
#'
#' @details
#' \code{Name} and \code{GWAS} are required arguments.
#' If \code{GWAS} is a path to a file (regular or .gz) or a data.frame, it should contain the following
#' columns : \cr
#' SNPID (rs numbers) should be : \code{rs}, \code{rsid}, \code{snp}, \code{snpid}, \code{rnpid} \cr
#' A1 should be : \code{a1}, \code{alt}, \code{alts} \cr
#' A2 should be : \code{a2}, \code{a0}, \code{ref} \cr
#' Z should be : \code{z}, \code{Z}, \code{zscore} \cr
#' If Z is not present, it can be calculated from BETA and SE, in this case, a temporary
#' gzipped file is created and removed after the analysis. \cr
#' BETA should be : \code{b}, \code{beta}, \code{beta1} \cr
#' SE should be : \code{se}, \code{std} \cr
#' Z-Matrix files, containing Z-scores for all prior GWASs should be downloaded separately
#' and stored in \code{"~/ZMatrices"} or in the folder specified with the argument
#' \code{Z_matrices}. \cr
#' ## MORE INFO NEEDED HERE / HOW TO DOWNLOAD ID (see GitHub README)
#'
#' Use \code{\link{list_priorGWASs}()} to see all the prior GWASs available.
#' Using one of them as your conventionnal GWAS (argument \code{GWAS} = numeric ID) will automatically
#' remove it from the list of prior GWASs used to build the prior.
#'
#' Use \link{select_priorGWASs}() to automatically select the studies to
#' be included/excluded when building the prior (argument \code{prior_studies}).
#'
#' Be careful, in the results, all the GWAS (conventional + prior) are aligned with UK10K
#' data for the analysis (some alleles might be swapped when comparing with the initial data)
#'
#' @return \code{bGWAS}() returns an object of class "bGWAS". \cr
#' Aditionnaly, if \code{save_files=T}, several files are created... \cr
#' ... in your working directory :
#' \itemize{
#' \item  "\code{name}.log" - log file
#' }
#' ... in the folder \code{./name/} :
#' \itemize{
#' \item "PriorGWASs.tsv" - contains Prior GWASs information
#' (general info + status (used/removed) + MR coefficients)
#' \item "CoefficientsByChromosome.csv" - contains the MR estimates when masking the focal
#' chromosome (22 coefficients / study)
#' \item "PriorBFp.csv" - contains BF and p-values estimated for all SNPs
#' \item "SignificantSNPs.csv" - contains BF and p-values estimated for a subset of SNPs
#' }
#'
#' @examples
#' # Permorm bGWAS, using a conventional GWAS from the list of prior GWASs
#' MyGWAS = 5
#' list_priorGWASs(MyGWAS)
#' \dontrun{
#' A = bGWAS(name = "Test_UsingGWASfromPriorGWASs",
#'          GWAS = MyGWAS,
#'          verbose=T)
#'          }
#'
#'# Permorm bGWAS, using a small conventional GWAS included in data (file) and selecting a subset of
#'# studies for the prior
#' MyGWAS = system.file("Data/SmallGWAS_Pilling2017.csv", package="bGWAS")
#' MyStudies = select_priorGWASs(include_traits=c("Type 2 diabetes", "Smoking"),
#'                          include_files=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))
#' \dontrun{
#' B = bGWAS(name = "Test_UsingSmallGWAS",
#'           GWAS = MyGWAS,
#'           prior_studies=MyStudies,
#'           verbose=T)
#'          }
#'
#'#'# Permorm bGWAS, using a small conventional GWAS included in data (data.frame) and selecting a subset of
#'# studies for the prior
#'\dontrun{
#' data("SmallGWAS_Pilling2017")
#' C = bGWAS(name="Test_UsingSmallDataFrame",
#'           GWAS = SmallGWAS_Pilling2017,
#'           prior_studies=MyStudies,
#'           verbose=T,
#'           save_files=T)
#'           }
#'
#'# Note that B and C are using the same data (stored differently) and give the same results.

#' @export



bGWAS <- function(name,
                  GWAS,
                  Z_matrices = NULL,
                  prior_studies = NULL,
                  SNPs_list = NULL,
                  MR_threshold = 1e-5,
                  sign_method = "p",
                  sign_thresh = 5e-8,
                  prune_res = FALSE,
#                  OutPath = getwd(),
                  save_files = FALSE,
                  verbose = FALSE) {


  InitPath = getwd()
### TO BE DONE
  OutPath = getwd() # to be cleaned
  # depends if we allow the user to run the analysis in another directory
  # of if we force the use of getwd()

  StartTime =  proc.time()

  # platform identification : used in the main function
  # automatically re-detected when needed by other sub-functions
  platform = c("Linux", "macOS", "W")[c(grepl("Linux", sessionInfo()$running)
                                        , grepl("macOS", sessionInfo()$running)
                                        , grepl("Windows", sessionInfo()$running))]
  if(platform=="W") stop("Windows is not supported yet")

  # initialization of log_info file
  log_info = c()

  tmp = paste0("<<< Preparation of analysis >>> \n")
  log_info = update_log(log_info, tmp, verbose)

  ### check the parameters ###
  tmp = paste0("> Checking parameters \n")
  log_info = update_log(log_info, tmp, verbose)


  ## Be chatty ?
  if(!is.logical(verbose)) stop("verbose : should be logical")


  ## Name of analysis
  if(!is.character(name)) stop("name : non-character argument") # should be a string

  tmp = paste0("The name of your analysis is: \"", name, "\". \n")
  log_info = update_log(log_info, tmp, verbose)

  ### create the directory to store the results ###
### TO BE DONE
  # only if save_file=T
  Dir = file.path(OutPath, name)

  #  if the directory already exists : error
  if(file.exists(file.path(OutPath, paste0(name, ".log")))) stop("You already run an analysis with the same name in that directory,
                                                    please specify another name or choose another directory to run the analysis")
  if(save_files){
    ifelse(!dir.exists(Dir), dir.create(Dir), stop("You already run an analysis with the same name in that directory,
                                                 please specify another name or choose another directory to run the analysis"))
  }


  ## ZMatrices
  if(is.null(Z_matrices)){
    if(!file.exists("~/ZMatrices/ZMatrix_Imputed.csv.gz")) stop("No ZMatrix_Imputed.csv.gz file in ~/ZMatrices")
    if(!file.exists("~/ZMatrices/ZMatrix_NotImputed.csv.gz")) stop("No ZMatrix_NotImputed.csv.gz file in ~/ZMatrices")
    Z_matrices = "~/ZMatrices"
  } else if (is.character(Z_matrices)){
    if(!file.exists(paste0(Z_matrices, "/ZMatrix_Imputed.csv.gz"))) stop("No ZMatrix_Imputed.csv.gz file in specified Z_matrices folder")
    if(!file.exists(paste0(Z_matrices, "/ZMatrix_NotImputed.csv.gz"))) stop("No ZMatrix_NotImputed.csv.gz file in specified Z_matrices folder")
    # Define if the path is relative or not
    if(!substr(Z_matrices,1,1) %in% c("/", "~")){
      Z_matrices=paste0(InitPath, "/", Z_matrices)
    }
  } else stop("Z_matrices : wrong format")

  tmp = paste0("The Z-Matrix files are stored in \"", Z_matrices, "\".  \n")
  log_info = update_log(log_info, tmp, verbose)




  ## GWAS of interest, should be a path to a GWAS file (format ? .tar.gz or file ?), a data.frame or an ID
  TMP_FILE = F # flag : is a temporary file with Z-scores created ??
  if(is.numeric(GWAS)) { # if it is an ID
    if(!GWAS %in% c(1:length(list_files()))) stop("The ID specified as a conventional GWAS is not in the list")
    tmp = paste0("The conventional GWAS used as input is:",
                 list_files(IDs=GWAS), " (ID = ",  GWAS,").  \n")
    log_info = update_log(log_info, tmp, verbose)
  } else if(is.character(GWAS)) { # if it is a file
    # First, does the file exists ?
    if(!file.exists(GWAS)) stop("GWAS : the file does not exist")
    # Then, check if it is an absolute or a relative path to the file
    if(!substr(GWAS,1,1) %in% c("/", "~")){ # We should work with absolute path to avoid errors
      GWAS=paste0(InitPath, "/", GWAS)
    }
    tmp = paste0("The conventional GWAS used as input is: \"",
                 strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n")
    log_info = update_log(log_info, tmp, verbose)

    # Check colnames...
    if(!grepl(".gz", GWAS)){
      # ...if regular file
      HeaderGWAS = colnames(data.table::fread(GWAS, nrows = 0, showProgress = FALSE))
    } else if(grepl(".gz", GWAS)) {
      # ...if tar.gz
      if(platform %in% c("Linux", "macOS")) HeaderGWAS = colnames(data.table::fread(paste0("zcat < ", GWAS), nrows = 0, showProgress = FALSE))
    }

    if(all(!HeaderGWAS %in% c("rsid", "snpid", "snp", "rnpid", "rs"))) stop("GWAS : no SNPID column")
    # how to deal with multiple rsid / snpid columns ???
    # here, we don't care, we need at least one
    tmp = paste0("   SNPID column, ok")
    if(all(!HeaderGWAS %in% c("a1", "alts", "alt"))) stop("GWAS : no ALT column")
    tmp = c(tmp, paste0("ALT column, ok"))
    if(all(!HeaderGWAS %in% c("a2", "a0", "ref"))) stop("GWAS : no REF column")
    tmp = c(tmp, paste0("REF column, ok"))
    if(all(!HeaderGWAS %in% c("z", "Z", "zscore"))){
      # allow for beta + se to calculate Z ???
      if(!all(!HeaderGWAS %in% c("b", "beta", "beta1")) & !all(!HeaderGWAS %in% c("se", "std"))){
        # if beta + se : read the data
        if(!grepl(".gz", GWAS)){
          # ...if regular file
          DataGWAS = data.table::fread(GWAS, showProgress = FALSE)
        } else if(grepl(".gz", GWAS)) {
          # ...if tar.gz
          if(platform %in% c("Linux", "macOS"))  DataGWAS = data.table::fread(paste0("zcat < ", GWAS), showProgress = FALSE)
        }
        # calculate Z
        DataGWAS$Z = DataGWAS[,HeaderGWAS[HeaderGWAS %in% c("b", "beta", "beta1")], with=F] /
          DataGWAS[,HeaderGWAS[HeaderGWAS %in% c("se", "std")], with=F]
        # keep only relevant column to save space
        DataGWAS = DataGWAS[,c(HeaderGWAS[HeaderGWAS %in% c("rsid", "snpid", "snp", "rnpid", "rs")][1], # SNPID, order by relevance :
                               # if there is an "rsid" and a "snpid" -> choose "rsid" !
                               # order potential column names from the most likely to
                               # the least likely, so we chose the "right one" here
                               HeaderGWAS[HeaderGWAS %in% c("a1", "alts")], # ALT
                               HeaderGWAS[HeaderGWAS %in% c("a2", "a0", "ref")], # REF
                               "Z"), # Z
                            with=F]
        # write the data (as tar.gz) and change GWAS name to the created file
        # but save it in the current folder, no initial one
        TMP_Name = paste0(gsub(".gz", "",  paste0(getwd(), "/", strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])])), "_withZ.gz")
        write.table(DataGWAS, file=gzfile(TMP_Name), sep="\t",
                    quote=F, row.names=F)
        GWAS = TMP_Name
        # flag the created file
        TMP_FILE = T
      } else {
        stop("GWAS : no Z-SCORE column")
      }
    } else {
      tmp = c(tmp, paste0("Z column, ok  \n"))
    }
    tmp = paste(tmp, collapse= " - ")
    log_info = update_log(log_info, tmp, verbose)

    if(TMP_FILE){
      tmp = paste0("A temporary file with a Z column has been created : \"",
                   strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n")
      log_info = update_log(log_info, tmp, verbose)
    }
  } else if(is.data.frame(GWAS)){ # if data.frame
    # add attribute GName to the data.frame, to be re-used in other subfunctions
    attributes(GWAS)$GName =  deparse(substitute(GWAS)) # get the "name" of the object used as an argument in the function
    # transform into data.table (needed when creating ZMat for data manipulation)
    GWAS = data.table::as.data.table(GWAS)
    tmp = paste0("The conventional GWAS used as input the object: \"",
                 attributes(GWAS)$GName, "\".  \n")
    log_info = update_log(log_info, tmp, verbose)

    HeaderGWAS = colnames(GWAS)
    if(all(!HeaderGWAS %in% c("rsid", "snpid", "snp", "rnpid", "rs"))) stop("GWAS : no SNPID column")
    # how to deal with multiple rsid / snpid columns ???
    # here, we don't care, we need at least one
    tmp = paste0("   SNPID column, ok")
    if(all(!HeaderGWAS %in% c("a1", "alts", "alt"))) stop("GWAS : no ALT column")
    tmp = c(tmp, paste0("ALT column, ok"))
    if(all(!HeaderGWAS %in% c("a2", "a0", "ref"))) stop("GWAS : no REF column")
    tmp = c(tmp, paste0("REF column, ok"))
    if(all(!HeaderGWAS %in% c("z", "Z", "zscore"))){
      # allow for beta + se to calculate Z
      if(!all(!HeaderGWAS %in% c("b", "beta", "beta1")) & !all(!HeaderGWAS %in% c("se", "std"))){
        # if beta + se : calculate Z
        GWAS$Z = GWAS[,HeaderGWAS[HeaderGWAS %in% c("b", "beta", "beta1")], with=F] /
          GWAS[,HeaderGWAS[HeaderGWAS %in% c("se", "std")], with=F]
      } else {
        stop("GWAS : no Z-SCORE column")
      }
    } else {
      tmp = c(tmp, paste0("Z column, ok  \n"))
    }
    tmp = paste(tmp, collapse= " - ")
    log_info = update_log(log_info, tmp, verbose)
  } else stop("GWAS : unrecognized format")


### TO BE DONE
  ## OutPath, check that the directory exist. Create it if necessary ?
  if(is.null(OutPath)) OutPath = getwd()
  if(!dir.exists(OutPath)) stop("OutPath : the directory does not exist")

  tmp = paste0("The analysis will be run in the folder: \"", OutPath, "\".  \n")
  log_info = update_log(log_info, tmp, verbose)

  ## save_files
  if(!is.logical(save_files)) stop("save_files should be logical")
  if(save_files){
    tmp = paste0("Files will be saved in: \"", OutPath, "/", name, "\".  \n")
    log_info = update_log(log_info, tmp, verbose)
  }



### TO BE DONE
  ## ListOfSNPs
  # We should have at least XX SNPs / check Linux-MAC for Zcat
  # Also check the number or SNPs in common in the file if ListOfSNPs not specified ?
  if(!is.null(SNPs_list)){
    if(!is.character(SNPs_list)) stop("SNPs_list : non character")
    OurSNPsMR = data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=1, showProgress = FALSE)
    commonMR = table(ListOfSNPs %in% OurSNPsMR$rs)["TRUE"]
    if(commonMR<20) stop("ListOfSNPs : You should provide at least 20 SNPs that can be used as strong instruments for MR")
    OurSNPsAll = data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=1, showProgress = FALSE)
    commonAll = table(ListOfSNPs %in% OurSNPsAll$rs)["TRUE"]
    tmp = paste0(length(ListOfSNPs), " SNPs provided \n")
    tmp = c(tmp, paste0(commonMR, " can be used for MR \n"))
    tmp = c(tmp, paste0(commonAll, " can be used to compute prior \n"))
    log_info = update_log(log_info, tmp, verbose)
  }


  ## prior_studies
  # check that all the files required exist in our list of studies
  # should be specified as "File ID"
  if(is.null(prior_studies)) prior_studies = c(1:length(list_files()))
  if(!all(prior_studies %in% c(1:length(list_files())))) stop("prior_studies : all the IDs provided should belong to the ones available")
  # if GWAS from data, make sure to remove it
  if(is.numeric(GWAS) && GWAS %in% prior_studies){
    prior_studies = prior_studies[-GWAS]
    tmp = paste0("The study ", list_files(ID=GWAS), " (ID=", GWAS, ") has been removed from the studies used to build the prior since it is used as conventionnal GWAS. \n")
    log_info = update_log(log_info, tmp, verbose)
### TO BE DONE
    # check that the user did not use exactly the same study for GWAS and for Prior,
    # in this case the prior study is removed here... and we don't have any study left
  }


  ## MR_threshold -> should not be larger than 10-5, can only be more stringent
  if(!is.numeric(MR_threshold)) stop("MR_threshold : non-numeric argument")
  if(MR_threshold>10^5) stop("MR_threshold : superior to the threshold limit (10^-5)")

  tmp = paste0("The p-value threshold used for selecting MR instruments is: ", format(MR_threshold, scientific = T), ".  \n")
  log_info = update_log(log_info, tmp, verbose)


  ## sign_method -> should not be "p" or "fdr"
  if(!sign_method %in% c("p", "fdr")) stop("sign_method : method not accepted, should be p or fdr")

  ## sign_thresh -> should be numeric and lower (or equal) to 1
  if(!is.numeric(sign_thresh)) stop("sign_thresh : non numeric threshold")
  if(sign_thresh>1) stop("sign_thresh : a threshold higher than 1 does not make sense")

  if(sign_method=="p"){
    tmp = paste0("Significant SNPs will be identified according to p-value. The threshold used is :",
                 format(sign_thresh, scientific = T), ".  \n")
  } else if(sign_method=="fdr"){
    tmp = paste0("Significant SNPs will be identified according to FDR. The threshold used is :",
                 format(sign_thresh, scientific = T), ".  \n")
  }
  log_info = update_log(log_info, tmp, verbose)

  ## prune_res
  if(!is.logical(prune_res)) stop("prune_res : should be log_infoical")
  if(prune_res){
    tmp = "Significant SNPs will be pruned (distance pruning, 500kb) \n"
    log_info = update_log(log_info, tmp, verbose)
  }


  # Go into the analysis' directory
  if(save_files){
    setwd(Dir)
  }

  ### 1 : create "Summarize_file" ###

  if(save_files){
    tmp = paste0("# Initializing the summary information file\n")
    log_info = update_log(log_info, tmp, verbose)

    Files_Info = list_priorGWASs()
    # keep only interesting columns + add "Status" column
    Files_Info = Files_Info[, c(1,3:6)]
    Files_Info$Status= "Exluded by user"
    Files_Info[prior_studies, "Status"] = "USED"
    if(is.numeric(GWAS))  Files_Info[GWAS, "Status"] = "Conventionnal GWAS"

    write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )

    tmp = paste0("List of files : ", Dir, "/PriorGWASs.csv has been successfully created.  \n")
    log_info = update_log(log_info, tmp, verbose)
  }

  # 2 : Z-Matrix for MR
  log_info = c(log_info, "", "")
  tmp = paste0("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n",
               "<<< Identification of significant studies for MR >>>  \n")
  log_info = update_log(log_info, tmp, verbose)

  # We should keep the Z-Matrix creation outside of the study identification function
  # so that we can quickly re-run the second part using a file containing the Z-Matrix
  tmp = "> Creating the Z-Matrix of strong instruments \n"
  log_info = update_log(log_info, tmp, verbose)

  # the "global z_matrix" for all GWAS should already be done, just select the studies kept for the prior + prune + add the GWAS of interest
  # makeMR_ZMatrix() create a ZMatrix file and returns the log_info
  matrix_MR = makeMR_ZMatrix(prior_studies, GWAS, MR_threshold, Z_matrices, save_files, verbose)
  log_info = c(log_info,matrix_MR$log_info)


  log_info=c(log_info,"")
  tmp = paste0("> Performing MR  \n")
  log_info = update_log(log_info, tmp, verbose)

  res_MR = identify_studiesMR(matrix_MR$mat, save_files, verbose)
  log_info = c(log_info, res_MR$log_info)

  # 4 : Compute Prior
  log_info = c(log_info, "", "")
  tmp = paste0("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n",
                       "<<< Estimation of the prior >>>  \n")
  log_info = update_log(log_info, tmp, verbose)

  tmp = paste0("> Creating the full Z-Matrix  \n")
  log_info = update_log(log_info, tmp, verbose)
  Studies = select_priorGWASs(include_files=res_MR$studies$study_selected)
  matrix_all = makeFull_ZMatrix(Studies, GWAS, Z_matrices, save_files, verbose)
  log_info = c(log_info,matrix_all$log_info)

  tmp = paste0("> Computing prior  \n")
  log_info = update_log(log_info, tmp, verbose)

  Prior = compute_prior(res_MR$studies,matrix_MR$mat, matrix_all$mat, save_files, verbose)
  log_info = c(log_info, Prior$log_info)


  ##### COMPUTE THE BAYES FACTOR AND THE P-VALUE #####
  log_info = c(log_info, "", "")
  tmp = paste0("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n",
               "<<< Calculation of Bayes Factors and p-values >>>  \n")
  log_info = update_log(log_info, tmp, verbose)

  # Compute BFs for 100 nulls + for our GWAS of Interest
  # Calculate the p-values from these BFs (comparison with the nulls)
  # This script create a file containing all SNPs in common between prior file / imputed files
  tmp = paste0("> Calculating them for all SNPs  \n")
  log_info = update_log(log_info, tmp, verbose)

  PriorWithBF = request_BFandP(Prior$prior, save_files, verbose)
  log_info = c(log_info, PriorWithBF$log_info)



  ##### IDENTIFY SIGNIFICANT SNPS AT 5% FDR + PRUNING #####
  tmp = paste0("> Pruning and identifying significant SNPs \n")
  log_info = update_log(log_info, tmp, verbose)


  Results = get_significantSNPs(PriorWithBF$SNPs, sign_method, sign_thresh, prune_res, save_files, verbose)
  log_info = c(log_info, Results$log_info)



  ### go back to inital folder ###
  log_info = c(log_info, "", "")
  tmp = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"
  log_info = update_log(log_info, tmp, verbose)

  if(save_files){
    setwd(InitPath)
  }
  rm(ZMatrix, envir = .GlobalEnv)

  if(TMP_FILE){
    system(paste0("rm ", GWAS))
    tmp = paste0("The temporary file \"", GWAS, "\" has been deleted \n")
    log_info = update_log(log_info, tmp, verbose)
  }

  ### write log_info File ###
  Time = as.integer((proc.time()-StartTime)[3])
  minutes <- as.integer(trunc(Time/60))
  seconds <- Time - minutes * 60

  tmp = paste0("Time of the analysis: ", minutes, " minute(s) and ", seconds, " second(s).  \n")
  log_info = update_log(log_info, tmp, verbose)

  log_info = apply(as.array(log_info), 1,function(x) gsub("\n", "", x, fixed=T))
  write(log_info, paste0(name,".log"))

  results=list()
  results$log_info_info = log_info
  results$significant_SNPs = Results$SNPs
  results$all_BFs = PriorWithBF$SNPs
  results$significant_studies = res_MR$coeffs
  results$all_MRcoeffs = Prior$all_coeffs

  class(results) = "bGWAS"
 return(results)
}
