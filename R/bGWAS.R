###### Main Function ######



#' bGWAS
#'
#' Performs a bayesian GWAS from Summary Statistics, using publicly available results
#' to calculate the prior effects of the SNPs and compare it to observed z-scores
#'
#'
#' @param Name The name of the analysis (character)
#' @param GWAS The path to the conventionnal GWAS of interest or ID of the GWAS from the
#'        list of studies availables (prior GWASs) (character or numeric)
#' @param ZMatrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @param PriorStudies The IDs of prior GWASs to use for the analysis, \code{default=NULL},
#'        will include all the prior GWASs available (numeric vector)
#' @param ListOfSNPs, The path to a file containing the rsids of the SNPs to use ,
#'        \code{default=NULL}, will use all the SNPs in common between prior GWASs and the
#'        conventionnal GWAS (character)
#'              # NOT IMPLEMENTED YET
#' @param MRthreshold The threshold used to select strong instruments for MR, should be lower
#'        than 1e-5, \code{default=1e-5} (numeric)
#' @param SignMethod The method used to identify significant SNPs, should be \code{"p"} for
#'        p-value or \code{"fdr"} for false discovery rate, \code{default="p"} (character)
#' @param SignThresh The threshold used to identify significant SNPs, \code{default="5.10e-8"}
#'        (numeric)
#' @param pruneRes A logical indicating if the results should be pruned (by distance, 500kb),
#'        \code{default=FALSE}
#' @param saveFiles A logical indicating if the results should be saved as files,
#'        \code{default=FALSE}
# #' @param OutPath character, path to the outputs, needed if saveFiles is TRUE, by default,
# #' current working dictory
#' @param verbose  A logical indicating if information on progress should be reported,
#'        \code{default=FALSE}
#' @details
#' \code{Name} and \code{GWAS} are required arguments.
#'
#' If \code{GWAS} is a path to a file (regular or .gz), this file should contain the following
#' columns : \cr
#' SNPID (rs numbers) should be : \code{rs}, \code{rsid}, \code{snp}, \code{snpid}, \code{rnpid} \cr
#' A1 should be : \code{a1}, \code{alt}, \code{alts} \cr
#' A2 should be : \code{a2}, \code{a0}, \code{ref} \cr
#' Z should be : \code{z}, \code{Z}, \code{zscore} \cr
#' If Z is not present, it can be calculated from BETA and SE, in this case, a temporary
#' gzipped file is created and removed after the analysis. \cr
#' BETA should be : \code{b}, \code{beta}, \code{beta1} \cr
#' SE should be : \code{se}, \code{std} \cr
#'
#' Z-Matrix files, containing Z-scores for all prior GWASs should be downloaded separately
#' and stored in \code{"~/ZMatrices"} or in the folder specified with the argument
#' \code{ZMatrices}. \cr
#' ## MORE INFO NEEDED HERE / HOW TO DOWNLOAD ID (see GitHub README)
#'
#' Use \code{\link{availableStudies}()} to see all the prior GWASs available.
#' Using one of them as your conventionnal GWAS (argument \code{GWAS} = numeric ID) will automatically
#' remove it from the list of prior GWASs used to build the prior.
#'
#' Use \code{\link{selectStudies}()} to automatically select the studies to
#' be included/excluded when building the prior (argument \code{PriorStudies}).
#'
#' Be careful, in the results, all the GWAS (conventional + prior) are aligned with UK10K
#' data for the analysis (some alleles might be swapped when comparing with the initial data)
#'
#' @return An object containing the significant SNPs + Files created if saveFiles=T
#' ## DESCRIBE THE FILES
#'
#'
#' @examples
#' # Permorm bGWAS, using a conventional GWAS from the list of prior GWASs
#' MyGWAS = 1
#' listFiles(MyGWAS)
#' \dontrun{
#' A = bGWAS(Name = "Test_UsingGWASfromList",
#'          GWAS = MyGWAS
#'          verbose=T)
#'          }
### TO BE DONE
#'# Permorm bGWAS, using a small conventional GWAS included in data and selecting a subset of
#'# studies for the prior
#' MyGWAS = system.file("Data/SmallGWAS_Pilling2017.csv", package="bGWAS")
#' MyStudies = selectStudies(includeTraits=c("Type 2 diabetes", "Smoking"),
#'                          includeFiles=c("jointGwasMc_HDL.txt.gz","jointGwasMc_LDL.txt.gz"))
#' \dontrun{
#' B = bGWAS(Name = "Test_UsingSmallGWAS",
#'          GWAS = MyGWAS,
#'          PriorStudies=MyStudies,
#'          verbose=T)
#'          }
#' @export



bGWAS <- function(Name,
                  GWAS,
                  ZMatrices = NULL,
                  PriorStudies = NULL,
                  ListOfSNPs = NULL,
                  MRthreshold = 1e-5,
                  SignMethod = "p",
                  SignThresh = 5*10e-8,
                  pruneRes = FALSE,
#                  OutPath = getwd(),
                  saveFiles = FALSE,
                  verbose = FALSE) {


  InitPath = getwd()
### TO BE DONE
  OutPath = getwd() # to be cleaned
  # depends if we allow the user to run the analysis in another directory
  # of if we force the use of getwd()

  StartTime =  proc.time()

  # used in the main function
  # automatically re-detected when needed by other sub-functions
  platform = c("Linux", "macOS", "W")[c(grepl("Linux", sessionInfo()$running)
                                        , grepl("macOS", sessionInfo()$running)
                                        , grepl("Windows", sessionInfo()$running))]
  if(platform=="W") stop("Windows is not supported yet")

  # initialization of log file
  Log = c()

  tmp = paste0("<<< Preparation of analysis >>> \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  ### check the parameters ###
  tmp = paste0("> Checking parameters \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  if(!is.logical(verbose)) stop("verbose : should be logical")


  ## Name of analysis
  if(!is.character(Name)) stop("Name : non-character argument") # should be a string

  tmp = paste0("The name of your analysis is: \"", Name, "\". \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  ### create the directory to store the results ###
  Dir = file.path(OutPath, Name)
  #  if the directory already exists : error
  if(file.exists(file.path(OutPath, paste0(Name, ".Log")))) stop("You already run an analysis with the same name in that directory,
                                                    please specify another name or choose another directory to run the analysis")
  ifelse(!dir.exists(Dir), dir.create(Dir), stop("You already run an analysis with the same name in that directory,
                                                 please specify another name or choose another directory to run the analysis"))



  ## ZMatrices
  if(is.null(ZMatrices)){
    if(!file.exists("~/ZMatrices/ZMatrix_Imputed.csv.gz")) stop("No ZMatrix_Imputed.csv.gz file in ~/ZMatrices")
    if(!file.exists("~/ZMatrices/ZMatrix_NotImputed.csv.gz")) stop("No ZMatrix_NotImputed.csv.gz file in ~/ZMatrices")
    ZMatrices = "~/ZMatrices"
  } else if (is.character(ZMatrices)){
    if(!file.exists(paste0(ZMatrices, "/ZMatrix_Imputed.csv.gz"))) stop("No ZMatrix_Imputed.csv.gz file in specified ZMatrices folder")
    if(!file.exists(paste0(ZMatrices, "/ZMatrix_NotImputed.csv.gz"))) stop("No ZMatrix_NotImputed.csv.gz file in specified ZMatrices folder")
    # Define if the path is relative or not
    if(!substr(ZMatrices,1,1) %in% c("/", "~")){
      ZMatrices=paste0(InitPath, "/", ZMatrices)
    }
  } else stop("ZMatrices : wrong format")

  tmp = paste0("The Z-Matrix files are stored in \"", ZMatrices, "\".  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)



  ## GWAS of interest, should be a path to a GWAS file (format ? .tar.gz or file ?) or an ID
  if(is.numeric(GWAS) && !GWAS %in% c(1:length(listFiles()))) { # if it is an ID
    stop("The ID specified as a conventional GWAS is not in the list")
    tmp = paste0("The conventional GWAS used as input is:",
                 listFiles(IDs=GWAS), " (ID = ",  GWAS,").  \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  } else if(is.character(GWAS)) { # if it is a file
    # First, does the file exists ?
    if(!file.exists(GWAS)) stop("GWAS : the file does not exist")
    # Then, check if it is an absolute or a relative path to the file
    if(!substr(GWAS,1,1) %in% c("/", "~")){ # We should work with absolute path to avoid errors
      GWAS=paste0(InitPath, "/", GWAS)
    }
    tmp = paste0("The conventional GWAS used as input is: \"",
                 strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

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
    TMP_FILE = F # flag : is a temporary file with Z-scores created ??
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
        print(TMP_Name)
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
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

    if(TMP_FILE){
      tmp = paste0("A temporary file with a Z column has been created : \"",
                   strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n")
      Log = c(Log, tmp)
      if(verbose) cat(tmp)
    }
  }

### TO BE DONE
  ## OutPath, check that the directory exist. Create it if necessary ?
  if(is.null(OutPath)) OutPath = getwd()
  if(!dir.exists(OutPath)) stop("OutPath : the directory does not exist")

  tmp = paste0("The analysis will be run in the folder: \"", OutPath, "\".  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  ## saveFiles
  if(!is.logical(saveFiles)) stop("saveFiles should be logical")
  if(saveFiles){
    tmp = paste0("Files will be saved in: \"", OutPath, "/", Name, "\".  \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  } else {
    tmp = paste0("Temporary files will be removed after the analysis.  \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }



### TO BE DONE
  ## ListOfSNPs
  # We should have at least XX SNPs / check Linux-MAC for Zcat
  # Also check the number or SNPs in common in the file if ListOfSNPs not specified ?
  if(!is.null(ListOfSNPs)){
    if(!is.character(ListOfSNPs)) stop("ListOfSNPs : non character")
    OurSNPsMR = data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=1, showProgress = FALSE)
    commonMR = table(ListOfSNPs %in% OurSNPsMR$rs)["TRUE"]
    if(commonMR<20) stop("ListOfSNPs : You should provide at least 20 SNPs that can be used as strong instruments for MR")
    OurSNPsAll = data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=1, showProgress = FALSE)
    commonAll = table(ListOfSNPs %in% OurSNPsAll$rs)["TRUE"]
    tmp = paste0(length(ListOfSNPs), " SNPs provided \n")
    tmp = c(tmp, paste0(commonMR, " can be used for MR \n"))
    tmp = c(tmp, paste0(commonAll, " can be used to compute prior \n"))
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }


  ## PriorStudies
  # check that all the files required exist in our list of studies
  # should be specified as "File ID"
  if(is.null(PriorStudies)) PriorStudies = c(1:length(listFiles()))
  if(!all(PriorStudies %in% c(1:length(listFiles())))) stop("PriorStudies : all the IDs provided should belong to the ones available")
  # if GWAS from data, make sure to remove it
  if(is.numeric(GWAS) && GWAS %in% PriorStudies){
    PriorStudies = PriorStudies[-GWAS]
    tmp = paste0("The study ", listFiles(ID=GWAS), " (ID=", GWAS, ") has been removed from the studies used to build the prior since it is used as conventionnal GWAS. \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
### TO BE DONE
    # check that the user did not use exactly the same study for GWAS and for Prior,
    # in this case the prior study is removed here...
  }


  ## MRthreshold -> should not be larger than 10-5, can only be more stringent
  if(!is.numeric(MRthreshold)) stop("MRthreshold : non-numeric argument")
  if(MRthreshold>10^5) stop("MRthreshold : superior to the threshold limit (10^-5)")

  tmp = paste0("The p-value threshold used for selecting MR instruments is: ", format(MRthreshold, scientific = T), ".  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  ## SignMethod -> should not be "p" or "fdr"
  if(!SignMethod %in% c("p", "fdr")) stop("SignMethod : method not accepted, should be p or fdr")

  ## SignThresh -> should be numeric and lower (or equal) to 1
  if(!is.numeric(SignThresh)) stop("SignThresh : non numeric threshold")
  if(SignThresh>1) stop("SignThresh : a threshold higher than 1 does not make sense")

  if(SignMethod=="p"){
    tmp = paste0("Significant SNPs will be identified according to p-value. The threshold used is :",
                 format(SignThresh, scientific = T), ".  \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  } else if(SignMethod=="fdr"){
    tmp = paste0("Significant SNPs will be identified according to FDR. The threshold used is :",
                 format(SignThresh, scientific = T), ".  \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  ## pruneRes
  if(!is.logical(pruneRes)) stop("pruneRes : should be logical")
  if(pruneRes){
    tmp = "Significant SNPs will be pruned (distance pruning, 500kb) \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }


  # Go into the analysis' directory
  setwd(Dir)


  ### 1 : create "Summarize_file" ###

  tmp = paste0("# Removing unused GWAS (if necessary?) and reading the summary information files of the studies used  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  Files_Info = availableStudies()
  Files_Info = Files_Info[Files_Info$ID %in% PriorStudies, ]

  write.table(Files_Info, file=paste0(Dir, "/Summarize_file.csv"), sep=",", quote=F, row.names=F )
  # -> make it nicer, reusable by the user
  # wait to write it and add info ?
  # one column name / one column trait / one column Ref / one column Cohort ?

  tmp = paste0("List of files : ", Dir, "/SummarizeFiles.csv has been successfully created.  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # 2 : Z-Matrix for MR
  Log = c(Log, "", "")
  tmp = paste0("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n",
               "<<< Identification of significant studies for MR >>>  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # We should keep the Z-Matrix creation outside of the study identification function
  # so that we can quickly re-run the second part using a file containing the Z-Matrix
  tmp = "> Creating the Z-Matrix of strong instruments \n"
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  # the "global z_matrix" for all GWAS should already be done, just select the studies kept for the prior + prune + add the GWAS of interest
  # makeMR_ZMatrix() create a ZMatrix file and returns the log
  MR_ZMatrix = makeMR_ZMatrix(PriorStudies, GWAS, MRthreshold, ZMatrices, saveFiles, verbose)
  Log = c(Log, MR_ZMatrix$Log)

  Log=c(Log,"")
  tmp = paste0("> Performing MR  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  MR_Res = identify_StudiesMR(ZMatrix=MR_ZMatrix$Mat, saveFiles, verbose)
  Log = c(Log, MR_Res$Log)

  # 4 : Compute Prior
  Log = c(Log, "", "")
  tmp = paste0("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n",
                       "<<< Estimation of the prior >>>  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  tmp = paste0("> Creating the full Z-Matrix  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  Studies = selectStudies(includeFiles=MR_Res$Studies$study_selected)
  Full_ZMatrix = makeFull_ZMatrix(Studies, GWAS, ZMatrices, saveFiles, verbose)
  Log = c(Log, Full_ZMatrix$Log)

  tmp = paste0("> Computing prior  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  Prior = compute_Prior(MR_Res$Studies, MR_ZMatrix$Mat, Full_ZMatrix$Mat, saveFiles, verbose)
  Log = c(Log, Prior$Log)

  ##### COMPUTE THE BAYES FACTOR AND THE P-VALUE #####
  Log = c(Log, "", "")
  tmp = paste0("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n",
               "<<< Calculation of Bayes Factors and p-values >>>  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # Compute BFs for 100 nulls + for our GWAS of Interest
  # Calculate the p-values from these BFs (comparison with the nulls)
  # This script create a file containing all SNPs in common between prior file / imputed files
  tmp = paste0("> Calculating them for all SNPs  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  PriorWithBF = request_BFandP(Prior$Prior, saveFiles, verbose)
  Log = c(Log, PriorWithBF$Log)



  ##### IDENTIFY SIGNIFICANT SNPS AT 5% FDR + PRUNING #####
#  r2PruningPost=0.0
#  FDRthreshold=0.05
#  system(paste(".//Scripts/GetResults", GWASofInterest, Name, r2PruningPost, PriorName, FDRthreshold))
  tmp = paste0("> Pruning and identifying significant SNPs \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

#  print(paste0("Results File : Results/",Name, "/",
#               GWASofInterest, "-r2post-",
#               r2PruningPost, "-fdr-", FDRthreshold, ".csv has been successfully created."))

  # 5 : Bayes Factors and p-values
  # save it in ouptut/Name/...

  Results = get_SignificantSNPs(PriorWithBF$SNPs, SignMethod, SignThresh, pruneRes, saveFiles, verbose)
  Log = c(Log, Results$Log)



  ### go back to inital folder ###
  Log = c(Log, "", "")
  tmp = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  setwd(InitPath)
  rm(ZMatrix, envir = .GlobalEnv)
  # if saveFiles=F, remove all the created files (normally no file created) but temporary folder !
  if(!saveFiles)  system(paste0("rm -rf ", Name))
  if(TMP_FILE){
    system(paste0("rm ", GWAS))
    tmp = paste0("The temporary file \"", GWAS, "\" has been deleted \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  ### write Log File ###
  Time = as.integer((proc.time()-StartTime)[3])
  minutes <- as.integer(trunc(Time/60))
  seconds <- Time - minutes * 60

  tmp = paste0("Time of the analysis: ", minutes, " minute(s) and ", seconds, " second(s).  \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  Log = apply(as.array(Log), 1,function(x) gsub("\n", "", x, fixed=T))
  write(Log, paste0(Name,".log"))

  return(Results$SNPs)

}
