###### Main Function - Using pre-computed prior ######



#' bGWAS - Main Function, using pre-computed prior
#'
#' Performs a bayesian GWAS
#' @param Name character, name of your analysis
#' @param GWAS character, path to your conventionnal GWAS of interest - columns needed :
##' @param GWAS character or numeric, path to your conventionnal GWAS of interest or ID of the GWAS
##' @param MRthreshold numeric, threshold used to select strong instruments for MR, should be lower than 1e-5, by default=1e-5
##' @param IncludeForPrior vector, list of files to include when creating the prior, by default=NULL, will include all the studies available.
#'                        See XXX for the full list.
#' @param saveFiles boolean,
#' @param OutPath character, path to the outputs, needed if saveFiles is TRUE, by default, current working dictory
#' @param verbose boolean,
#' @return Value Res
#' @export
#'
#'
#'


bGWASfromPrior <- function(Name,
                  GWAS,
                  IncludeForPrior=NULL,
                  ListOfSNPs=NULL,
                  MRthreshold=10e-5,
                  OutPath=getwd(),
                  saveFiles=T,
                  verbose=F) {


  InitPath = getwd()
  StartTime =  proc.time()
  platform = c("Linux", "macOS", "W")[c(grepl("Linux", sessionInfo()$running)
                                        , grepl("macOS", sessionInfo()$running)
                                        , grepl("Windows", sessionInfo()$running))]
  # can be useful in the main function ?
  # automatically detected when needed by other sub-functions

  Log = c()

  tmp = paste0("### Preparation of analysis ###")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  ### check the parameters ###
  tmp = paste0("# Checking parameters")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  ## Name of analysis
  if(!is.character(Name)) stop("Name : non-character argument") # should be a string

  tmp = paste0("The name of your analysis is: ", Name, ".")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  ## GWAS of interest, should be a path to a GWAS file (format ? .tar.gz or file ?)
  # First, does the file exists ?
  if(!file.exists(GWAS)) stop("GWAS : the file does not exist")
  # Then, check if it is an absolute or a relative path to the file
  if(!grepl(InitPath, GWAS)){ # We should work with absolute path to avoid errors
    GWAS=paste0(InitPath, "/", GWAS)
  }

  tmp = paste0("The conventional GWAS used as input is: \"",
               strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # Check colnames...
  if(!grepl(".gz", GWAS)){
    # ...if regular file
    HeaderGWAS = colnames(data.table::fread(GWAS, nrows = 0))
  } else if(grepl(".gz", GWAS)) {
    # ...if tar.gz
    HeaderGWAS = colnames(data.table::fread(paste0("zcat < ", GWAS), nrows = 0))
  }

  if(all(!HeaderGWAS %in% c("rsid", "snpid", "snp", "rnpid"))) stop("GWAS : no SNPID column")
  # how to deal with multiple rsid / snpid columns ???
  # here, we don't care, we need at least one
  tmp = paste0("SNPID column, ok")
  if(all(!HeaderGWAS %in% c("a1", "alts"))) stop("GWAS : no ALT column")
  tmp = c(tmp, paste0("ALT column, ok"))
  if(all(!HeaderGWAS %in% c("a2", "a0", "ref"))) stop("GWAS : no REF column")
  tmp = c(tmp, paste0("REF column, ok "))
  TMP_FILE = F # flag : is a temporary file with Z-scores created ??
  if(all(!HeaderGWAS %in% c("z", "Z"))){
    # allow for beta + se to calculate Z ???
    if(!all(!HeaderGWAS %in% c("b", "beta", "beta1")) & !all(!HeaderGWAS %in% c("se", "std"))){
      # if beta + se : read the data
      if(!grepl(".gz", GWAS)){
        # ...if regular file
        DataGWAS = data.table::fread(GWAS)
      } else if(grepl(".gz", GWAS)) {
        # ...if tar.gz
        DataGWAS = data.table::fread(paste0("zcat < ", GWAS))
      }
      # calculate Z
      DataGWAS$Z = DataGWAS[,HeaderGWAS[HeaderGWAS %in% c("b", "beta", "beta1")], with=F] /
        DataGWAS[,HeaderGWAS[HeaderGWAS %in% c("se", "std")], with=F]
      # keep only relevant column to save space
      DataGWAS = DataGWAS[,c(HeaderGWAS[HeaderGWAS %in% c("rsid", "snpid", "snp", "rnpid")][1], # SNPID, order by relevance :
                             # if there is an "rsid" and a "snpid" -> choose "rsid" !
                             HeaderGWAS[HeaderGWAS %in% c("a1", "alts")], # ALT
                             HeaderGWAS[HeaderGWAS %in% c("a2", "a0", "ref")], # REF
                             "Z"), # Z
                          with=F]
      # write the data (as tar.gz) and change GWAS name to the created file
      TMP_Name = paste0(gsub(".gz", "", GWAS), "_withZ.gz")
      write.table(DataGWAS, file=gzfile(TMP_Name), sep="\t",
                  quote=F, row.names=F)
      GWAS = TMP_Name
      # flag the created file
      TMP_FILE = T
      tmp = c(tmp, paste0("Z column, created"))
    } else {
      stop("GWAS : no Z-SCORE column")
    }
  } else {
    tmp = c(tmp, paste0("Z column, ok "))
  }
  tmp = paste(tmp, collapse= "-- ")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  if(TMP_FILE){
    tmp = paste0("A temporary file with a Z column has been created : ",
                 strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], ".")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  ## OutPath, check that the directory exist. Create it if necessary ?
  if(is.null(OutPath)) OutPath = getwd()
  if(!dir.exists(OutPath)) stop("OutPath : the directory does not exist")

  tmp = paste0("The analysis will be run in the folder: ", OutPath, ".")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  if(!is.logical(saveFiles)) stop("saveFiles should be logical")
  if(saveFiles){
    tmp = paste0("Files will be saved in: ", OutPath, "/", Name, ".")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  } else {
    tmp = paste0("Temporary files will be removed after the analysis.")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  ## MRthreshold -> should not be larger than 10-5, can only be more stringent
  if(!is.numeric(MRthreshold)) stop("MRthreshold : non-numeric argument")
  if(MRthreshold>10^5) stop("MRthreshold : superior to the threshold limit (10^-5)")

  tmp = paste0("The p-value threshold used for selecting MR instruments is: ", MRthreshold, ".")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  ## IncludeForPrior
  # check that all the files required exist in our list of studies
  # should be specified as "File names"
  if(is.null(IncludeForPrior)) IncludeForPrior = listFiles()
  if(!all(IncludeForPrior %in% listFiles())) stop("IncludeForPrior : all the files provided should be in the ones available")

  ## ListOfSNPs
  # check that all the files required exist in our list of studies
  # should be specified as "File names"
  if(is.null(IncludeForPrior)) IncludeForPrior = listFiles()
  if(!all(IncludeForPrior %in% listFiles())) stop("IncludeForPrior : all the files provided should be in the ones available")



  # Scheme ???
  # Not needed
  # -> always use :
  # - StepUpAIC (and not StepUpBIC / StepDnAIC / StepDnBIC)
  # - INpXXtozero / OUTpXXtozero ?
  # - RmDupeStudies ? (really needed when using AIC ?)
  # - pAIC-pdt (t : redundant with RmDupeStudies???)
  # - PrioVar ?
  # - dropChinMDD ? (just remove it from our list !!!)
  # - keepfit0 ?


  ### create the directory to store the results ###
  Dir = file.path(OutPath, Name)
  #  if the directory already exists : error
  ifelse(!dir.exists(Dir), dir.create(Dir), stop("You already run an analysis with the same name in that directory,
                                                 please specify another name or choose another directory to run the analysis"))
  setwd(Dir)





  ### 1 : create "Summarize_file" ###

  tmp = paste0("# Removing unused GWAS (if necessary?) and reading the summary information files of the studies used")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  Files_Info = availableStudies()
  Files_Info = Files_Info[Files_Info$File %in% IncludeForPrior, ]

  write.table(Files_Info, file=paste0(Dir, "/Summarize_file.csv"), sep=",", quote=F, row.names=F )
  # -> make it nicer, reusable by the user
  # wait to write it and add info ?
  # one column name / one column trait / one column Ref / one column Cohort ?

  tmp = paste0("List of files : ", Dir, "/SummarizeFiles.csv has been successfully created.")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # 2 : Z-Matrix for MR
  Log = c(Log, "", "")
  tmp = paste0("### Identification of significant studies for MR ###")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # We should keep the Z-Matrix creation outside of the study identification function
  # so that we can quickly re-run the second part using a file containing the Z-Matrix
  tmp = paste0("# Creating the Z-Matrix")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  # the "global z_matrix" for all GWAS should already be done, just select the studies kept for the prior + prune + add the GWAS of interest
  # save it in ouptut/Name/...
  # makeMR_ZMatrix() create a ZMatrix file and returns the log
  MR_ZMatrix = makeMR_ZMatrix(IncludeForPrior, GWAS, MRthreshold, verbose)
  Log = c(Log, MR_ZMatrix$Log)

  tmp = paste0("# Performing MR")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  MR_Res = identify_StudiesMR(ZMatrix=MR_ZMatrix$Mat, verbose)
  Log = c(Log, MR_Res$Log)

  # 4 : Compute Prior
  # save it in ouptut/Name/...

  # At this point, we need to create a Z matrix of all SNPs of interest
  #  ListOfSNPs = "whichSNPs"
  #  ListOfSelectedStudies = paste0(ZMatrix, "-scheme-", scheme)
  #  system(paste(".//Scripts/MakeFullMatrix", Name, "Useless", GWASofInterest, ListOfSelectedStudies, ListOfSNPs))

  #  print(paste0("Z Matrix for all SNPs : ZMatrix/ZMatrix_", ListOfSNPs, "_",
  #               Name, "-prune", as.integer(pruningDistance), "-r2", pruningR2, ".csv has been successfully created."))
  Log = c(Log, "", "")
  tmp = paste0("### Estimation of the prior ###")
  Log = c(Log, tmp)
  if(verbose) print(tmp)


  Full_ZMatrix = makeFull_ZMatrix(Dir, MR_Res$Studies, GWAS, verbose)

  # Build the prior
  # function that needs the ZMatrix and that's it ?
  #  system(paste(".//Scripts/ComputePrior", ZMatrix, scheme, Name, ListOfSNPs, GWASofInterest))


  #  print(paste0("Prior File : Prior/", ZMatrix, "-scheme-",
  #               scheme, ".csv has been successfully created."))
  #  print(paste0("Additional Files (with coefficients) have also been created in Prior/"))

  #  PriorName <- paste0( ZMatrix, "-scheme-",
  #                       scheme, ".csv")

  ##### COMPUTE THE BAYES FACTOR AND THE P-VALUE #####
  Log = c(Log, "", "")
  tmp = paste0("### BF and p ###")
  Log = c(Log, tmp)
  if(verbose) print(tmp)

  # Compute BFs for 100 nulls + for our GWAS of Interest
  # Calculate the p-values from these BFs (comparison with the nulls)
  # This script create a file containing all SNPs in common between prior file / imputed files:
  # rs logBF null_count BF_p

  #  system(paste(".//Scripts/RequestBFandPValue", PriorName, GWASofInterest, Name))


  #  print(paste0("BF_p File : CompareToNulls/",Name, "/",
  #               PriorName, ".csv has been successfully created."))


  ##### IDENTIFY SIGNIFICANT SNPS AT 5% FDR + PRUNING #####
  #  r2PruningPost=0.0
  #  FDRthreshold=0.05
  #  system(paste(".//Scripts/GetResults", GWASofInterest, Name, r2PruningPost, PriorName, FDRthreshold))


  #  print(paste0("Results File : Results/",Name, "/",
  #               GWASofInterest, "-r2post-",
  #               r2PruningPost, "-fdr-", FDRthreshold, ".csv has been successfully created."))

  # 5 : Bayes Factors and p-values
  # save it in ouptut/Name/...





  ### go back to inital folder ###
  setwd(InitPath)
  # if saveFiles=F, remove all the created files
  if(!saveFiles)  system(paste0("rm -rf ", Name))

  ### write Log File ###
  Time = as.integer((proc.time()-StartTime)[3])
  minutes <- as.integer(trunc(Time/60))
  seconds <- Time - minutes * 60

  tmp = paste0("Time of the analysis: ", minutes, " minute(s) and ", seconds, " second(s).")
  Log = c(Log, tmp)
  if(verbose) print(tmp)

  write(Log, paste0(Name,".log"))

}
