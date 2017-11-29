###### Main Function ######



#' bGWAS - Main Function
#'
#' Performs a bayesian GWAS
#' @param Name character, name of your analysis
#' @param GWAS character, path to your conventionnal GWAS of interest
#' @param ExludeFromPrior vector, list of studies to exclude when creating the prior, see XXX, by default=NULL
#' @param OutPath character, path to the outputs, by default, current working dictory
#' @param keepFiles ?
#' @param verbose boolean
#' @return Value Resultat
#' @export
#'
#'
#'


bGWAS <- function(Name, GWAS, OutPath=getwd(), threshold=NULL, scheme=NULL, ExludeFromPrior=NULL, verbose=F) {
  InitPath = getwd()
  Log = c()

  # check the parameters
  # Name of analysis
  if(!is.character(Name)) stop("Name : non-character argument") # should be a string
  Log = c(Log, paste0("The name of your analysis is: ", Name, ". \n"))
  # GWAS of interest, should be a path to a GWAS file (format ? .tar.gz or file ?), or a name of one of our GWAS, if so, should be automatically exluded from the prior
  Log = c(Log, paste0("The conventional GWAS used as input is: ", GWAS, ". \n"))
  # OutPath, check that the directory exist. Create it if necessary ?

  # create the directory to store the results / if the directory already exists : error "please specify a new name or delete old analysis results"
  Dir = file.path(OutPath, Name)
  ifelse(!dir.exists(Dir), dir.create(Dir), stop("You already run an analysis with the same name in that directory, please specify another name"))
  setwd(Dir)
  # Threshold

  # Scheme


  #stopif()
  # 1 : create "Summarize_file" -> make it nicer, reusable by the user
  # one column name / one column trait / one column Ref / one column Cohort ?

  if(verbose) print(paste0("List of files : SummarizeFiles_", Name, ".csv has been successfully created."))

  # 2 : Z-Matrix
  # the "global z_matrix" for all GWAS should already be done, just select the studies kept for the prior + prune + add the GWAS of interest
  # save it in ouptut/Name/...

  # 3 : Identify significant studies

  # 4 : Compute Prior

  # 5 : Bayes Factors and p-values

  # write Log File
  setwd(InitPath)
  }
