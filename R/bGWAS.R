###### Main Function ######



#' bGWAS - Main Function
#'
#' Performs a bayesian GWAS
#' @param Name character, name of your analysis
#' @param GWAS character, path to your conventionnal GWAS of interest
#' @param ExludeFromPrior vector, list of studies to exclude when creating the prior, see XXX, by default=NULL
#' @param out.path character, path to the outputs, by default, current working dictory
#' @param keepFiles ?
#' @param verbose boolean
#' @return Value Resultat
#' @export
#'
#'
#'


bGWAS <- function(Name, GWAS, threshold=NULL, scheme=NULL, ExludeFromPrior=NULL, verbose=F) {
  # check the parameters / unit test



  # create the directory to store the results / if the directory already exists : error "please specify a new name or delete old analysis results"

  #stopif()
  # 1 : create "Summarize_file" -> make it nicer, reusable by the user

  if(verbose) print(paste0("List of files : SummarizeFiles_", Name, ".csv has been successfully created."))

  # 2 : Z-Matrix
  # the "global z_matrix" for all GWAS should already be done, just select the studies kept for the prior + prune + add the GWAS of interest
  # save it in ouptut/Name/...

  # 3 : Identify significant studies

  # 4 : Compute Prior

  # 5 : Bayes Factors and p-values

  }
