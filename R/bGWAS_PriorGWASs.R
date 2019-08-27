###### Functions for study selection for the prior ######



#' List prior GWASs
#'
#' Lists the studies that can be used as prior GWASs
#' @param IDs the IDs of the studies to print, \code{default="~/ZMatrices/"} will list
#'        all of them (numeric),
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @return data.frame containing prior GWASs information
#' @export

list_priorGWASs <- function(IDs=NULL, Z_matrices = "~/ZMatrices/") {
  if(!is.null(IDs) && !is.numeric(IDs)) stop("ID : should be numeric")
  Studies = data.table::fread(file.path(Z_matrices, "AvailableStudies.tsv"), showProgress = FALSE, data.table = F)
  if(!"Name" %in% colnames(Studies)) stop("Z_matrices : please use the most recent version")
  
  if(!is.null(IDs)){
    Studies %>%
      slice(IDs) -> Studies
  }
  return(Studies)
}


#' List prior GWASs files
#' 
#' Lists the filenames of the prior GWASs
#' @param IDs the IDs of the studies to print, \code{default="~/ZMatrices/"} will list
#'        all of them (numeric),
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @return List of files
#' @export

list_files <- function(IDs=NULL, Z_matrices = "~/ZMatrices/") {
  if(is.null(IDs)) {
    Files = data.table::fread(file.path(Z_matrices, "AvailableStudies.tsv"), select = "File", showProgress = FALSE, data.table=F)
  } else {
    # check that the IDs exists
    Studies = list_priorGWASs(Z_matrices=Z_matrices)
    if(!all(IDs %in% Studies$ID)) print("Please check the IDs, some of them do not match")
    Files = subset(Studies, .data$ID %in% IDs)
  }
  return(Files$File)
}


#' List prior GWASs traits
# 
#' Lists the traits of the prior GWASs
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @return List of traits
#' @export

list_traits <- function(Z_matrices = "~/ZMatrices/") {
  Traits = data.table::fread(file.path(Z_matrices, "AvailableStudies.tsv"), select = "Trait", showProgress = FALSE, data.table = F)$Trait
  return(unique(Traits))
}



# #' Get names from file for a set of Prior GWASs
# #'
# #' Get names from file for a set of Prior GWASs (useful for more user friendly reports, plots...)
# #' @param Files 
# #'        (character)
# #' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
# #'        (character)
# #' @return Names
# NOT EXPORTED


get_names <- function(Files, Z_matrices = "~/ZMatrices/") {
  Studies = data.table::fread(file.path(Z_matrices, "AvailableStudies.tsv"), showProgress = FALSE, data.table = F)
  if(any(!Files %in% Studies$File)) stop("Some files are not part of Prior GWASs")
  Studies %>%
    slice(match(Files, .data$File)) %>%
    pull(.data$Name)  -> Names
  return(Names)
}



#' Select prior GWASs
#'
#' Allow the quick selection of a subset of prior GWASs based on 2 criteria. First, include all the files
#' specified (if all including parameters are NULL, include all studies), and then remove all the files specified
#' (if all excluding parameters are NULL, keep all studies included at the step before)
#' @param include_files list of file names (see \code{\link{list_files}()}) (character)
#' @param include_traits list of trait (see \code{\link{list_traits}()}) (character)
#' @param exclude_files list of file names (see \code{\link{list_files}()}) (character)
#' @param exclude_traits list of trait (see \code{\link{list_traits}()}) (character)
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @param verbose boolean, default = FALSE
#' @return IDs (numeric) of studies that meet the criteria
#' @examples
#'   \dontrun{
#'   AllStudies = list_priorGWASs()
#'   list_traits()
#'   MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body Mass Index", "Smoking"))
#'   AllStudies[AllStudies$ID %in% MyStudies, c("ID", "Name", "Trait", "File")]}
#' @export


select_priorGWASs <- function(include_files=NULL, include_traits=NULL,
                          exclude_files=NULL, exclude_traits=NULL, 
                          Z_matrices = "~/ZMatrices/", verbose=F) {
  # Check parameters
  if(is.null(c(include_files,include_traits, exclude_files,exclude_traits))) stop("You did not specify any criteria for the selection.")
  Studies = list_priorGWASs(Z_matrices = Z_matrices)

  if(!all(include_files %in% list_files(Z_matrices = Z_matrices))) stop("Some files specified in Files are not correct")
  if(!all(exclude_files %in% list_files(Z_matrices = Z_matrices))) stop("Some files specified in Files are not correct")
  if(!all(include_traits %in% list_traits(Z_matrices))) stop("Some traits in Traits specified are not correct")
  if(!all(exclude_traits %in% list_traits(Z_matrices))) stop("Some traits in Traits specified are not correct")
 
  # Should not be possible to include / exclude the same file given different criteria
  # or it could...
#  allExcluded = Studies[File %in% exclude_files | Trait %in% exclude_traits | Consortium %in% excludeConsortia,"ID"]
#  allIncluded = Studies[File %in% include_files | Trait %in% include_traits | Consortium %in% includeConsortia,"ID"]
#  if(length(intersect(allExcluded$ID, allIncluded$ID))!=0) stop(paste0("You can't include and exclude the same file ! Problematic ID(s) = ", intersect(allExcluded$ID, allIncluded$ID)))

  n = nrow(Studies)
  Crit = c()

  # if inclusion criteria
  if(!is.null(include_traits) | !is.null(include_files)){
    AllStudies = Studies
    Studies = Studies[0,] # create empty data table
    n = nrow(Studies)

    ## A : do all inclusion
    # 1st : selection based on Study Name
    if(!is.null(include_files)){
      Studies <- rbind(Studies, AllStudies[AllStudies$File %in% include_files,])
      Crit = c(Crit, c("inclusion of files"))
      nS = nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nS, " studies have been included when selection using include_files"))
    }
    # 2nd : selection based on Trait
    if(!is.null(include_traits)){
      Studies <- rbind(Studies, AllStudies[AllStudies$Trait %in% include_traits,])
      Crit = c(Crit, c("inclusions of traits"))
      nT =  nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nT, " studies have been included when selection using include_traits"))
    }
  }

  ## B : do all exclusion
  # 1st : selection based on Study Name
  if(!is.null(exclude_files)){
    Studies = Studies[!(Studies$File %in% exclude_files),]
    Crit = c(Crit, c("exclusion of files"))
    nS = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nS, " studies have been removed when selection using exclude_files"))
  }
  # 2nd : selection based on Trait
  if(!is.null(exclude_traits)){
    Studies = Studies[!(Studies$Trait %in% exclude_traits),]
    Crit = c(Crit, c("exclusion of traits"))
    nT = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nT, " studies have been removed when selection using exclude_traits"))
  }
 
  if(verbose) print(paste0(n, " studies left after selection on ", paste(Crit, collapse=", ")))
  return(Studies$ID)
}


