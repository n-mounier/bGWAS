###### Functions for study selection for the prior ######



#' Get list of studies that can be used to build the prior
#'
#' Get list of studies with available GWAS summary statistics
#' @param verbose boolean, default = FALSE
#' @return Data.Frame of details for all available studies
#' @export

list_priorGWASs <- function(verbose=F) {
  Studies = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), showProgress = FALSE)
  return(Studies)
}


#' Get file names of studies that can be used to build the prior
#'
#' Get list of studies with available GWAS summary statistics
#' @param IDs The IDs of the files to be listed, by default=NULL, list all files (numeric)
#' @return List of files
#' @export

list_files <- function(IDs=NULL, verbose=F) {
  if(is.null(IDs)) {
    Files = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "File", showProgress = FALSE)$File
  } else {
    # check that the IDs exists
    Studies = list_priorGWASs()
    if(!all(IDs %in% Studies$ID)) print("Please check the IDs, some of them do not match")
    Files = Studies$File[match(IDs, Studies$ID)]
  }
  return(Files)
}


#' Get traits from studies that can be used to build the prior
#'
#' Get list of traits with available GWAS summary statistics
#' @return List of traits
#' @export

list_traits <- function(verbose=F) {
  Traits = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "Trait", showProgress = FALSE)$Trait
  return(unique(Traits))
}


#' Get consortia from studies that can be used to build the prior
#'
#' Get list of consortia with available GWAS summary statistics
#' @return List of consortia
#' @export

list_consortia <- function(verbose=F) {
  Consortia = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "Consortium", showProgress = FALSE)$Consortium
  return(unique(Consortia))
}


#' select_priorGWASs
#'
#' Allow the quick selection of a subset of studies for prior based on 3 criteria. First, include all the files
#' specified (if all including parameters are NULL, include all studies), and then remove all the files specified
#' (if all excluding parameters are NULL, keep all studies included at the step before)
#' @param include_files list of file names
#' @param include_traits list of trait
#' @param include_consortium vector, list of consortium ### TO BE DONE
#' @param exclude_files list of file names
#' @param exclude_traits list of trait
#' @param exclude_consortium vector, list of consortium ### TO BE DONE
#' @param verbose boolean
#' @return IDs of studies that meet the criteria
#' @examples
#'   AllStudies = list_priorGWASs()
#'   list_traits()
#'   MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body mass index", "Smoking"))
#'   #AllStudies[AllStudies$ID %in% MyStudies, ]
#' @export


select_priorGWASs <- function(include_files=NULL, include_traits=NULL, includeConsortia=NULL,
                          exclude_files=NULL, exclude_traits=NULL, excludeConsortia=NULL, verbose=F) {
  # Check parameters
  if(is.null(c(include_files,include_traits,includeConsortia, exclude_files,exclude_traits,excludeConsortia))) stop("You did not specify any criteria for the selection.")
  Studies = list_priorGWASs()

  if(!all(include_files %in% list_files())) stop("Some files specified in Files are not correct")
  if(!all(exclude_files %in% list_files())) stop("Some files specified in Files are not correct")
  if(!all(include_traits %in% list_traits())) stop("Some traits in Traits specified are not correct")
  if(!all(exclude_traits %in% list_traits())) stop("Some traits in Traits specified are not correct")
  if(!all(includeConsortia %in% list_consortia())) stop("Some consortia specified in Consortia are not correct")
  if(!all(excludeConsortia %in% list_consortia())) stop("Some consortia specified in Consortia are not correct")

  # Should not be possible to include / exclude the same file given different criteria
  # or it could...
#  allExcluded = Studies[File %in% exclude_files | Trait %in% exclude_traits | Consortium %in% excludeConsortia,"ID"]
#  allIncluded = Studies[File %in% include_files | Trait %in% include_traits | Consortium %in% includeConsortia,"ID"]
#  if(length(intersect(allExcluded$ID, allIncluded$ID))!=0) stop(paste0("You can't include and exclude the same file ! Problematic ID(s) = ", intersect(allExcluded$ID, allIncluded$ID)))

### TO BE DONE
  includeConsortia = NULL
  excludeConsortia = NULL

  n = nrow(Studies)
  Crit = c()

  # if inclusion criteria
  if(!is.null(includeConsortia) | !is.null(include_traits) | !is.null(include_files)){
    AllStudies = Studies
    Studies = Studies[Trait=="",] # create empty data table
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
    # 3rd : selection based on Consortium
    if(!is.null(includeConsortia)){
      Studies <- rbind(Studies, AllStudies[AllStudies$Consortium %in% includeConsortia,])
      Crit = c(Crit, c("inclusion of consortia"))
      nC = nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nC, " studies have been included when selection using includeConsortia"))
    }
  }

  ## B : do all exclusion
  # 1st : selection based on Study Name
  if(!is.null(exclude_files)){
    Studies = Studies[!(File %in% exclude_files),]
    Crit = c(Crit, c("exclusion of files"))
    nS = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nS, " studies have been removed when selection using exclude_files"))
  }
  # 2nd : selection based on Trait
  if(!is.null(exclude_traits)){
    Studies = Studies[!(Trait %in% exclude_traits),]
    Crit = c(Crit, c("exclusion of traits"))
    nT = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nT, " studies have been removed when selection using exclude_traits"))
  }
  # 3rd : selection based on Consortium
  if(!is.null(excludeConsortia)){
    Studies = Studies[!(Consortium %in% excludeConsortia),]
    Crit = c(Crit, c("exclusion of consortia"))
    nC = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nC, " studies have been removed when selection using excludeConsortia"))
  }
  if(verbose) print(paste0(n, " studies left after selection on ", paste(Crit, collapse=", ")))
  return(Studies$ID)
}


