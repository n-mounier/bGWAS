###### Functions for study selection for the prior ######



#' Get list of studies that can be used to build the prior
#'
#' Get list of studies with available GWAS summary statistics
#' @param verbose boolean, default = FALSE
#' @return Data.Frame of details for all available studies
#' @export

availableStudies <- function(verbose=F) {
  Studies = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), showProgress = FALSE)
  return(Studies)
}


#' Get file names of studies that can be used to build the prior
#'
#' Get list of studies with available GWAS summary statistics
#' @param IDs The IDs of the files to be listed, by default=NULL, list all files (numeric)
#' @return List of files
#' @export

listFiles <- function(IDs=NULL, verbose=F) {
  if(is.null(IDs)) {
    Files = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "File", showProgress = FALSE)$File
  } else {
    # check that the IDs exists
    Studies = availableStudies()
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

listTraits <- function(verbose=F) {
  Traits = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "Trait", showProgress = FALSE)$Trait
  return(unique(Traits))
}


#' Get consortia from studies that can be used to build the prior
#'
#' Get list of consortia with available GWAS summary statistics
#' @return List of consortia
#' @export

listConsortia <- function(verbose=F) {
  Consortia = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "Consortium", showProgress = FALSE)$Consortium
  return(unique(Consortia))
}


#' selectStudies
#'
#' Allow the quick selection of a subset of studies for prior based on 3 criteria. First, include all the files
#' specified (if all including parameters are NULL, include all studies), and then remove all the files specified
#' (if all excluding parameters are NULL, keep all studies included at the step before)
#' @param includeFiles list of file names
#' @param includeTraits list of trait
#' @param includeConsortium vector, list of consortium ### TO BE DONE
#' @param excludeFiles list of file names
#' @param excludeTraits list of trait
#' @param excludeConsortium vector, list of consortium ### TO BE DONE
#' @param verbose boolean
#' @return IDs of studies that meet the criteria
#' @examples
#'   AllStudies = availableStudies()
#'   listTraits()
#'   MyStudies = selectStudies(includeTraits=c("Heart Rate", "Body mass index", "Smoking"))
#'   #AllStudies[AllStudies$ID == MyStudies, ]
#' @export


selectStudies <- function(includeFiles=NULL, includeTraits=NULL, includeConsortia=NULL,
                          excludeFiles=NULL, excludeTraits=NULL, excludeConsortia=NULL, verbose=F) {
  # Check parameters
  if(is.null(c(includeFiles,includeTraits,includeConsortia, excludeFiles,excludeTraits,excludeConsortia))) stop("You did not specify any criteria for the selection.")
  Studies = availableStudies()

  if(!all(includeFiles %in% listFiles())) stop("Some files specified in Files are not correct")
  if(!all(excludeFiles %in% listFiles())) stop("Some files specified in Files are not correct")
  if(!all(includeTraits %in% listTraits())) stop("Some traits in Traits specified are not correct")
  if(!all(excludeTraits %in% listTraits())) stop("Some traits in Traits specified are not correct")
  if(!all(includeConsortia %in% listConsortia())) stop("Some consortia specified in Consortia are not correct")
  if(!all(excludeConsortia %in% listConsortia())) stop("Some consortia specified in Consortia are not correct")

  # Should not be possible to include / exclude the same file given different criteria
  # or it could...
#  allExcluded = Studies[File %in% excludeFiles | Trait %in% excludeTraits | Consortium %in% excludeConsortia,"ID"]
#  allIncluded = Studies[File %in% includeFiles | Trait %in% includeTraits | Consortium %in% includeConsortia,"ID"]
#  if(length(intersect(allExcluded$ID, allIncluded$ID))!=0) stop(paste0("You can't include and exclude the same file ! Problematic ID(s) = ", intersect(allExcluded$ID, allIncluded$ID)))

### TO BE DONE
  includeConsortia = NULL
  excludeConsortia = NULL

  n = nrow(Studies)
  Crit = c()

  # if inclusion criteria
  if(!is.null(includeConsortia) | !is.null(includeTraits) | !is.null(includeFiles)){
    AllStudies = Studies
    Studies = Studies[Trait=="",] # create empty data table
    n = nrow(Studies)

    ## A : do all inclusion
    # 1st : selection based on Study Name
    if(!is.null(includeFiles)){
      Studies <- rbind(Studies, AllStudies[AllStudies$File %in% includeFiles,])
      Crit = c(Crit, c("inclusion of files"))
      nS = nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nS, " studies have been included when selection using includeFiles"))
    }
    # 2nd : selection based on Trait
    if(!is.null(includeTraits)){
      Studies <- rbind(Studies, AllStudies[AllStudies$Trait %in% includeTraits,])
      Crit = c(Crit, c("inclusions of traits"))
      nT =  nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nT, " studies have been included when selection using includeTraits"))
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
  if(!is.null(excludeFiles)){
    Studies = Studies[!(File %in% excludeFiles),]
    Crit = c(Crit, c("exclusion of files"))
    nS = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nS, " studies have been removed when selection using excludeFiles"))
  }
  # 2nd : selection based on Trait
  if(!is.null(excludeTraits)){
    Studies = Studies[!(Trait %in% excludeTraits),]
    Crit = c(Crit, c("exclusion of traits"))
    nT = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nT, " studies have been removed when selection using excludeTraits"))
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





#' excludeStudies
#'
#' @param ListOfStudies ID of the studies, by default, use all the studies in availableStudies()
#' @param File list of file names
#' @param Trait list of trait
#' @param Consortium vector, list of consortium
#' @param verbose boolean
#' @return List of studies that meet the criteria
#' @examples
#'   AllStudies = availableStudies()
#'   listTraits()
#'   MyStudies = exludeStudies(Traits=c("Smoking", "Personality"))
#'   nrow(AllStudies)
#'   nrow(AllStudies[AllStudies$ID == MyStudies, ])
#' @export


excludeStudies <- function(ListOfStudies=NULL, Files=NULL, Traits=NULL, Consortia=NULL, verbose=F) {
  # Check parameters
  if(is.null(c(Files,Traits,Consortia))) stop("You did not specify any criteria for the exclusion.")
  Studies = availableStudies()
  if(!all(ListOfStudies %in% Studies$ID)) stop("Some IDs specified in ListOfStudies are not correct")
  if(!all(Files %in% listFiles())) stop("Some files specified in Files are not correct")
  if(!all(Traits %in% listTraits())) stop("Some traits in Traits specified are not correct")
  if(!all(Consortia %in% listConsortia())) stop("Some consortia specified in Consortia are not correct")
  if(!is.null(ListOfStudies)) Studies = Studies[Studies$ID %in% ListOfStudies,]
  n = nrow(Studies)
  Crit = c()
  # 1st : selection based on Study Name
  if(!is.null(Files)){
    Studies <- Studies[!Studies$File %in% Files,]
    Crit = c(Crit, c("files"))
    nS = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nS, " studies have been removed when selection using Files"))
  }
  # 2nd : selection based on Trait
  if(!is.null(Traits)){
    Studies <- Studies[!Studies$Trait %in% Traits,]
    Crit = c(Crit, c("traits"))
    nT = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nT, " studies have been removed when selection using Traits"))
  }
  # 3rd : selection based on Consortium
  if(!is.null(Consortia)){
    Studies <- Studies[!Studies$Consortium %in% Consortia,]
    Crit = c(Crit, c("consortia"))
    nC = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nC, " studies have been removed when selection using Consortia"))
  }
  if(verbose) print(paste0(n, " studies left after selection on ", paste(Crit, collapse=", ")))
  return(Studies$File)
}
