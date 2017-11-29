###### Functions for study selection for the prior ######



#' availableStudies
#'
#' Get list of studies with available GWAS summary statistics
#' @param verbose boolean, default = FALSE
#' @return Data.Frame of details for all available studies
#' @export

availableStudies <- function(verbose=F) {
  Studies = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"))
  return(Studies)
}


#' listFiles
#'
#' Get list of studies with available GWAS summary statistics
#' @return List of files
#' @export

listFiles <- function(verbose=F) {
  Files = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "File")
  return(Files)
}


#' listTraits
#'
#' Get list of traits with available GWAS summary statistics
#' @return List of traits
#' @export

listTraits <- function(verbose=F) {
  Traits = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "Trait")
  return(unique(Traits$Trait))
}


#' listConsortia
#'
#' Get list of consortia with available GWAS summary statistics
#' @return List of consortia
#' @export

listConsortia <- function(verbose=F) {
  Consortia = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), select = "Consortium")
  return(unique(Consortia))
}


#' selectStudies
#'
#' Allow the quick selection of a subset of studies for prior based on 3 criteria
#' @param ListOfStudies ID of the studies, by default, use all the studies in availableStudies()
#' @param File list of file names
#' @param Trait list of trait
#' @param Consortium vector, list of consortium
#' @param verbose boolean
#' @return List of studies that meet the criteria
#' @examples
#'   AllStudies = availableStudies()
#'   listTraits()
#'   MyStudies = selectStudies(Traits=c("Heart Rate", "Body mass index", "Smoking"))
#'   AllStudies[AllStudies$ID == MyStudies, ]
#' @export



selectStudies <- function(ListOfStudies=NULL, Files=NULL, Traits=NULL, Consortia=NULL, verbose=F) {
  # Check parameters
  if(is.null(c(Files,Traits,Consortia))) stop("You did not specify any criteria for the selection.")
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
    Studies <- Studies[Studies$File %in% Files,]
    Crit = c(Crit, c("files"))
    nS = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nS, " studies have been removed when selection using Files"))
  }
  # 2nd : selection based on Trait
  if(!is.null(Traits)){
    Studies <- Studies[Studies$Trait %in% Traits,]
    Crit = c(Crit, c("traits"))
    nT = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nT, " studies have been removed when selection using Traits"))
  }
  # 3rd : selection based on Consortium
  if(!is.null(Consortia)){
    Studies <- Studies[Studies$Consortium %in% Consortia,]
    Crit = c(Crit, c("consortia"))
    nC = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nC, " studies have been removed when selection using Consortia"))
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
  return(Studies$ID)
}
