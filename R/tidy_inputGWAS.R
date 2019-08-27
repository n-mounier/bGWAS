###### Function to create a nice data.frame for any type of input GWAS######



# #' Tidy input GWAS
# #'
# #' From the GWAS argument of the main bGWAS function, create a nice/tidy data.frame
# #' that can be used by all other functions.
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED



tidy_inputGWAS <- function(GWAS, Z_matrices = "~/ZMatrices/", verbose=FALSE){
  
  Log = c()
  
  
  tmp = paste0("# Preparation of the data... \n")
  Log = update_log(Log, tmp, verbose)
  
  # returned and used by other functions later to know if rescaled effects should be calculated
  rescaling <- FALSE

  # only to know how to handle data (calculate Z from beta/se or not)
  getZ <- FALSE

  
  GWASnames = list(SNPID = c("rsid", "snpid", "snp", "rnpid", "rs"),
                   ALT = c("a1", "alt", "alts"),
                   REF = c("ref", "a0", "a2"),
                   BETA = c("beta", "b", "beta1"),
                   SE = c("se", "std"),
                   Z = c("z", "zscore"))
  
  ## GWAS of interest, should be a path to a GWAS file (format ? .tar.gz or file ?), a data.frame or an ID
  
  if(is.numeric(GWAS)) { # if it is an ID
    
    rescaling = FALSE # no rescaling, no AF in ZMatrices files
    
    if(!GWAS %in% c(1:length(list_files(Z_matrices = Z_matrices)))) stop("The ID specified as a conventional GWAS is not in the list", call. = FALSE)
    
    GWASData <- data.table::fread(file.path(Z_matrices, "ZMatrix_Imputed.csv.gz"), select=c(1:5, GWAS+5), showProgress = FALSE, data.table=F)
    attributes(GWASData)$GName <- colnames(GWASData)[6]
    colnames(GWASData)[6] = "Z"
    
    tmp = paste0("The conventional GWAS used as input is:",
                 list_files(IDs=GWAS, Z_matrices = Z_matrices), " (ID = ",  GWAS,").  \n")
    Log = update_log(Log, tmp, verbose)
    
    colNames = c("rsid", "alt", "ref", "z_obs")

    GWASData %>%       
      select(1,4,5,6) %>%
      stats::setNames(colNames) -> GWASData_clean
    
  } else if(is.character(GWAS) | is.data.frame(GWAS)) { # if it is a file or a data.frame
    
    # file
    if(is.character(GWAS)) {
      
      # First, does the file exists ?
      if(!file.exists(GWAS)) stop("GWAS : the file does not exist", call. = FALSE)
      tmp = paste0("The conventional GWAS used as input is: \"",
                   strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n")
      Log = update_log(Log, tmp, verbose)
      
      # Check colnames...
      HeaderGWAS = colnames(data.table::fread(GWAS, nrows = 1, showProgress = FALSE, data.table=F))

      
    } else if(is.data.frame(GWAS)){ # if data.frame
      # add attribute GName to the data.frame, to be re-used in other subfunctions
      attributes(GWAS)$GName =  deparse(substitute(GWAS)) # get the "name" of the object used as an argument in the function
      # we want a data.frame (tidyr), not a data.table
      if(data.table::is.data.table(GWAS)) GWAS = as.data.frame(GWAS)
      tmp = paste0("The conventional GWAS used as input the object: \"",
                   attributes(GWAS)$GName, "\".  \n")
      Log = update_log(Log, tmp, verbose)
      HeaderGWAS = colnames(GWAS)
      
    }
    
    
    
    
    HeaderGWAS = tolower(HeaderGWAS)

    
    
    if(all(!HeaderGWAS %in% GWASnames[["SNPID"]])) stop("GWAS : no SNPID column", call. = FALSE)
    # how to deal with multiple rsid / snpid columns ???
    # here, we don't care, we need at least one
    tmp = paste0("   SNPID column, ok")
    
    if(all(!HeaderGWAS %in% GWASnames[["ALT"]])) stop("GWAS : no ALT column", call. = FALSE)
    tmp = c(tmp, paste0("ALT column, ok"))
    
    if(all(!HeaderGWAS %in% GWASnames[["REF"]])) stop("GWAS : no REF column", call. = FALSE)
    tmp = c(tmp, paste0("REF column, ok"))
    
    
    # if beta + se
    if(!all(!HeaderGWAS %in% GWASnames[["BETA"]]) & !all(!HeaderGWAS %in% GWASnames[["SE"]])){
      tmp = c(tmp, paste0("BETA column, ok"))
      tmp = c(tmp, paste0("SE column, ok"))
      rescaling = TRUE
      tmp2 = paste0("\nPosterior effects will be rescaled using BETA and SE.")
      if(all(!HeaderGWAS %in% GWASnames[["Z"]])) getZ=TRUE
      # if z, ok too
    } else if(!all(!HeaderGWAS %in% GWASnames[["Z"]])){
      tmp = c(tmp, paste0("Z column, ok"))
    } else {
      stop("GWAS : no effect (BETA/SE or Z) column(s)", call. = FALSE)
    }
    
    
    
    
    tmp = paste(tmp, collapse= " - ")
    Log = update_log(Log, tmp, verbose)
    if(exists("tmp2"))  Log = update_log(Log, tmp2, verbose)

    
    
    
    # if headers ok, and file or data.frame get GWASData
    if(is.character(GWAS)) {
      # Get the full data
      GWASData = data.table::fread(GWAS, showProgress = FALSE, data.table=F)
      attributes(GWASData)$GName = basename(GWAS)
      
    } else if(is.data.frame(GWAS)){ # if data.frame
      # add attribute GName to the data.frame, to be re-used in other subfunctions
      GWASData=GWAS
      rm(GWAS)
    }
    
    # by default, always use the first column of input GWAS matching the names
    # use col numbers because of different lower/upper possibilities
    SNPID = match(HeaderGWAS, GWASnames[["SNPID"]])
    SNPID = which(!is.na(SNPID))[1]
    ALT = match(HeaderGWAS, GWASnames[["ALT"]])
    ALT = which(!is.na(ALT))[1]
    REF = match(HeaderGWAS, GWASnames[["REF"]])
    REF = which(!is.na(REF))[1]
    BETA = match(HeaderGWAS, GWASnames[["BETA"]])
    BETA = which(!is.na(BETA))[1]
    SE = match(HeaderGWAS, GWASnames[["SE"]])
    SE = which(!is.na(SE))[1]
    ZSTAT = match(HeaderGWAS, GWASnames[["Z"]])
    ZSTAT = which(!is.na(ZSTAT))[1]
    
    
    colNumbers = c(SNPID, ALT, REF, BETA, SE, ZSTAT)
    colNames = c("rsid", "alt", "ref", "beta", "se", "z_obs")
    colNames = colNames[!is.na(colNumbers)]
    colNumbers = colNumbers[!is.na(colNumbers)]
    
    
    GWASData %>%
      select(colNumbers) %>%
      stats::setNames(colNames) -> GWASData_clean
    
    #if no Z, but BETA and SE, calculate Z
    if(getZ){
      GWASData_clean %>%
        mutate(z_obs = .data$beta/.data$se) -> GWASData_clean
    }
    
    
  } else stop("GWAS : unrecognized format", call. = FALSE)
  
  res=list(log_info = Log,
           GWAS = GWASData_clean,
           GWAS_Name = attributes(GWASData)$GName,
           rescaling = rescaling)
  return(res)
}
