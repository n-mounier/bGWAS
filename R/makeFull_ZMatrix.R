###### Function to create full Z-matrix to compute prior ######



#' Create full Z-matrix for prior computation
#'
#' From a list of significant studies, create the full Z-matrix that can be used to compute
#' the prior.
#'
#' @inheritParams bGWAS
#' @param Studies The IDs of the significant studies selected by \code{identify_StudiesMR()}
#'        (numeric vector)
#'
#' @return An object containing Log file and pruned Z-Matrix of MR instrument (data table) + create a file if saveFiles=T
#'
#' @export




makeFull_ZMatrix <- function(Studies=NULL, GWAS,  ZMatrices="~/ZMatrices", saveFiles=F, verbose=F) {
  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  tmp = paste0("Selecting studies :\n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  if(grepl("macOS", sessionInfo()$running)) {
    ZMatrix=data.table::fread(paste0("zcat < ",paste0(ZMatrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, Studies+5 ), showProgress = FALSE)
  } else if(grepl("Linux", sessionInfo()$running)){
    ZMatrix=data.table::fread(paste0("zcat < ",paste0(ZMatrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, Studies+5 ), showProgress = FALSE)
  } else {
    stop("Only UNIL and MAC OS are supported")
  }
  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)


  # Add conventional GWAS column, at the end (make sure alleles are aligned)

  if(is.numeric(GWAS)){  # if GWAS from our data
    tmp = paste0("# Adding data from the conventional GWAS (ID=", GWAS, "): \n \"", listFiles(IDs = GWAS) , "\" \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

    if(grepl("macOS", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(ZMatrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5))
    } else if(grepl("Linux", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(ZMatrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5))
    } else {
      stop("Only UNIL and MAC OS are supported")
    }
    t
    # no need to check for alignment of alleles, just subset and rename the column
    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs,GWASData$rs),]
    ZMatrix[,  listFiles(IDs = GWAS)  := GWASData[,6]]

    tmp = "Done! \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

  } else if(is.character(GWAS)){  # if external GWAS
    tmp = paste0("# Adding data from the conventional GWAS : \n \"", GWAS, "\" \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

    if(!grepl(".gz", GWAS)){
      GWASData = data.table::fread(GWAS, showProgress = FALSE)
    } else if(grepl(".gz", GWAS)) {
      # if tar.gz
      GWASData = data.table::fread(paste0("zcat < ", GWAS), showProgress = FALSE)
    }

    rs = match(colnames(GWASData),c("snpid", "snp", "rnpid", "rs", "rsid"))
    rs = which(!is.na(rs))
    alt = match(colnames(GWASData),c("a1", "alts", "alt"))
    alt = which(!is.na(alt))
    ref = match(colnames(GWASData),c("a2", "ref", "a0"))
    ref = which(!is.na(ref))
    z = match(colnames(GWASData),c("z", "Z", "zscore"))
    z = which(!is.na(z))

    #zed = c("z", "Z")[z[!is.na(z)]]

    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs, unlist(GWASData[,..rs])),]
    # check alignment
    aligned = which(GWASData[,..alt, with=F] == ZMatrix$alt &
                      GWASData[,..ref, with=F] == ZMatrix$ref)
    swapped = which(GWASData[,..ref, with=F] == ZMatrix$alt &
                      GWASData[,..alt, with=F] == ZMatrix$ref)
    weird = c(1:nrow(GWASData))[!c(1:nrow(GWASData)) %in% c(aligned, swapped)]

    GWASData[aligned, myZ:= GWASData[aligned, ..z, with=F]]
    GWASData[swapped, myZ:= -GWASData[swapped, ..z, with=F]]
    GWASData[weird, myZ:= NA]

    ZMatrix[, strsplit(GWAS, "/")[[1]][ length(strsplit(GWAS, "/")[[1]])]] = GWASData$myZ

    tmp = "Done! \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }



  ZMatrix = ZMatrix[complete.cases(ZMatrix)]
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # write ZMatrix
  # only if write.file=T
  if(saveFiles){
    write.table(ZMatrix, file= "Full_ZMatrix.csv", sep=",", row.names=F, quote=F)
    tmp = "The file \"Full_ZMatrix.csv\" has been successfully created \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  res=list()
  res$Log = Log
  res$Mat = ZMatrix
  return(res)
}

