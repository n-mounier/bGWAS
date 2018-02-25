###### Function to create full Z-matrix to compute prior ######



# #' Create full Z-matrix for prior computation
# #'
# #' From a list of significant studies, create the full Z-matrix that can be used to compute
# #' the prior.
# #'
# #' @inheritParams bGWAS
# #' @param studies The IDs of the significant studies selected by \code{identify_studiesMR()}
# #'        (numeric vector)
# #'
# #' @return An object containing Log file and pruned Z-Matrix of MR instrument (data table) + create a file if saveFiles=T
# #'




makeFull_ZMatrix <- function(studies=NULL, GWAS,  Z_matrices="~/Z_matrices", save_files=F, verbose=F) {
  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = update_log(Log, tmp, verbose)


  tmp = paste0("Selecting studies :\n")
  Log = update_log(Log, tmp, verbose)

  if(grepl("macOS", sessionInfo()$running)) {
    ZMatrix=data.table::fread(paste0("zcat < ",paste0(Z_matrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, studies+5 ), showProgress = FALSE)
  } else if(grepl("Linux", sessionInfo()$running)){
    ZMatrix=data.table::fread(paste0("zcat < ",paste0(Z_matrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, studies+5 ), showProgress = FALSE)
  } else {
    stop("Only UNIL and MAC OS are supported")
  }
  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = update_log(Log, tmp, verbose)

  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)



  # Add conventional GWAS column, at the end (make sure alleles are aligned)

  if(is.numeric(GWAS)){  # if GWAS from our data
    tmp = paste0("# Adding data from the conventional GWAS (ID=", GWAS, "): \n \"", list_files(IDs = GWAS) , "\" \n")
    Log = update_log(Log, tmp, verbose)


    if(grepl("macOS", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(Z_matrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5))
    } else if(grepl("Linux", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(Z_matrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5))
    } else {
      stop("Only UNIL and MAC OS are supported")
    }
    t
    # no need to check for alignment of alleles, just subset and rename the column
    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs,GWASData$rs),]
    ZMatrix[,  list_files(IDs = GWAS)  := GWASData[,6]]

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  } else if(is.character(GWAS)){  # if external GWAS
    tmp = paste0("# Adding data from the conventional GWAS : \n \"", GWAS, "\" \n")
    Log = update_log(Log, tmp, verbose)


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

    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs, unlist(GWASData[,..rs])),]
    # check alignment
    aligned = which(GWASData[,..alt, with=F] == ZMatrix$alt &
                      GWASData[,..ref, with=F] == ZMatrix$ref)
    swapped = which(GWASData[,..ref, with=F] == ZMatrix$alt &
                      GWASData[,..alt, with=F] == ZMatrix$ref)
   # weird = c(1:nrow(GWASData))[!c(1:nrow(GWASData)) %in% c(aligned, swapped)]

    GWASData[, myZ:= numeric()]
    GWASData[aligned, myZ:= GWASData[aligned, ..z, with=F]]
    GWASData[swapped, myZ:= -GWASData[swapped, ..z, with=F]]


    ZMatrix[, strsplit(GWAS, "/")[[1]][ length(strsplit(GWAS, "/")[[1]])]] = GWASData$myZ

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  } else if(is.data.frame(GWAS)){  # if GWAS is data.frame
    GName = attributes(GWAS)$GName
    tmp = paste0("# Adding data from the conventional GWAS : \n \"", GName,
                 "\" \n")
    Log = update_log(Log, tmp, verbose)


    rs = match(colnames(GWAS),c("snpid", "snp", "rnpid", "rs", "rsid"))
    rs = which(!is.na(rs))
    alt = match(colnames(GWAS),c("a1", "alts", "alt"))
    alt = which(!is.na(alt))
    ref = match(colnames(GWAS),c("a2", "ref", "a0"))
    ref = which(!is.na(ref))
    z = match(colnames(GWAS),c("z", "Z", "zscore"))
    z = which(!is.na(z))

    # keep the SNPs in our pruned matrix and order them correctly
    GWAS = GWAS[match(ZMatrix$rs, unlist(GWAS[,..rs])),]
    # check alignment
    aligned = which(GWAS[,..alt, with=F] == ZMatrix$alt &
                      GWAS[,..ref, with=F] == ZMatrix$ref)
    swapped = which(GWAS[,..ref, with=F] == ZMatrix$alt &
                      GWAS[,..alt, with=F] == ZMatrix$ref)
    #weird = c(1:nrow(GWAS))[!c(1:nrow(GWAS)) %in% c(aligned, swapped)]

    GWAS[, myZ:= numeric()]
    GWAS[aligned, myZ:= GWAS[aligned, ..z, with=F]]
    GWAS[swapped, myZ:= -GWAS[swapped, ..z, with=F]]

    ZMatrix[, GName] = GWAS$myZ

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }



  ZMatrix = ZMatrix[complete.cases(ZMatrix)]
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS")
  Log = update_log(Log, tmp, verbose)


  # write ZMatrix
  # only if write.file=T
  # if(save_files){
  #   write.table(ZMatrix, file= "Full_ZMatrix.csv", sep=",", row.names=F, quote=F)
  #   tmp = "The file \"Full_ZMatrix.csv\" has been successfully created \n"
  #   Log = c(Log, tmp)
  #   if(verbose) cat(tmp)
  # }

  res=list()
  res$log_info = Log
  res$mat = ZMatrix
  return(res)
}

