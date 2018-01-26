###### Function to create pruned Z-matrix of MR instruments ######



#' Create pruned Z-matrix of MR instruments
#'
#' From a list of studies, create pruned Z-matrix of MR instruments (significant at
#' a specified threshold)
#'
#' @inheritParams bGWAS
#'
#' @return Log file and pruned Z-Matrix of MR instrument + create a file if saveFiles=T
#' @export

makeMR_ZMatrix <- function(PriorStudies=NULL, GWAS, MRthreshold=10e-5, path="~/ZMatrices", saveFiles=F, verbose=F) {
  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  if(!is.null(PriorStudies)){
    tmp = paste0("Selecting studies :\n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    if(grepl("macOS", sessionInfo()$running)){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=c(1:5, PriorStudies+5))
    } else if(grepl("Linux", sessionInfo()$running)){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=c(1:5, PriorStudies+5))
    } else {
      stop("Only UNIL and MAC OS are supported")
    }

  } else {
    if(grepl("macOS", sessionInfo()$running)) {
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")))
    } else if(grepl("Linux", sessionInfo()$running)){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")))
    } else {
      stop("Only UNIL and MAC OS are supported")
    }
}

  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # select based on threshold if different from 1e- 5and removed rows without any Z-Score ok
  if(MRthreshold != 10e-5 ){
    tmp = paste0("# Thresholding... \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    Zlimit = qnorm(MRthreshold, lower.tail = F)
    SNPsToKeep = apply(ZMatrix[,-c(1:5)], 1, function(x) any(abs(x)>Zlimit))
    ZMatrix=ZMatrix[SNPsToKeep,]
    tmp = paste0(ncol(ZMatrix)-5, " studies left after thresholding \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs left after thresholding \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  } else {
    tmp = paste0("No thresholding needed... \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }

  # pruning
  # by default, distance pruning, but enable the LD pruning from Aaron's function to compare
  DistancePruning=T
  dist = 500000
  if(DistancePruning){
    tmp = paste0("Distance pruning... \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    tmp = paste0("   distance : ", dist/1000, "kb \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    ZMatrixPruned = prune_ZMatrix(ZMatrix, prune.dist = dist)
  } else {
    r2 = 0.8
    tmp = paste0("LD pruning... \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    tmp = paste0("   distance : ", dist/1000, "kb", " - r2 threshold : ", r2, "\n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
    ZMatrixPruned = prune_ZMatrix(ZMatrix, prune.dist = dist, r2.limit=r2)
  }

  tmp = paste0(format(nrow(ZMatrixPruned), big.mark = ",", scientific = F), " SNPs left after pruning \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # Add conventional GWAS column, at the end (make sure alleles are aligned)
  if(is.numeric(GWAS)){  # if GWAS from our data
    tmp = paste0("# Adding data from the conventional GWAS (ID=", GWAS, "): \n \"", listFiles(IDs = GWAS) , "\" \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

    if(grepl("macOS", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5))
    } else if(grepl("Linux", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5))
    } else {
      stop("Only UNIL and MAC OS are supported")
    }

    # no need to check for alignment of alleles, just subset and rename the column
    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrixPruned$rs,GWASData$rs),]
    ZMatrixPruned[,  listFiles(IDs = GWAS)  := GWASData[,6]]

    tmp = "Done! \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

  } else if(is.character(GWAS)){  # if external GWAS
    tmp = paste0("# Adding data from the conventional GWAS : \n \"", GWAS, "\" \n")
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

    if(!grepl(".gz", GWAS)){
      GWASData = data.table::fread(GWAS)
    } else if(grepl(".gz", GWAS)) {
      # if tar.gz
      GWASData = data.table::fread(paste0("zcat < ", GWAS))
    }


    rs = match(colnames(GWASData),c("snpid", "snp", "rnpid", "rs"))
    rs = which(!is.na(rs))
    alt = match(colnames(GWASData),c("a1", "alts", "alt"))
    alt = which(!is.na(alt))
    ref = match(colnames(GWASData),c("a2", "ref"))
    ref = which(!is.na(ref))
    z = match(colnames(GWASData),c("z", "Z"))
    z = which(!is.na(z))


    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrixPruned$rs, unlist(GWASData[,..rs])),]
    # check alignment
    aligned = which(GWASData[,..alt, with=F] == ZMatrixPruned$alt &
                      GWASData[,..ref, with=F] == ZMatrixPruned$ref)
    swapped = which(GWASData[,..ref, with=F] == ZMatrixPruned$alt &
                      GWASData[,..alt, with=F] == ZMatrixPruned$ref)
    weird = c(1:nrow(GWASData))[!c(1:nrow(GWASData)) %in% c(aligned, swapped)]

    GWASData[aligned, myZ:= GWASData[aligned, ..z, with=F]]
    GWASData[swapped, myZ:= -GWASData[swapped, ..z, with=F]]
    GWASData[weird, myZ:= NA]


    ZMatrixPruned[, strsplit(GWAS, "/")[[1]][ length(strsplit(GWAS, "/")[[1]])]] = GWASData$myZ

    tmp = "Done! \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)

    }


  ZMatrixPruned = ZMatrixPruned[complete.cases(ZMatrixPruned)]
  tmp = paste0(format(nrow(ZMatrixPruned), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS \n")
  Log = c(Log, tmp)
  if(verbose) cat(tmp)

  # write Pruned ZMatrix
  if(saveFiles){
    write.table(ZMatrixPruned, file="MR_ZMatrix.csv", sep=",", row.names=F, quote=F)
    tmp = "The file \"MR_ZMatrix.csv\" has been successfully created \n"
    Log = c(Log, tmp)
    if(verbose) cat(tmp)
  }
  res=list()
  res$Log = Log
  res$Mat = ZMatrixPruned
  return(res)
}

