###### Function to create pruned Z-matrix of MR instruments ######



# #' Create pruned Z-matrix of MR instruments
# #'
# #' From a list of studies, create pruned Z-matrix of MR instruments (significant at
# #' a specified threshold)
# #'
# #' @inheritParams bGWAS
# #'
# #' @return Log file and pruned Z-Matrix of MR instrument + create a file if saveFiles=T
# #'

makeMR_ZMatrix <- function(prior_studies=NULL, GWAS,
                           MR_threshold=1e-5, path="~/ZMatrices", save_files=F, verbose=F) {
  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = update_log(Log, tmp, verbose)


  if(save_files) Files_Info = data.table::fread("PriorGWASs.tsv")

  if(!is.null(prior_studies)){
    tmp = paste0("Selecting studies :\n")
    Log = update_log(Log, tmp, verbose)
    if(grepl("macOS", sessionInfo()$running)){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=c(1:5, prior_studies+5), showProgress = FALSE)
    } else if(grepl("Linux", sessionInfo()$running)){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=c(1:5, prior_studies+5), showProgress = FALSE)
    } else {
      stop("Only UNIX and MAC OS are supported")
    }

  } else {
    if(grepl("macOS", sessionInfo()$running)) {
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), showProgress = FALSE)
    } else if(grepl("Linux", sessionInfo()$running)){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), showProgress = FALSE)
    } else {
      stop("Only UNIX and MAC OS are supported")
    }
  }

  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = update_log(Log, tmp, verbose)

  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)

  ## 1st step should be taking the SNPs in common, otherwise we might exclude all the SNPs
  # from the conventionnal GWAS when pruning
  # Add conventional GWAS column, at the end (make sure alleles are aligned)
  if(is.numeric(GWAS)){  # if GWAS from our data
    tmp = paste0("# Adding data from the conventional GWAS (ID=", GWAS, "): \n \"", list_files(IDs = GWAS) , "\" \n")
    Log = update_log(Log, tmp, verbose)

    if(grepl("macOS", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5), showProgress = FALSE)
    } else if(grepl("Linux", sessionInfo()$running)){
      GWASData=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5), showProgress = FALSE)
    } else {
      stop("Only UNIX and MAC OS are supported")
    }

    # no need to check for alignment of alleles, just subset and rename the column
    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs,GWASData$rs),]
    ZMatrix[,  list_files(IDs = GWAS)  := GWASData[,6]]

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  } else if(is.character(GWAS)){  # if GWAS is a file
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
    #  weird = c(1:nrow(GWAS))[!c(1:nrow(GWAS)) %in% c(aligned, swapped)]

    GWAS[, myZ:= numeric()]
    GWAS[aligned, myZ:= GWAS[aligned, ..z, with=F]]
    GWAS[swapped, myZ:= -GWAS[swapped, ..z, with=F]]



    ZMatrix[, GName] = GWAS$myZ

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }


  ZMatrix = ZMatrix[complete.cases(ZMatrix)]
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS \n")
  Log = update_log(Log, tmp, verbose)

  # select based on threshold if different from 1e-5 and removed rows without any Z-Score ok
  Zlimit = qnorm(MR_threshold/2, lower.tail = F)

  if(MR_threshold != 1e-5 ){
    tmp = paste0("# Thresholding... \n")
    Log = update_log(Log, tmp, verbose)
    # DO NOT USE THE LAST COLUMN!!
    SNPsToKeep = apply(ZMatrix[,-c(1:5,as.numeric(ncol(ZMatrixPruned)))], 1, function(x) any(abs(x)>Zlimit))
    ZMatrix=ZMatrix[SNPsToKeep,]

    # check that each study have at least one SNP surviving thresholding
    # not needed to check for the last column
    StudiesToKeep = apply(ZMatrix[,-c(1:5,as.numeric(ncol(ZMatrixPruned)))], 2, function(x) any(abs(x)>Zlimit))
    if(!all(StudiesToKeep)){
      tmp = paste0(paste0(colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep], collapse=" - "), " : removed (no strong instrument after thresholding) \n")
      if(save_files){
        Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
          "Excluded for MR: no strong instrument left after thresholding"
      }
      ZMatrix[,colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep] := NULL]
    }
    Log = update_log(Log, tmp, verbose)
    tmp = paste0(ncol(ZMatrix)-6, " studies left after thresholding \n") # -6 and not -5 because the last one is the conv. GWAS
    Log = update_log(Log, tmp, verbose)

    tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs left after thresholding \n")
    Log = update_log(Log, tmp, verbose)
  } else {
    tmp = paste0("No thresholding needed... \n")
    Log = update_log(Log, tmp, verbose)
  }

  # pruning
  # by default, distance pruning, but enable the LD pruning from Aaron's function to compare
  DistancePruning=T
  dist = 100000
  if(DistancePruning){
    tmp = paste0("Distance pruning... \n")
    Log = update_log(Log, tmp, verbose)
    tmp = paste0("   distance : ", dist/1000, "kb \n")
    Log = update_log(Log, tmp, verbose)
    ZMatrixPruned = prune_ZMatrix(ZMatrix, prune.dist = dist)
  } else {
    r2 = 0.8
    tmp = paste0("LD pruning... \n")
    Log = update_log(Log, tmp, verbose)
    tmp = paste0("   distance : ", dist/1000, "kb", " - r2 threshold : ", r2, "\n")
    Log = update_log(Log, tmp, verbose)
    ZMatrixPruned = prune_ZMatrix(ZMatrix, prune.dist = dist, r2.limit=r2)
  }

    # check that each study have at least one SNP surviving pruning
  StudiesToKeep = apply(ZMatrixPruned[,-c(1:5,as.numeric(ncol(ZMatrixPruned))), with=F], 2, function(x) any(abs(x)>Zlimit))
  if(!all(StudiesToKeep)){
    tmp = paste0(paste0(colnames(ZMatrixPruned[,-c(1:5)])[!StudiesToKeep], collapse=" - "), " : removed (no strong instrument after pruning) \n")
    if(save_files){
      Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
        "Excluded for MR: no strong instrument left after pruning"
    }
    ZMatrixPruned[,colnames(ZMatrixPruned[,-c(1:5)])[!StudiesToKeep] := NULL]
  }
  Log = update_log(Log, tmp, verbose)
  tmp = paste0(ncol(ZMatrixPruned)-6, " studies left after pruning \n")
  Log = update_log(Log, tmp, verbose)

  tmp = paste0(format(nrow(ZMatrixPruned), big.mark = ",", scientific = F), " SNPs left after pruning \n")
  Log = update_log(Log, tmp, verbose)


  if(save_files) write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )



  res=list()
  res$log_info = Log
  res$mat = ZMatrixPruned
  return(res)
}

