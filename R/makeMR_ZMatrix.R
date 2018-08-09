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
                           MR_threshold, MR_pruning_dist, MR_pruning_LD, path="~/ZMatrices", save_files=F, verbose=F) {
  platform = .Platform$OS.type
  if(platform=="windows") stop("Windows is not supported yet", call. = FALSE)

  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = update_log(Log, tmp, verbose)


  if(save_files) Files_Info = data.table::fread("PriorGWASs.tsv")

  if(!is.null(prior_studies)){
    tmp = paste0("Selecting studies :\n")
    Log = update_log(Log, tmp, verbose)
    if(platform == "unix"){
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), select=c(1:5, prior_studies+5), showProgress = FALSE, data.table=F)
    }

  } else {
    if(platform == "unix") {
      ZMatrix=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_NotImputed.csv.gz")), showProgress = FALSE, data.table=F)
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
    tmp = paste0("# Adding data from the conventional GWAS (ID=", GWAS, "): \n \"", list_files(IDs = GWAS, Z_matrices = path) , "\" \n")
    Log = update_log(Log, tmp, verbose)

    if(platform == "unix"){
      GWASData=data.table::fread(paste0("zcat < ",paste0(path, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, GWAS+5), showProgress = FALSE, data.table=F)
    }

    # no need to check for alignment of alleles, just subset and rename the column
    # keep the SNPs in our pruned matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs,GWASData$rs),]
    ZMatrix$outcome =  GWASData[,6]
    colnames(ZMatrix)[ncol(ZMatrix)]= list_files(IDs = GWAS, Z_matrices = path)

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  } else if(is.character(GWAS)){  # if GWAS is a file
    tmp = paste0("# Adding data from the conventional GWAS : \n \"", GWAS, "\" \n")
    Log = update_log(Log, tmp, verbose)

    if(!grepl(".gz", GWAS)){
      GWASData = data.table::fread(GWAS, showProgress = FALSE, data.table=F)
    } else if(grepl(".gz", GWAS)) {
      # if tar.gz
      GWASData = data.table::fread(paste0("zcat < ", GWAS), showProgress = FALSE, data.table = F)
    }

    SNPID = match(colnames(GWASData),c("snpid", "snp", "rnpid", "rs", "rsid"))
    SNPID = which(!is.na(SNPID))[1]
    ALT = match(colnames(GWASData),c("a1", "alts", "alt"))
    ALT = which(!is.na(ALT))[1]
    REF = match(colnames(GWASData),c("a2", "ref", "a0"))
    REF = which(!is.na(REF))[1]
    ZSTAT = match(colnames(GWASData),c("z", "Z", "zscore"))
    ZSTAT = which(!is.na(ZSTAT))[1]

    # keep the SNPs in our Z matrix and order them correctly
    GWASData = GWASData[match(ZMatrix$rs, GWASData[,SNPID]),]
    # check alignment
    aligned = which(GWASData[,ALT] == ZMatrix$alt &
                      GWASData[,REF] == ZMatrix$ref)
    swapped = which(GWASData[,REF] == ZMatrix$alt &
                      GWASData[,ALT] == ZMatrix$ref)
    #  weird = c(1:nrow(GWAS))[!c(1:nrow(GWAS)) %in% c(aligned, swapped)]

    GWASData$myZ = NA
    GWASData[aligned, "myZ"] =  GWASData[aligned, ZSTAT]
    GWASData[swapped, "myZ"] = -GWASData[swapped, ZSTAT]

    ZMatrix[, strsplit(GWAS, "/")[[1]][ length(strsplit(GWAS, "/")[[1]])]] = GWASData$myZ

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  } else if(is.data.frame(GWAS)){  # if GWAS is data.frame
    GName = attributes(GWAS)$GName
    tmp = paste0("# Adding data from the conventional GWAS : \n \"", GName,
                 "\" \n")
    Log = update_log(Log, tmp, verbose)

    SNPID = match(colnames(GWAS),c("snpid", "snp", "rnpid", "rs", "rsid"))
    SNPID = which(!is.na(SNPID))[1]
    ALT = match(colnames(GWAS),c("a1", "alts", "alt"))
    ALT = which(!is.na(ALT))[1]
    REF = match(colnames(GWAS),c("a2", "ref", "a0"))
    REF = which(!is.na(REF))[1]
    ZSTAT = match(colnames(GWAS),c("z", "Z", "zscore"))
    ZSTAT = which(!is.na(ZSTAT))[1]

    # keep the SNPs in our Z matrix and order them correctly
    GWAS = GWAS[match(ZMatrix$rs, GWAS[,SNPID]),]
    # check alignment
    aligned = which(GWAS[,ALT] == ZMatrix$alt &
                      GWAS[,REF] == ZMatrix$ref)
    swapped = which(GWAS[,REF] == ZMatrix$alt &
                      GWAS[,ALT] == ZMatrix$ref)
    #  weird = c(1:nrow(GWAS))[!c(1:nrow(GWAS)) %in% c(aligned, swapped)]

    GWAS$myZ = NA
    GWAS[aligned, "myZ"] =  GWAS[aligned, ZSTAT]
    GWAS[swapped, "myZ"] = -GWAS[swapped, ZSTAT]



    ZMatrix[, GName] = GWAS$myZ

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }


  ZMatrix = ZMatrix[complete.cases(ZMatrix),]
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS \n")
  Log = update_log(Log, tmp, verbose)

  # select based on threshold if different from 1e-5 and removed rows without any Z-Score ok

  # ZLimit should be define even if threshold == 1e-5 to remove studies without strong instruments after pruning
  Zlimit = qnorm(MR_threshold/2, lower.tail = F)

  if(MR_threshold != 1e-5 ){
    tmp = paste0("# Thresholding... \n")
    Log = update_log(Log, tmp, verbose)
    # DO NOT USE THE LAST COLUMN!!
    SNPsToKeep = apply(ZMatrix[,-c(1:5,as.numeric(ncol(ZMatrix)))], 1, function(x) any(abs(x)>Zlimit))
    ZMatrix=ZMatrix[SNPsToKeep,]

    # check that each study have at least one SNP surviving thresholding
    # not needed to check for the last column
    StudiesToKeep = apply(ZMatrix[,-c(1:5,as.numeric(ncol(ZMatrix)))], 2, function(x) any(abs(x)>Zlimit))
    if(!all(StudiesToKeep)){
      tmp = paste0(paste0(colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep], collapse=" - "), " : removed (no strong instrument after thresholding) \n")
      if(save_files){
        Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
          "Excluded for MR: no strong instrument left after thresholding"
      }
      ZMatrix[,colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] <- NULL
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
  tmp = paste0("Pruning MR instruments... \n")
  Log = update_log(Log, tmp, verbose)
  apply(ZMatrix[,-c(1:5, ncol(ZMatrix))], 1, function(ZMatrix){
    max(abs(ZMatrix))
  }) -> maxZ
  ToPrune = ZMatrix[,1:3]
  colnames(ToPrune) = c("SNP", "chr_name", "chr_start")
  if(MR_pruning_LD>0){# LD-pruning
    tmp = paste0("   distance : ", MR_pruning_dist, "Kb", " - LD threshold : ", MR_pruning_LD, "\n")
    Log = update_log(Log, tmp, verbose)
    # get max Z-score and convert to p-value
    ToPrune$pval.exposure = 2*pnorm(-abs(maxZ))
    # Do pruning, chr by chr
    SNPsToKeep = c()
    for(chr in unique(ToPrune$chr_name)){
      SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name==chr,], clump_kb = MR_pruning_dist, clump_r2 = MR_pruning_LD)$SNP))
      }
    }
  else {# distance pruning
    tmp = paste0("   distance : ", MR_pruning_dist, "Kb \n")
    Log = update_log(Log, tmp, verbose)
    ToPrune$Z = maxZ
    SNPsToKeep = prune_byDistance(ToPrune, prune.dist=MR_pruning_dist, byP=F)
  }


  ZMatrixPruned = ZMatrix[ZMatrix$rs %in% SNPsToKeep,]

  # check that each study have at least one SNP surviving pruning
  StudiesToKeep = apply(ZMatrixPruned[,-c(1:5,as.numeric(ncol(ZMatrixPruned)))], 2, function(x) any(abs(x)>Zlimit))
  if(!all(StudiesToKeep)){
    tmp = paste0(paste0(colnames(ZMatrixPruned[,-c(1:5)])[!StudiesToKeep], collapse=" - "), " : removed (no strong instrument after pruning) \n")
    if(save_files){
      Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
        "Excluded for MR: no strong instrument left after pruning"
    }
    ZMatrixPruned[,colnames(ZMatrixPruned[,-c(1:5)])[!StudiesToKeep]] <- NULL
  }
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

