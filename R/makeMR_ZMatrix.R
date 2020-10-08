###### Function to create pruned Z-matrix of MR instruments ######



# #' Create pruned Z-matrix of MR instruments
# #'
# #' From a list of studies, create pruned Z-matrix of MR instruments (significant at
# #' a specified threshold)
# #'
# #' @inheritParams bGWAS
# NOT EXPORTED


makeMR_ZMatrix <- function(prior_studies=NULL, GWASData, GName,
                           MR_threshold, MR_ninstruments, MR_pruning_dist, MR_pruning_LD, Z_matrices="~/ZMatrices", 
                           MR_shrinkage, save_files=F, verbose=F) {
  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = update_log(Log, tmp, verbose)
  
  
  if(save_files) Files_Info = readr::read_tsv("PriorGWASs.tsv", progress = FALSE, col_types = readr::cols())
  
  if(!is.null(prior_studies)){
    tmp = paste0("Selecting studies :\n")
    Log = update_log(Log, tmp, verbose)
    ZMatrix = as_tibble(data.table::fread(file.path(Z_matrices, "ZMatrix_MR.csv.gz"), select=c(1:5, prior_studies+5), showProgress = FALSE, data.table=F))
  } else {
    ZMatrix = as_tibble(data.table::fread(file.path(Z_matrices, "ZMatrix_MR.csv.gz"), showProgress = FALSE, data.table=F))
  }
  
  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = update_log(Log, tmp, verbose)
  
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)
  
  ## 1st step should be taking the SNPs in common, otherwise we might exclude all the SNPs
  # from the conventionnal GWAS when pruning
  # Add conventional GWAS column, at the end (make sure alleles are aligned)
  tmp = paste0("# Adding data from the conventional GWAS : \n \"", GName,
               "\" \n")
  Log = update_log(Log, tmp, verbose)
  
  # keep the SNPs in our Z matrix and order them correctly + check allele alignement
  ZMatrix %>%
    filter(.data$rs %in% GWASData$rsid) -> ZMatrix
  
  GWASData %>%
    slice(match(ZMatrix$rs, .data$rsid)) %>%
    mutate(ref = toupper(ref),
           alt = toupper(alt),
           ZMat_alt = toupper(ZMatrix$alt),
           Zmat_ref = toupper(ZMatrix$ref),
           aligned_Z = case_when(
             (.data$alt == .data$ZMat_alt &
                .data$ref == .data$Zmat_ref) ~ .data$z_obs,
             (.data$ref == .data$ZMat_alt &
                .data$alt == .data$Zmat_ref) ~ -.data$z_obs,
             TRUE ~ NA_real_))-> GWASData
  
  ZMatrix %>%
    mutate({{GName}} := GWASData$aligned_Z) -> ZMatrix
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  
  ZMatrix %>%
    tidyr::drop_na() -> ZMatrix
  
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS \n")
  Log = update_log(Log, tmp, verbose)
  
  # select based on threshold and remove rows without any Z-Score ok
  
  # ZLimit should be define even if threshold == 1e-5 to remove studies without strong instruments after pruning
  Zlimit = stats::qnorm(MR_threshold/2, lower.tail = F)
  
  tmp = paste0("# Thresholding... \n")
  Log = update_log(Log, tmp, verbose)
  # DO NOT USE THE LAST COLUMN!!
  ZMatrix %>%
    filter_at(vars(-c(1:5,as.numeric(ncol(ZMatrix)))), 
              any_vars(abs(.data$.) > Zlimit)) -> ZMatrix
  
  
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs left after thresholding \n")
  Log = update_log(Log, tmp, verbose)
  
  # before pruning, remove studies with less than MR_ninstruments instruments
  # check that each study have at least two SNPs surviving pruning+thresholding
  ZMatrix %>%
    select(-c(1:5,as.numeric(ncol(ZMatrix)))) %>%
    apply(MARGIN = 2, FUN = function(col){
      sum(abs(col)>Zlimit)>=MR_ninstruments})  -> StudiesToKeep

  
  if(!all(StudiesToKeep)){
    ZMatrix  %>%
      select(-c(1:5,as.numeric(ncol(ZMatrix)))) %>%
      select(names(StudiesToKeep[!StudiesToKeep])) %>%
      colnames -> StudiesToRemove
    tmp = paste0(paste0(get_names(StudiesToRemove, Z_matrices), collapse=" - "), " : removed (less than ", MR_ninstruments, " instrument after thresholding) \n")
    Log = update_log(Log, tmp, verbose)
    
    if(save_files){
      Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
        paste0("Excluded for MR: less than ",  MR_ninstruments, " strong instrument left after thresholding/pruning")
    }
    
    ZMatrix %>%
      select(-StudiesToRemove) -> ZMatrix
    
  }
  tmp = paste0(ncol(ZMatrix)-6, " studies left after thresholding \n")
  Log = update_log(Log, tmp, verbose)
  
  
  
  # pruning
  tmp = paste0("Pruning MR instruments... \n")
  Log = update_log(Log, tmp, verbose)
  
  ZMatrix %>%
    select(-c(1:5,as.numeric(ncol(ZMatrix)))) %>%
    apply(MARGIN = 1, FUN = function(row){
      max(abs(row))}) -> maxZ
 
  
  ZMatrix %>%
    transmute(SNP=.data$rs,
              chr_name=.data$chrm,
              chr_start=.data$pos,
              Z=maxZ) -> ToPrune

  if(MR_pruning_LD>0){# LD-pruning
    tmp = paste0("   distance : ", MR_pruning_dist, "Kb", " - LD threshold : ", MR_pruning_LD, "\n")
    Log = update_log(Log, tmp, verbose)
    # convert max Z-score to p-value
    ToPrune %>%
      mutate(pval.exposure = 2 * stats::pnorm(-abs(.data$Z)),
             Z=NULL) -> ToPrune
    # Do pruning, chr by chr
    SNPsToKeep = c()
    for(chr in unique(ToPrune$chr_name)){
      SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name==chr,], clump_kb = MR_pruning_dist, clump_r2 = MR_pruning_LD)$SNP))
    }
  } else{# distance pruning
    tmp = paste0("   distance : ", MR_pruning_dist, "Kb \n")
    Log = update_log(Log, tmp, verbose)
    SNPsToKeep = prune_byDistance(ToPrune, prune.dist=MR_pruning_dist, byP=F)
  }
  
  ZMatrix %>%
    filter(.data$rs %in% SNPsToKeep) -> ZMatrixPruned

  tmp = paste0(format(nrow(ZMatrixPruned), big.mark = ",", scientific = F), " SNPs left after pruning \n")
  Log = update_log(Log, tmp, verbose)
  
  NAllStudies = nrow(ZMatrixPruned)
  
  # check that each study have at least MR_ninstruments SNPs surviving pruning+thresholding
  ZMatrixPruned %>%
    select(-c(1:5,as.numeric(ncol(ZMatrixPruned)))) %>%
    apply(MARGIN = 2, FUN = function(col){
      sum(abs(col)>Zlimit)>=MR_ninstruments}) -> StudiesToKeep
    
  
  if(!all(StudiesToKeep)){
    ZMatrixPruned  %>%
      select(-c(1:5,as.numeric(ncol(ZMatrixPruned)))) %>%
      select(names(StudiesToKeep[!StudiesToKeep])) %>%
      colnames -> StudiesToRemove
    tmp = paste0(paste0(get_names(StudiesToRemove, Z_matrices), collapse=" - "), " : removed (less than ", MR_ninstruments, " strong instrument after pruning) \n")
    Log = update_log(Log, tmp, verbose)
    
    if(save_files){
      Files_Info$status[Files_Info$File %in% colnames(ZMatrix[,-c(1:5)])[!StudiesToKeep]] =
        paste0("Excluded for MR: less than", MR_ninstruments,  " strong instrument left after thresholding/pruning")
      }

    ZMatrixPruned %>%
      select(-StudiesToRemove) -> ZMatrixPruned
    
  }
  tmp = paste0(ncol(ZMatrixPruned)-6, " studies left after thresholding+pruning \n")
  Log = update_log(Log, tmp, verbose)
  
  
  # Further checking of the SNPs to remove SNPs associated with studies removed because only one SNP
  ZMatrixPruned %>%
    select(-c(1:5,as.numeric(ncol(ZMatrixPruned)))) %>%
    apply(MARGIN = 1, FUN = function(row){
      any(abs(row)>Zlimit)}) -> SNPsToKeep  
  
    
  if(sum(SNPsToKeep) != NAllStudies){
    ZMatrixPruned %>%
      filter(SNPsToKeep) -> ZMatrixPruned
    tmp = paste0(format(nrow(ZMatrixPruned), big.mark = ",", scientific = F), " SNPs left after removing studies with only one strong instrument \n")
    Log = update_log(Log, tmp, verbose)
  }
  
  

  push_extreme_zs_back_a_little_towards_zero <- function(d) { # Some z-scores are just too far from zero
    maxAllowed_z = abs(stats::qnorm(1e-300 / 2)) # p=1e-300 is the max allowed now, truncate z-scores accordingly
    names(d)[!names(d) %in% c("rs","chrm","pos","alt","ref")] -> studies_here
    for(n in studies_here) {
      d %>%
        mutate(!!n := case_when(
          abs(eval(parse(text=n))) > maxAllowed_z  ~  maxAllowed_z,
          abs(eval(parse(text=n))) < -maxAllowed_z ~ -maxAllowed_z,
          TRUE ~ eval(parse(text=n)))) -> d
    }
    return(d)
  }
  # Truncate Z-scores
  ZMatrixPruned %>%
    push_extreme_zs_back_a_little_towards_zero() -> ZMatrixPruned
  
  
  
  # Set the z-scores to 0 for the regression if shrinkage
  if(MR_shrinkage < 1.0) {
    names(ZMatrixPruned)[!names(ZMatrixPruned) %in% c('rs','chrm','pos','alt','ref', GName)] -> Prior_study_names
    threshold = abs(stats::qnorm(MR_shrinkage/2))
    for(column_of_zs in Prior_study_names) { 
      ZMatrixPruned %>%
        mutate(!!column_of_zs := case_when(
          abs(eval(parse(text=column_of_zs))) < threshold ~ 0,
          TRUE ~ eval(parse(text=column_of_zs)))) -> ZMatrixPruned
    }
    tmp = paste0("Applying shrinkage (threshold = ", MR_shrinkage, ") before performing MR. \n")
    Log = update_log(Log, tmp, verbose)
  }
  
  
  if(save_files) utils::write.table(Files_Info, file="PriorGWASs.tsv", sep="\t", quote=F, row.names=F )
  
  
  
  res=list(log_info = Log,
           mat = ZMatrixPruned)
  return(res)
}


