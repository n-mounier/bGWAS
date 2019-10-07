###### Function to identify significant SNPs ######



# #' Identify significant SNPs
# #'
# #' From BFs and p-values, identify significant SNPs given a specified threshold and prune them
# #' if necessary
# #' @inheritParams bGWAS
# NOT EXPORTED


get_significantSNPs <- function(Prior, sign_method="p", sign_thresh=5e-8, res_pruning_dist, res_pruning_LD, 
                                selected_studies, All_ZMatrix, save_files=F, verbose=F) {
  Log = c()

  
  
  if(!sign_method %in% c("p", "fdr")) stop("method not accepted",  call. = FALSE)
  if(!is.numeric(sign_thresh)) stop("non numeric threshold",  call. = FALSE)
  
  
  
  #### BFs ####
  tmp = "Identification based on BFs \n"
  Log = update_log(Log, tmp, verbose)
  
  
  tmp = paste0("   Starting with ", format(nrow(Prior), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)
  
  
  if(sign_method == "fdr"){
    tmp = "# Selecting significant SNPs according to FDR (Benjamini-Hochberg procedure)... \n"
    Log = update_log(Log, tmp, verbose)
    
    Prior %>%
      filter(.data$BF_fdr<sign_thresh) -> PriorThr

    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  if(sign_method == "p"){
    tmp = "# Selecting significant SNPs according to p-values... \n"
    Log = update_log(Log, tmp, verbose)
    
    
    Prior %>%
      filter(.data$BF_p<sign_thresh) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  
  if(nrow(PriorThr)>0){
  # pruning results only if res_pruning_dist!=0
  if(!is.null(res_pruning_dist)){
    tmp = "# Pruning significant SNPs... \n"
    Log = update_log(Log, tmp, verbose)
    # create rsid, chr, pos, p data.frame (all columns renamed in tidyinput, ok to use rsid directly)
    PriorThr %>% 
      transmute(SNP = .data$rsid,
                chr_name = .data$chrm_UK10K,
                chr_start = .data$pos_UK10K,
                pval.exposure = .data$BF_p) ->ToPrune 
    if(res_pruning_LD>0){# LD-pruning
      tmp = paste0("   distance : ", res_pruning_dist, "Kb", " - LD threshold : ", res_pruning_LD, "\n")
      Log = update_log(Log, tmp, verbose)
      # Do pruning, chr by chr
      SNPsToKeep = c()
      for(chr in unique(ToPrune$chr_name)){
        SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name==chr,], clump_kb = res_pruning_dist, clump_r2 = res_pruning_LD)$SNP))
      }
    } else {# distance pruning
      tmp = paste0("   distance : ", res_pruning_dist, "Kb \n")
      Log = update_log(Log, tmp, verbose)
      SNPsToKeep = prune_byDistance(ToPrune, prune.dist=res_pruning_dist, byP=T)
    }
    PriorThr %>%
      filter(.data$rsid %in% SNPsToKeep) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
  }
  
  
  BF_SNPs = pull(PriorThr, .data$rsid)
  
  if(save_files){
    readr::write_csv(path="SignificantSNPs.csv",  x=PriorThr)
    tmp = "The file \"SignificantSNPs.csv\" has been successfully created \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  }
  
  #### posterior ####
  tmp = "Identification based on posterior effects \n"
  Log = update_log(Log, tmp, verbose)
  
  
  tmp = paste0("   Starting with ", format(nrow(Prior), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)
  
  
  if(sign_method == "fdr"){
    tmp = "# Selecting significant SNPs according to FDR (Benjamini-Hochberg procedure)... \n"
    Log = update_log(Log, tmp, verbose)
  
    
    Prior %>%
      filter(.data$fdr_posterior<sign_thresh) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  if(sign_method == "p"){
    tmp = "# Selecting significant SNPs according to p-values... \n"
    Log = update_log(Log, tmp, verbose)
    
    
    Prior %>%
      filter(.data$p_posterior<sign_thresh) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  
  
  # pruning results only if res_pruning_dist!=0
  print(nrow(PriorThr))
 if(nrow(PriorThr)>0){
  if(!is.null(res_pruning_dist)){
    tmp = "# Pruning significant SNPs... \n"
    Log = update_log(Log, tmp, verbose)
    # create rsid, chr, pos, p data.frame (all columns renamed in tidyinput, ok to use rsid directly)
    PriorThr %>% 
      transmute(SNP = .data$rsid,
                chr_name = .data$chrm_UK10K,
                chr_start = .data$pos_UK10K,
                pval.exposure = .data$p_posterior) ->ToPrune 
    if(res_pruning_LD>0){# LD-pruning
      tmp = paste0("   distance : ", res_pruning_dist, "Kb", " - LD threshold : ", res_pruning_LD, "\n")
      Log = update_log(Log, tmp, verbose)
      # Do pruning, chr by chr
      SNPsToKeep = c()
      for(chr in unique(ToPrune$chr_name)){
        SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name==chr,], clump_kb = res_pruning_dist, clump_r2 = res_pruning_LD)$SNP))
      }
    } else {# distance pruning
      tmp = paste0("   distance : ", res_pruning_dist, "Kb \n")
      Log = update_log(Log, tmp, verbose)
      SNPsToKeep = prune_byDistance(ToPrune, prune.dist=res_pruning_dist, byP=T)
    }
    PriorThr %>%
      filter(.data$rsid %in% SNPsToKeep) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
  }
 
  
  posterior_SNPs = pull(PriorThr, .data$rsid)
 }
  
  #### direct ####
  tmp = "Identification based on direct effects \n"
  Log = update_log(Log, tmp, verbose)
  
  
  tmp = paste0("   Starting with ", format(nrow(Prior), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)
  
  
  if(sign_method == "fdr"){
    tmp = "# Selecting significant SNPs according to FDR (Benjamini-Hochberg procedure)... \n"
    Log = update_log(Log, tmp, verbose)
    
    
    Prior %>%
      filter(.data$fdr_direct<sign_thresh) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  if(sign_method == "p"){
    tmp = "# Selecting significant SNPs according to p-values... \n"
    Log = update_log(Log, tmp, verbose)
    
    
    Prior %>%
      filter(.data$p_direct<sign_thresh) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
    
  }
  
  
  
  # pruning results only if res_pruning_dist!=0
  if(nrow(PriorThr)>0){
  if(!is.null(res_pruning_dist)){
    tmp = "# Pruning significant SNPs... \n"
    Log = update_log(Log, tmp, verbose)
    # create rsid, chr, pos, p data.frame (all columns renamed in tidyinput, ok to use rsid directly)
    PriorThr %>% 
      transmute(SNP = .data$rsid,
                chr_name = .data$chrm_UK10K,
                chr_start = .data$pos_UK10K,
                pval.exposure = .data$p_direct) ->ToPrune 
    if(res_pruning_LD>0){# LD-pruning
      tmp = paste0("   distance : ", res_pruning_dist, "Kb", " - LD threshold : ", res_pruning_LD, "\n")
      Log = update_log(Log, tmp, verbose)
      # Do pruning, chr by chr
      SNPsToKeep = c()
      for(chr in unique(ToPrune$chr_name)){
        SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name==chr,], clump_kb = res_pruning_dist, clump_r2 = res_pruning_LD)$SNP))
      }
    } else {# distance pruning
      tmp = paste0("   distance : ", res_pruning_dist, "Kb \n")
      Log = update_log(Log, tmp, verbose)
      SNPsToKeep = prune_byDistance(ToPrune, prune.dist=res_pruning_dist, byP=T)
    }
    PriorThr %>%
      filter(.data$rsid %in% SNPsToKeep) -> PriorThr
    
    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)
    
    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
  }
  
  
  direct_SNPs = pull(PriorThr, .data$rsid)
  }
  
  
  
  # To create heatmap, need ZMat for significant SNPs
  All_ZMatrix %>%
    select(c(1:5), selected_studies) %>%
    filter(.data$rs %in% BF_SNPs) -> Matrix_Heatmap
  
  res=list()
  res$log_info = Log
  res$SNPs = ifelse(exists("BF_SNPs"), BF_SNPs, NA)
  res$posterior= ifelse(exists("posterior_SNPs"), posterior_SNPs, NA)
  res$direct = ifelse(exists("direct_SNPs"), direct_SNPs, NA)
  res$mat = Matrix_Heatmap
  return(res)
}
