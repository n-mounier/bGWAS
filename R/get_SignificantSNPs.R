###### Function to identify significant SNPs ######



# #' Identify significant SNPs
# #'
# #' From BFs and p-values, identify significant SNPs given a specified threshold and prune them
# #' if necessary
# #' @inheritParams bGWAS
# #' @param Prior from the function request_BFandP()
# #'
# #' @return Log + data.table containing rs-chr-pos-alt-ref-obs-fit-se-...
# #' @export
# # Function not exported, no need of extended documentation?


get_significantSNPs <- function(Prior, sign_method="p", sign_thresh=5e-8, res_pruning_dist, res_pruning_LD, save_files=F, verbose=F) {
  Log = c()

  "%S>%" <- function(x,to.be.ignored){ # like a 'sink' version of magrittr's  %T>%
    x
  }


  if(!sign_method %in% c("p", "fdr")) stop("method not accepted",  call. = FALSE)
  if(!is.numeric(sign_thresh)) stop("non numeric threshold",  call. = FALSE)

  tmp = paste0("   Starting with ", format(nrow(Prior), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)


  if(sign_method == "fdr"){
    tmp = "# Selecting significant SNPs according to FDR (Benjamini-Hochberg procedure)... \n"
    Log = update_log(Log, tmp, verbose)


    Prior$fdr = p.adjust(Prior$BF_P, method='fdr')
    PriorThr = subset(Prior, fdr<sign_thresh)

    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)



    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }

  if(sign_method == "p"){
    tmp = "# Selecting significant SNPs according to p-values... \n"
    Log = update_log(Log, tmp, verbose)


    PriorThr = subset(Prior, BF_P<sign_thresh)

    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)


    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }



  # pruning results only if res_pruning_dist!=0
  if(!is.null(res_pruning_dist)){
    tmp = "# Pruning significant SNPs... \n"
    Log = update_log(Log, tmp, verbose)
    ToPrune = PriorThr[,c(1,2,3,12)]
    colnames(ToPrune) = c("SNP", "chr_name", "chr_start", "pval.exposure")
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
    PriorThr =  PriorThr[PriorThr$rs %in% SNPsToKeep,]

    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)

    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)
  }




  if(save_files){
    write.table(PriorThr, file= "SignificantSNPs.csv", sep=",", row.names=F, quote=F)
    tmp = "The file \"SignificantSNPs.csv\" has been successfully created \n"
    Log = update_log(Log, tmp, verbose)

  }


  res=list()
  res$log_info = Log
  res$SNPs = PriorThr$rs
  return(res)
}
