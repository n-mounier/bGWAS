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


get_significantSNPs <- function(Prior, sign_method="p", sign_thresh=5e-8, prune_res=T, save_files=F, verbose=F) {
  Log = c()

  "%S>%" <- function(x,to.be.ignored){ # like a 'sink' version of magrittr's  %T>%
    x
  }


  if(!sign_method %in% c("p", "fdr")) stop("method not accepted")
  if(!is.numeric(sign_thresh)) stop("non numeric threshold")

  tmp = paste0("   Starting with ", format(nrow(Prior), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)


  if(sign_method == "fdr"){
    tmp = "# Selecting significant SNPs according to FDR (Benjamini-Hochberg procedure)... \n"
    Log = update_log(Log, tmp, verbose)


    Prior[, fdr := p.adjust(BF_p, method='fdr')]
    PriorThr = Prior[fdr<sign_thresh]

    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)



    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }

  if(sign_method == "p"){
    tmp = "# Selecting significant SNPs according to p-values... \n"
    Log = update_log(Log, tmp, verbose)


    PriorThr = Prior[BF_p<sign_thresh]

    tmp = paste0(format(nrow(PriorThr), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)


    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  }

  if(prune_res){
    tmp = "# Pruning significant SNPs... \n"
    Log = update_log(Log, tmp, verbose)


    # snps should be ordered by significance
    accepted.snps = ""[0] # an empty character vector of the 'rs' names that will survive pruning
    i = 1
    while(i <= nrow(PriorThr)) {
      # 'accept' the current SNP
      PriorThr[i] -> accepted
      accepted.rs = accepted$rs

      accepted.snps = c(accepted.snps, accepted.rs)


      # now, skip on until we find the next one to accept (or we reach the end)
      while(T) {

        i = i+1

        if(i >  nrow(PriorThr)) {
          break # all done
        }

        # check if the next (candidate) 'current' SNP is close to earlier accepted SNPs
        cand = PriorThr[i]

        cand$rs -> cand.rs
        cand$chr -> cand.chrm
        cand$pos -> cand.pos

        PriorThr[PriorThr$chr == cand.chrm & PriorThr$rs %in% accepted.snps ] %S>% print -> same.chrm.already.accepted


        if(nrow(same.chrm.already.accepted) == 0) { break } # no problem with this one, just restart and allow it to be accepted
        same.chrm.already.accepted[, pos.delta := pos-cand.pos ]
        same.chrm.already.accepted = same.chrm.already.accepted[ abs(pos.delta) < 500000 ] # no point looking far away
        if(nrow(same.chrm.already.accepted) == 0) { break }
        next

      }
    }
    # select Significant SNPs
    SignifSNPs = PriorThr[PriorThr$rs %in% accepted.snps,]

    tmp = paste0(format(nrow(SignifSNPs), big.mark = ",", scientific = F), " SNPs left \n")
    Log = update_log(Log, tmp, verbose)


    tmp = "Done! \n"
    Log = update_log(Log, tmp, verbose)

  } else {
    SignifSNPs = PriorThr
  }




  if(save_files){
    write.table(SignifSNPs, file= "SignificantSNPs.csv", sep=",", row.names=F, quote=F)
    tmp = "The file \"SignificantSNPs.csv\" has been successfully created \n"
    Log = update_log(Log, tmp, verbose)

  }


  res=list()
  res$log_info = Log
  res$SNPs = SignifSNPs$rs
  return(res)
}
