###### bGWAS class functions ######



#' Print a bGWAS object
#' @param x an object of class bGWAS
#' @param ... further arguments
#' 
#' @return print
#' @export
print.bGWAS <- function(x,...) {
  cat("-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n ")
  cat (paste0(" Analysis : \"",strsplit(strsplit(x$log_info[
    grep("The name of your analysis is: ", x$log_info)],
    "The name of your analysis is: \"", fixed=T)[[1]][2], "\"", fixed=T)[[1]][1],"\" \n"))
  cat( "bGWAS performed on", format( nrow(x$all_BFs) , big.mark = "," , scientific = F) ,
                                    "SNPs \n \n" )

  cat("-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n")

  if(nrow(x$significant_studies)>1){
    cat(nrow(x$significant_studies), "studies used to build the prior : \n")
    print(x$significant_studies[,1:3], row.names=F)
  } else {
    cat(nrow(x$significant_studies), "study used to build the prior : \n")
    print(x$significant_studies[,1:3], row.names=F)
  }

  cat("\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n")

  if(length(x$significant_SNPs)==0){
    cat("No significant SNP identified, because the analysis has been limited to prior estimation")
  } else
    if(length(x$significant_SNPs)==1){
      cat(length(x$significant_SNPs), "significant SNP identified : \n ")
      cat(x$significant_SNPs, sep=", ")
    } else {
      cat(length(x$significant_SNPs), "significant SNPs identified : \n ")
      cat(x$significant_SNPs, sep=", ")
    }
  cat("\n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ ")

  }




#' Equality test for bGWAS objects
#' @param obj1 an object of class bGWAS
#' @param obj2 an object of class bGWAS
#'
#' @return all.equal
#' @export
all.equal.bGWAS <- function(obj1, obj2) {
  if(length(obj1)!= length(obj2)) return(FALSE)
  # log - we don't care.
  # significant SNPs
  if(!all.equal(obj1$significant_SNPs,obj2$significant_SNPs)) return(FALSE)
  # all BFs
  if(!all.equal(obj1$all_BFs,obj2$all_BFs, tolerance=0.0001)) return(FALSE)
  # significant studies
  if(!all.equal(obj1$significant_studies,obj2$significant_studies, tolerance=0.0001)) return(FALSE)
  # all MR coeffs
  if(!all.equal(obj1$all_MRcoeffs,obj2$all_MRcoeffs, tolerance=0.0001)) return(FALSE)

  return(TRUE)
}




#' Manhattan Plot from bGWAS results
#'
#' Creates a Manhattan Plot from bGWAS results (for performance, only SNPs with 
#' p-value or FDR < 0.05 are plotted)
#'
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#' @param save_file A logical indicating if the graphic should be saved,
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device
#' @param file_name The name of the file saved (is \code{save_file} is \code{TRUE})
#'        \code{default=NULL}, will used NameOfYourAnalysis_ManhattanPlot.png
#' @param annotate A logical indicating if the significant SNPs identified in the
#'        analysis should be annotated on the plot, \code{default=TRUE}
#'        If your results are not pruned or if you have a high number of significant SNPs,
#'        be aware that \code{annotate=TRUE} might decrease readability of the figure.
#' @param results, "BF" / "posterior" / "direct",  \code{default="BF"}
#' 
#' @details
#' If \code{results = "BF"}, BF p-values / fdr-values will be used. \cr
#' If \code{results = "direct"}, direct effect p-values / fdr-values will be used. \cr
#' If \code{results = "posterior"}, posterior effect p-values / fdr-values will be used. \cr

#' @return a Manhattan Plot
#'
#' @export

manhattan_plot_bGWAS <- function(obj, save_file=F, file_name=NULL,
                                  annotate=T, results="BF") {

  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(length(obj$significant_SNPs)==0){
    stop("Manhattan plot can't be displayed, because the analysis has been limited to prior estimation", call. = F)
  }
  if(!is.logical(save_file)) stop("save_file : should be logical")
  # if no name, use the one from the analysis (in log file)
  if(save_file && is.null(file_name)){
    file_name = paste0(strsplit(strsplit(obj$log_info[
      grep("The name of your analysis is: ", obj$log_info)],
      "The name of your analysis is: \"", fixed=T)[[1]][2], "\"", fixed=T)[[1]][1],
    "_ManhattanPlot.png")
  }
  if(save_file && !is.character(file_name)) stop("file_name : should be a character")
  if(!is.logical(annotate)) stop("annotate : shoud be logical")
  method = strsplit(strsplit(obj$log_info[
    grep("Significant SNPs will be identified according to ", obj$log_info)],
    "Significant SNPs will be identified according to ")[[1]][2], ".", fixed=T)[[1]][1]
  threshold = as.numeric(strsplit(strsplit(obj$log_info[
    grep("Significant SNPs will be identified according to ", obj$log_info)],
    "The threshold used is :")[[1]][2], ".", fixed=T)[[1]][1])


  
  if(results == "BF"){
    value = ifelse(method=="FDR", "BF_fdr", "BF_p")
    obj$all_BFs %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, {{value}}) %>%
      set_names(c("rsid", "chrm_UK10K", "pos_UK10K", "p")) %>%
      filter(.data$p<0.05) -> ToPlot
    my_ylab=ifelse( method == "FDR" ,
                 expression(-log[10](italic(fdr))) ,
                 expression(-log[10](italic(p))) )
  } else if(results == "posterior"){
    value = ifelse(method=="FDR", "fdr_posterior", "p_posterior")
    obj$all_BFs %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, {{value}}) %>%
      set_names(c("rsid", "chrm_UK10K", "pos_UK10K", "p")) %>%
      filter(.data$p<0.05) -> ToPlot
    my_ylab=ifelse( method == "FDR" ,
                    expression(-log[10](italic(fdr))~-~posterior~effects) ,
                    expression(-log[10](italic(p))~-~posterior~effects))
  } else if(results == "direct"){
    value = ifelse(method=="FDR", "fdr_direct", "p_direct")
    obj$all_BFs %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, {{value}}) %>%
      set_names(c("rsid", "chrm_UK10K", "pos_UK10K", "p")) %>%
      filter(.data$p<0.05) -> ToPlot
    my_ylab=ifelse( method == "FDR" ,
                    expression(-log[10](italic(fdr))~-~direct~effects) ,
                    expression(-log[10](italic(p))~-~direct~effects))
  }


  if(save_file) grDevices::png(file_name, width = 20, height = 12, units = "cm", res = 500)

  qqman::manhattan(ToPlot ,
                   chr = "chrm_UK10K" ,
                   bp = "pos_UK10K" ,
                   p = "p" ,
                   snp = "rsid",
                   col = c("gray10", "gray60"),
                   cex = 0.6,
                   suggestiveline = FALSE,
                   genomewideline = -log10(threshold),
                   highlight = NULL,
                   logp = TRUE,
                   annotatePval = NULL,
                   annotateTop = FALSE,
                   ylab= my_ylab)



  if(annotate){ # significant SNPs from the analysis
    # extract them
    if(results=="BF"){
      ToPlot %>%
        filter( .data$rsid %in% obj$significant_SNPs) -> SNPs_to_plot
    } else if(results=="posterior"){
      ToPlot %>%
        filter( .data$rsid %in% obj$posterior_SNPs) -> SNPs_to_plot
    } else if(results=="direct"){
      ToPlot %>%
        filter( .data$rsid %in% obj$direct_SNPs) -> SNPs_to_plot
    } 
    
    
    
    ToPlot %>%
      mutate(p=NULL) %>%
      arrange(.data$chrm_UK10K,.data$pos_UK10K) -> all

    # y = -log10(p) ou -log10(fdr)
    SNPs_to_plot %>%
      mutate(y = -log10(.data$p) - 0.3) -> SNPs_to_plot # to make sure all SNPs names are in plotting windows
     # x = have to look at all SNPs chr/pos
    get_posx <- function(snp, all){
      chr = as.numeric(snp[2])
      pos = as.numeric(snp[3])
      if(chr==1){
        posx = pos
      }
      else{
        posx=0
        for(ch in 1:(chr-1)){
          posx = posx + utils::tail(all %>% filter(.data$chrm_UK10K==ch) %>% pull(.data$pos_UK10K),1)
        }
        posx = posx + pos
      }
      return(posx)
    }

    SNPs_to_plot$x = apply(SNPs_to_plot, 1, function(x) get_posx(x, all))

    calibrate::textxy(SNPs_to_plot$x, SNPs_to_plot$y,
                      labs=SNPs_to_plot$rsid,
                      offset = 0.625, cex=0.45)
  }


  if(save_file) grDevices::dev.off()
}




#' Coefficients Plot from bGWAS results
#'
#' Creates a Coefficients Plot (causal effect of each Prior GWASs) 
#'
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#' @param save_file A logical indicating if the graphic should be saved,
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device
#' @param file_name The name of the file saved (is \code{save_file} is \code{TRUE})
#'        \code{default=NULL}, will used NameOfYourAnalysis_CoefficientsPlot.png

#' @return a Coefficients Plot
#' @export

coefficients_plot_bGWAS <- function(obj, save_file=F, file_name=NULL){
  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(is.null(obj$all_MRcoeffs)) stop("The prior has not been created using Prior GWASs, there are no coefficients to plot")
  if(!is.logical(save_file)) stop("save_file : should be logical")
  # if no name, use the one from the analysis (in log file)
  if(save_file && is.null(file_name)){
    file_name = paste0(strsplit(strsplit(obj$log_info[
      grep("The name of your analysis is: ", obj$log_info)],
      "The name of your analysis is: \"", fixed=T)[[1]][2], "\"", fixed=T)[[1]][1],
      "_CoefficientsPlot.png")
  }
  if(save_file && !is.character(file_name)) stop("file_name : should be a character")


  coeffs = obj$significant_studies
  # add the trait name (if multiple studies for a same trait, add " - X" )
  Z_matrices = strsplit(obj$log_info[stringr::str_detect(obj$log_info, "The Z-Matrix files are stored in \"")], "\"")[[1]][2]
  all_studies =  list_priorGWASs(Z_matrices = Z_matrices)
  coeffs %>%
    mutate(Trait=get_names(.data$study, Z_matrices)) -> coeffs
  


  # add CI
  coeffs$Lower = coeffs$estimate - 1.96 * coeffs$std_error
  coeffs$Upper = coeffs$estimate + 1.96 * coeffs$std_error


  # theme
  apatheme=ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
          panel.grid.minor=ggplot2::element_blank(),
          panel.border=ggplot2::element_blank(),
          axis.line=ggplot2::element_line(),
          text=ggplot2::element_text(family='Times'),
          legend.title=ggplot2::element_blank())

  # use with to deal with R CMD check (because Trait / Estimate are not defined)
  P= with(coeffs,{ ggplot2::ggplot(data=coeffs, ggplot2::aes(x=Trait, y=estimate,
                     ymin=Lower,
                     ymax=Upper)) +
   # "global estimates"
   ggplot2::geom_pointrange(col="black", shape=21, fill="lightgray")+
   ggplot2::geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip

   #geom_segment(x=3, xend=3, y=3.5, yend=3.7, lty=2, col="#169D74", lwd = 0.45)
   ggplot2::labs(title = "") +
   ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
   ggplot2::xlab("") + ggplot2::ylab("Multivariate MR causal effect estimates (95% CI)") +
   apatheme  # use a white background
  })
 # "coefficient estimates")
 # ggplot order the studies by "Trait"
 coeffs = coeffs[order(coeffs$Trait),]

 for(i in 1:nrow(coeffs)){
   t = coeffs$study[i]
   chrm_estimates = unlist(obj$all_MRcoeffs[obj$all_MRcoeffs$study==t, "estimate"])
   for(c in chrm_estimates){
     P = P +
       ggplot2::geom_segment(x=i-0.1, xend=i+0.1, y=as.numeric(c), yend=as.numeric(c),
                              col="darkgrey", lwd = 0.8)
   }
 }
if(save_file) grDevices::png(file_name, width = 20, height = 12, units = "cm", res = 500)
print(P)
if(save_file) grDevices::dev.off()
}



#' Extract SNPs results from bGWAS results
#'
#' Extracts SNPs results from bGWAS results (BFs, p-value, prior, posterior and 
#' direct effects, depending on the value of the parameter \code{results})
#'
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#' @param SNPs, "all" / "significant", \code{default="significant"}
#' @param results, "BF" / "posterior" / "direct", \code{default="BF"}
#' 
#' @details
#' For all value of \code{results}, basic informations about the SNPs will be returned: \cr
#' \code{rsid} : rs number \cr
#' \code{chrm_UK10K} : chromosome (obtained from UK10K data)\cr
#' \code{pos_UK10K} : position (obtained from UK10K data) \cr
#' \code{alt} : alternative (effect) allele \cr
#' \code{ref} : reference allele \cr
#' \code{z_obs} : observed Z-score \cr
#' 
#' In addition, if \code{results = "BF"} the following information will be returned: \cr
#' \code{mu_prior_estimate} : prior effect estimate (z-score scale) \cr
#' \code{mu_prior_std_error} : prior effect standard error (z-score scale) \cr
#' \code{BF} : Bayes Factor\cr
#' \code{BF_p} : Bayes Factor p-value  \cr
#' \code{BF_fdr} : Bayes Factor FDR (only if FDR used to identify significant SNPs) \cr
#' 
#' Alternatively, if \code{results = "posterior"} the following information will be returned: \cr
#' \code{mu_posterior_estimate} : posterior effect estimate (z-score scale)  \cr
#' \code{mu_posterior_std_error} : posterior effect standard error (z-score scale) \cr
#' \code{z_posterior} : posterior Z-score \cr
#' \code{p_posterior} : posterior effect p-value \cr
#' \code{fdr_posterior} :  posterior effect FDR (only if FDR used to identify significant SNPs) \cr
#' 
#' Alternatively, if \code{results = "direct"} the following information will be returned: \cr
#' \code{mu_direct_estimate} : direct effect estimate (z-score scale)\cr
#' \code{mu_direct_std_error} : direct effect standard error (z-score scale) \cr
#' \code{z_direct} : direct Z-score\cr
#' \code{p_direct} : direct effect p-value \cr
#' \code{fdr_direct} : direct effect FDR (only if FDR used to identify significant SNPs) \cr
#'
#' @return a \code{tibble} containing the results for all / significant SNPs
#' @export

extract_results_bGWAS <- function(obj, SNPs="significant", results="BF"){
  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(!SNPs %in% c("all", "significant")) stop("SNPs : should be \"all\" or \"significant\".")
  if(!results %in% c("BF", "posterior", "direct")) stop("results : should be \"BF\", \"posterior\" or \"direct\".")
  
  if(SNPs=="all"){
    Res = obj$all_BFs
  } else {
    if(results=="BF"){
      if(length(obj$significant_SNPs)==0){
        stop("You can't extract \"significant\" results , there is no significant SNPs", call. = F)
      }
      Res = obj$all_BFs %>% filter(.data$rsid %in% obj$significant_SNPs)
    } else if(results=="posterior"){
      if(length(obj$posterior_SNPs)==0){
        stop("You can't extract \"significant\" results , there is no significant SNPs according to posterior p-value", call. = F)
      }
      Res = obj$all_BFs %>% filter(.data$rsid %in% obj$posterior_SNPs)
    } else if(results=="direct"){
      if(length(obj$direct_SNPs)==0){
        stop("You can't extract \"significant\" results , there is no significant SNPs according to direct p-value", call. = F)
      }
      Res = obj$all_BFs %>% filter(.data$rsid %in% obj$direct_SNPs)
    }
  }
  
  if(results=="BF"){
    Res %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, .data$alt, .data$ref, .data$z_obs,
             .data$mu_prior_estimate, .data$mu_prior_std_error, .data$BF, .data$BF_p, 
             if("BF_fdr" %in% names(.data$.)) .data$BF_fdr else NULL) -> Res
  } else if(results=="posterior"){
    Res %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, .data$alt, .data$ref, .data$z_obs,
             .data$mu_posterior_estimate, .data$mu_posterior_std_error, .data$z_posterior, .data$p_posterior, 
             if("fdr_posterior" %in% names(.data$.)) .data$fdr_posterior else NULL) -> Res
  } else if(results=="direct"){
    Res %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, .data$alt, .data$ref, .data$z_obs,
             .data$mu_direct_estimate, .data$mu_direct_std_error, .data$z_direct, .data$p_direct, 
             if("fdr_direct" %in% names(.data$.)) .data$fdr_direct else NULL) -> Res
  }
  
  
  return(Res)
}




#' Extract MR coefficients from bGWAS results
#'
#' Extracts MR coefficients (multivariate genome-wide and per-chromosome estimates)
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#'
#' @return a \code{tibble} containing the MR coefficients (1 estimate using
#' all chromosomes + 22 estimates with 1 chromosome masked)
#' @export

extract_MRcoeffs_bGWAS <- function(obj){
  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(is.null(obj$all_MRcoeffs)) stop("The prior has not been created using Prior GWASs, there are no coefficients to extract")

  # For each study :
  # global coeff / chromosomes coeff
  Res = obj$significant_studies
  Z_matrices = strsplit(obj$log_info[stringr::str_detect(obj$log_info, "The Z-Matrix files are stored in \"")], "\"")[[1]][2]
  
  Res %>%
    transmute(name = get_names(.data$study, Z_matrices)) -> WithNames
  bind_cols(WithNames, Res) -> Res
  for(ch in 1:22){
    obj$all_MRcoeffs %>%
      filter(.data$chrm == ch) -> CHRM
    Res[,paste0("chrm", ch, "_estimate")] = CHRM$estimate[match(Res$study, CHRM$study)]
    Res[,paste0("chrm", ch, "_std_error")] = CHRM$std_error[match(Res$study, CHRM$study)]
    Res[,paste0("chrm", ch, "_P")] = CHRM$P[match(Res$study, CHRM$study)]
    }
  return(Res)
}



#' Heatmap of SNP effects on prior traits from bGWAS results
#'
#' Creates a heatmap of SNP effects on prior traits
#'
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#' @param save_file A logical indicating if the graphic should be saved,
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device
#' @param file_name The name of the file saved (is \code{save_file} is \code{TRUE})
#'        \code{default=NULL}, will used NameOfYourAnalysis_Heatmap.png
#' @return a Heatmap
#'
#' @importFrom rlang :=
#' @export


heatmap_bGWAS <- function(obj, save_file=F, file_name=NULL) {
  # can only be run on significant SNPs?
  

  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(!is.logical(save_file)) stop("save_file : should be logical")
  # if no name, use the one from the analysis (in log file)
  if(save_file && is.null(file_name)){
    file_name = paste0(strsplit(strsplit(obj$log_info[
      grep("The name of your analysis is: ", obj$log_info)],
      "The name of your analysis is: \"", fixed=T)[[1]][2], "\"", fixed=T)[[1]][1],
      "_Heatmap.png")
  }
  if(save_file && !is.character(file_name)) stop("file_name : should be a character")

  
  if(nrow(obj$significant_studies)<2) stop("heatmap can only be created if at least 2 prior GWASs have been used to create the prior")

  if(length(obj$significant_SNPs)<2)  stop("heatmap can only be created if at least 2 significant hits have been identified")
    
  if(save_file) grDevices::png(file_name, width = 20, height = 12, units = "cm", res = 500)

  Matrix = obj$matrix_heatmap
  Res_signif = extract_results_bGWAS(obj)
  all_causalestimates = obj$all_MRcoeffs

  # 1) align to have obs_effect >0
  # sign = POS if the risk factor has a positive effect on our trait
  #        NEG if the risk factor has a negative effect on our trait
  Res_signif %>%
    mutate(alt = case_when(.data$z_obs<0 ~ .data$ref,
                           TRUE ~ .data$alt),
           z_obs = abs(.data$z_obs),
           z_prior_estimate = abs(.data$mu_prior_estimate)) -> Res_signif_aligned
  
  
  
  # 2) re-align ZMat with GWAS results (using ALT)
  Matrix %>%
    slice(match(Res_signif_aligned$rsid, .data$rs)) -> Matrix
  
  Matrix %>%
    select(-c(1:5))  %>%
    colnames -> RFs
  
  Res_signif_aligned %>%
    mutate( ZMat_alt = Matrix$alt,
            Zmat_ref = Matrix$ref) ->  Res_signif_aligned
  

  for(RF in RFs){
    Res_signif_aligned %>%
      mutate({{RF}} := case_when(
                (alt == .data$ZMat_alt) ~ Matrix[,RF],
                (alt == .data$Zmat_ref) ~ -Matrix[,RF],
                TRUE ~ NA_real_))-> Res_signif_aligned
  }
  
  Res_signif_aligned %>%
    select(RFs) %>%
    as.matrix -> RF_SNPs
  
  RF_contribution = RF_SNPs
  
  Res_signif_aligned %>%
    pull(.data$rsid) -> SNPs
  for(snp in 1:length(SNPs)){
    all_causalestimates %>%
      filter(.data$chrm==as.numeric(Res_signif_aligned[snp, "chrm_UK10K"])) -> causalestimates
    # re-order the RFs before multiplying
    causalestimates = causalestimates[match(colnames(RF_SNPs), causalestimates$study),] 
    RF_contribution[snp,] = RF_SNPs[snp,] * causalestimates$estimate
  }
  
  # check that our sum(RF_contributions) match the prior reported
  # cbind(apply(RF_contribution, 1, function(x) sum(x)) , Res_signif_aligned$z_prior_estimate)
  
  # plot
  significance = RF_SNPs
  limit = stats::qnorm(5e-8*0.5, lower.tail = F)
  significance[abs(significance)>=limit] = "*"
  significance[significance!="*"] = ""
  
  Res_signif_aligned %>%
    tidyr::unite(SNPs_Name, .data$rsid, .data$alt, sep = " - ") %>%
    pull(.data$SNPs_Name) -> SNPs_Name
  
  Z_matrices = strsplit(obj$log_info[stringr::str_detect(obj$log_info, "The Z-Matrix files are stored in \"")], "\"")[[1]][2]
  
  RFs_Names = get_names(RFs, Z_matrices)
  
  my_colors=grDevices::colorRampPalette(c("blue", "white", "red"))(20*20)
  
  
  gplots::heatmap.2(RF_contribution, 
            cellnote = significance,
            notecol="black", 
            notecex = 2,
            Rowv=F, Colv=F,
            dendrogram="none", col=my_colors, trace="none",
            margins=c(10, 10),labCol=RFs_Names, labRow=SNPs_Name , cexCol=0.8, cexRow=0.8,
            lwid=c(0.5, 3.5), keysize=0.75, key.ylab="", 
            key.par = list(cex=0.5, cex.main=0.0001,  cex.axis=1, cex.lab=1, cex.sub=1), srtCol=45,
            density.info='none',
            key.xlab="Contribution to prior effects")
  
  
  if(save_file) grDevices::dev.off()
}



#' Get squared correlation between observed and prior effects from bGWAS results
#'
#' Returns squared correlation between observed and prior effects,
#' for different subsets of SNPs (all, the ones having at least a moderate effects - p-value < 0.001 -,
#' MR instruments) 
#'
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#' @param SNPs, "all" / "moderate" / "instruments"
#' @return a squared correlation
#'
#' @export

get_RSquared_bGWAS <- function(obj, SNPs="all"){ # obj should be a bGWAS object
  # SNPs can be "all", "moderate", "instruments"
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(!SNPs %in% c("all", "moderate", "instruments")) stop("SNPs : should be \"all\", \"moderate\" or \"instruments\".")
  
  R2=NA
  if(SNPs=="all"){
    Line =  obj$log_info[stringr::str_detect(obj$log_info, "Correlation between prior and observed effects for all SNPs is ")]
    R2 = as.numeric(stringr::str_split(Line, "is ")[[1]][2])^2
  } else if(SNPs=="moderate"){
    Line =  obj$log_info[stringr::str_detect(obj$log_info, "Correlation between prior and observed effects for SNPs with GWAS p-value < 0.001 is ")]
    R2 = as.numeric(stringr::str_split(Line, "is ")[[1]][2])^2
  } else if(SNPs=="instruments"){
    Line =  obj$log_info[stringr::str_detect(obj$log_info, "Out-of-sample squared correlation for MR instruments across all chromosome is ")]
    R2 = as.numeric(stringr::str_split(Line, "is ")[[1]][2])
  }
  
  return(R2)
}



#' Print log from bGWAS results
#'
#' Prints the log (everything that is printed during a \code{bGWAS} analysis)
#' with \code{verbose=TRUE})
#'
#'
#' @param obj an object of class bGWAS created using \code{\link{bGWAS}()}
#'
#' @export

print_log_bGWAS <- function(obj){ # obj should be a bGWAS object
  # SNPs can be "all", "moderate", "instruments"
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")

  
  cat(obj$log_info, sep = "\n")
}


### Useful small function (re-used by other functions)

update_log <- function(log_obj, text, verbose=F){
  log_obj = c(log_obj, text)
  if(verbose) cat(text, sep = "")
  return(log_obj)
}







