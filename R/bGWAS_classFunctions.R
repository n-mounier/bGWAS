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
  
  if(any(stringr::str_detect(x$log_info, "Analysis failed"))){
     cat("Analysis failed")
  }else{
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
    
    if(any(is.na(x$significant_SNPs))){
      cat ("No significant SNP identified.")
    } else if(length(x$significant_SNPs)==0){
      cat("No significant SNP identified, because the analysis has been limited to prior estimation")
    } else
      if(length(x$significant_SNPs)==1){
        cat(length(x$significant_SNPs), "significant SNP identified : \n ")
        cat(x$significant_SNPs, sep=", ")
      } else {
        cat(length(x$significant_SNPs), "significant SNPs identified : \n ")
        cat(x$significant_SNPs, sep=", ")
      }
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
#'        be aware that \code{annotate=TRUE} might decrease readability of the figure. You
#'        could define a set of SNPs to annotate using \code{SNPs}.
#' @param SNPs A data.frame containing the SNPs (rsid) to annotate in the first column, and 
#'        optionnally the text that should be plotted in the second column, and the color
#'        in the third column, \code{default=NULL}, only evaluated if \code{annotate=TRUE}.
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
                                  annotate=T, SNPs = NULL, results="BF") {

  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")

  
  if(any(stringr::str_detect(obj$log_info, "Analysis failed"))){
    stop("Manhattan plot can not be displayed, because the analysis failed", call. = F)
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
  if(annotate && !(is.data.frame(SNPs) | is.null(SNPs))) stop("SNPs : should be a data.frame")
  if(is.data.frame(SNPs) && ncol(SNPs)>3)  stop("SNPs : should not have more than 2 columns")
  
  
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
  
  
  # for SNPs with super low p, shrink it and add a red point
  ToPlot %>%
    filter(.data$p<1e-20) -> SuperSignif
  
  ToPlot %>% 
    mutate(p = case_when(
      .data$p<10e-20 ~ 1e-20,
      TRUE ~ .data$p)) -> ToPlot


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
                   ylab= my_ylab,
                   ylim=c(1.3, ifelse(nrow(SuperSignif)>1, 21, -log10(min(ToPlot$p))+1)),
                   cex.lab=1.2,
                   cex.axis=0.8,
                   yaxp  = c(0, 
                             ifelse(nrow(SuperSignif)>1, 20, ceiling(-log10(min(ToPlot$p)))), 
                             ifelse(nrow(SuperSignif)>1, 4, 5)))
  graphics::abline(h=1.3, lwd=2)
  
  if(nrow(SuperSignif)>0){
    ToPlot %>%
      mutate(p=NULL) %>%
      arrange(.data$chrm_UK10K,.data$pos_UK10K) -> all
    
    # y = -log10(1e20)
    SuperSignif %>%
      mutate(y = -log10(1e-20)) -> SuperSignif 
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
    
    SuperSignif$x = apply(SuperSignif, 1, function(x) get_posx(x, all))
    
    graphics::points(x=SuperSignif$x, y=SuperSignif$y, col="blue")
  }



  if(annotate){ # significant SNPs from the analysis
    use_color=F
    if(is.data.frame(SNPs)){ 
      SNPs %>%
        stats::setNames( c("rs", "name", "color")[1:ncol(SNPs)]) -> SNPs
      SNPs %>%
        pull(1) -> rsids
      ToPlot %>%
        slice(match(rsids,.data$rsid)) -> SNPs_to_plot
      if(ncol(SNPs)>1){
        SNPs %>%
          slice(match(SNPs_to_plot$rsid, .data$rs)) %>%
          pull(2) -> labels
      } else {
        SNPs_to_plot %>% pull(.data$rsid) -> labels
      }
      if(ncol(SNPs)>2){
        use_color = T
        SNPs %>%
          slice(match(SNPs_to_plot$rsid, .data$rs)) %>%
          pull(3) -> my_colors
      }
    } else {
    # extract them
      if(results=="BF"){
        ToPlot %>%
          filter( .data$rsid %in% obj$significant_SNPs) -> SNPs_to_plot
        SNPs_to_plot %>% pull(.data$rsid) -> labels
      } else if(results=="posterior"){
        ToPlot %>%
          filter( .data$rsid %in% obj$posterior_SNPs) -> SNPs_to_plot
        SNPs_to_plot %>% pull(.data$rsid) -> labels
      } else if(results=="direct"){
        ToPlot %>%
          filter( .data$rsid %in% obj$direct_SNPs) -> SNPs_to_plot
        SNPs_to_plot %>% pull(.data$rsid) -> labels
      } 
      
    }
    
    ToPlot %>%
      mutate(p=NULL) %>%
      arrange(.data$chrm_UK10K,.data$pos_UK10K) -> all

    # y = -log10(p) ou -log10(fdr)
    SNPs_to_plot %>%
      mutate(y = -log10(.data$p)) -> SNPs_to_plot # to make sure all SNPs names are in plotting windows
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
    
    if(use_color){
      graphics::points(x=SNPs_to_plot$x, y=SNPs_to_plot$y, col=my_colors, cex=.7)
      
      calibrate::textxy(SNPs_to_plot$x, SNPs_to_plot$y+0.5,
                        labs= labels,
                        col = my_colors,
                        offset = 0.45, cex=0.6)
    } else {
      calibrate::textxy(SNPs_to_plot$x, SNPs_to_plot$y+0.5,
                        labs= labels,
                        offset = 0.45, cex=0.6)    }
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
  
  if(any(stringr::str_detect(obj$log_info, "Analysis failed"))){
    stop("Coefficients plot can not be displayed, because the analysis failed", call. = F)
  }
  
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
  P= ggplot2::ggplot(data=coeffs, ggplot2::aes(x=.data$Trait, y=.data$estimate,
                     ymin=.data$Lower,
                     ymax=.data$Upper)) +

   ggplot2::geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip

   ggplot2::labs(title = "") +
   ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
   ggplot2::xlab("") + ggplot2::ylab("Multivariate MR causal effect estimates (95% CI)") +
   apatheme  # use a white background
  
 # "coefficient estimates"
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
 
 P = P +
   # "global estimates" : add them at end to make sure they are in front
   ggplot2::geom_pointrange(col="black", shape=21, fill="black")
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
#' @param results, "BF" / "posterior" / "direct" / "everything", \code{default="BF"}
#' 
#' @details
#' For all value of \code{results}, basic informations about the SNPs will be returned: \cr
#' \code{rsid} : rs number \cr
#' \code{chrm_UK10K} : chromosome (obtained from UK10K data)\cr
#' \code{pos_UK10K} : position (obtained from UK10K data) \cr
#' \code{alt} : alternative (effect) allele \cr
#' \code{ref} : reference allele \cr
#' \code{beta} : observed effect size (if possible) \cr
#' \code{se} : observed effect size (if possible) \cr
#' \code{z_obs} : observed Z-score \cr
#' 
#' In addition, if \code{results = "BF"} the following information will be returned: \cr
#' \code{mu_prior_estimate} : prior effect estimate (z-score scale) \cr
#' \code{mu_prior_std_error} : prior effect standard error (z-score scale) \cr
#' \code{beta_prior_estimate} : prior effect estimate (beta scale, if possible) \cr
#' \code{beta_prior_std_error} : prior effect standard error (beta scale, if possible) \cr
#' \code{BF} : Bayes Factor\cr
#' \code{BF_p} : Bayes Factor p-value  \cr
#' \code{BF_fdr} : Bayes Factor FDR (only if FDR used to identify significant SNPs) \cr
#' 
#' Alternatively, if \code{results = "posterior"} the following information will be returned: \cr
#' \code{mu_posterior_estimate} : posterior effect estimate (z-score scale)  \cr
#' \code{mu_posterior_std_error} : posterior effect standard error (z-score scale) \cr
#' \code{beta_posterior_estimate} : posterior effect estimate (beta scale, if possible) \cr
#' \code{beta_posterior_std_error} : posterior effect standard error (beta scale, if possible) \cr
#' \code{z_posterior} : posterior Z-score \cr
#' \code{p_posterior} : posterior effect p-value \cr
#' \code{fdr_posterior} :  posterior effect FDR (only if FDR used to identify significant SNPs) \cr
#' 
#' Alternatively, if \code{results = "direct"} the following information will be returned: \cr
#' \code{mu_direct_estimate} : direct effect estimate (z-score scale)\cr
#' \code{mu_direct_std_error} : direct effect standard error (z-score scale) \cr
#' \code{beta_direct_estimate} : direct effect estimate (beta scale, if possible) \cr
#' \code{beta_direct_std_error} : direct effect standard error (beta scale, if possible) \cr
#' \code{z_direct} : direct Z-score\cr
#' \code{p_direct} : direct effect p-value \cr
#' \code{fdr_direct} : direct effect FDR (only if FDR used to identify significant SNPs) \cr
#' \code{CRR} : corrected to raw ratio (ratio between direct effect and observed effect) \cr
#' 
#' Alternatively, if \code{results = "everything"} all the results described above will be returned 
#' (possible only if \code{results . \cr
#'
#' @return a \code{tibble} containing the results for all / significant SNPs
#' @export

extract_results_bGWAS <- function(obj, SNPs="significant", results="BF"){
  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(any(stringr::str_detect(obj$log_info, "Analysis failed"))){
    stop("Results can not be extracted, because the analysis failed", call. = F)
  }
  
  
  if(!SNPs %in% c("all", "significant")) stop("SNPs : should be \"all\" or \"significant\".")
  if(!results %in% c("BF", "posterior", "direct", "everything")) stop("results : should be \"BF\", \"posterior\", \"direct\" or \"everything\".")
  if(SNPs=="significant" && results == "everything") stop("Not possible to return \"everything\" only for \"significant\" SNPs.")
  
  
  
  if(SNPs=="all" && results == "everything"){
    return(obj$all_BFs)
  }
  
  
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
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, .data$alt, .data$ref, 
             ends_with("beta"), ends_with("se"),
             .data$z_obs,
             .data$mu_prior_estimate, .data$mu_prior_std_error, 
             matches("beta_prior_estimate"),  matches("beta_prior_std_error"),
             .data$BF, .data$BF_p, 
             matches("BF_fdr")) -> Res
  } else if(results=="posterior"){
    Res %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, .data$alt, .data$ref,
             ends_with("beta"), ends_with("se"),
             .data$z_obs,
             .data$mu_posterior_estimate, .data$mu_posterior_std_error, 
             matches("beta_posterior_estimate"),  matches("beta_posterior_std_error"),
             .data$z_posterior, .data$p_posterior, 
             matches("fdr_posterior")) -> Res
  } else if(results=="direct"){
    Res %>%
      select(.data$rsid, .data$chrm_UK10K, .data$pos_UK10K, .data$alt, .data$ref,
             ends_with("beta"), ends_with("se"),
             .data$z_obs,
             .data$mu_direct_estimate, .data$mu_direct_std_error, 
             matches("beta_direct_estimate"),  matches("beta_direct_std_error"),
             .data$z_direct, .data$p_direct, 
             matches("fdr_direct"),
             matches("CRR")) -> Res
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
  
  if(any(stringr::str_detect(obj$log_info, "Analysis failed"))){
    stop("MR coefficients can not be extracted, because the analysis failed", call. = F)
  }
  
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
#' @param SNPs A data.frame containing the SNPs (rsid) to use in the first column, and 
#'        optionnally the text that should be plotted in addition to rsid in the second 
#'        column \code{default=NULL}.
#' @return a Heatmap
#'
#' @importFrom rlang :=
#' @export


heatmap_bGWAS <- function(obj, SNPs=NULL, save_file=F, file_name=NULL) {
  # can only be run on significant SNPs?
  

  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(any(stringr::str_detect(obj$log_info, "Analysis failed"))){
    stop("Heatmap can not be displayed, because the analysis failed", call. = F)
  }
  
  if(any(is.na(obj$significant_SNPs))) stop("Heatmap can not be displayed, because there is no significant SNP", call. = F) 
  
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
    
  if(!is.null(SNPs) && !is.data.frame(SNPs)) stop("SNPs : should be a data.frame")
  if(!is.null(SNPs)){
    SNPs %>%
      stats::setNames( c("rs", "name")[1:ncol(SNPs)]) -> SNPs
  }
  if(is.data.frame(SNPs) && !all(SNPs$rs %in% extract_results_bGWAS(obj)$rsid)) stop("SNPs : the rsids provided do not match the ones of significant SNPs")
  if(is.data.frame(SNPs) && ncol(SNPs)>2) stop("SNPs : should not have more than 2 columns")
  
  if(save_file) grDevices::png(file_name, width = 20, height = 12, units = "cm", res = 500)

  Matrix = obj$matrix_heatmap
  Res_signif = extract_results_bGWAS(obj)
  if(!is.null(SNPs)){
    SNPs %>%
      stats::setNames( c("rs", "name")[1:ncol(SNPs)]) -> SNPs
    Res_signif %>%
      slice(match(SNPs$rs, .data$rsid)) -> Res_signif
  }
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
                (.data$alt == .data$ZMat_alt) ~ pull(Matrix, RF),
                (.data$alt == .data$Zmat_ref) ~ -pull(Matrix,RF),
                TRUE ~ NA_real_))-> Res_signif_aligned
  }
  
  Res_signif_aligned %>%
    select(RFs) %>%
    as.matrix -> RF_SNPs
  
  RF_contribution = RF_SNPs
  
  Res_signif_aligned %>%
    pull(.data$rsid) -> all_SNPs
  for(snp in 1:length(all_SNPs)){
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
  
  if(!is.null(SNPs) && ncol(SNPs)>1){
    SNPs_Name = paste0(SNPs_Name, " (", SNPs$name, ")")
  }
  
  Z_matrices = strsplit(obj$log_info[stringr::str_detect(obj$log_info, "The Z-Matrix files are stored in \"")], "\"")[[1]][2]
  
  RFs_Names = get_names(RFs, Z_matrices)
  
  my_colors=grDevices::colorRampPalette(c("blue", "white", "red"))(20*20)
  
  
  gplots::heatmap.2(RF_contribution, 
            cellnote = significance,
            notecol="black", 
            notecex = 2,
            Rowv=F, Colv=F,
            dendrogram="none", col=my_colors, trace="none",
            margins=c(10, 10),labCol=RFs_Names, labRow=SNPs_Name , cexCol=0.8, cexRow=0.7,
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
  if(any(stringr::str_detect(obj$log_info, "Analysis failed"))) return(NA)
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







