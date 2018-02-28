###### bGWAS class functions ######


# SignificantSNPs = data.table::fread("../Tests/Test_UsingSmallDataFrame/SignificantSNPs.csv")
# AllResults = data.table::fread("../Tests/Test_UsingSmallDataFrame/PriorBFp.csv")
# SignificantStudies = data.table::fread("../Tests/Test_UsingSmallDataFrame/Multivariate_coefs.csv")
# AllCoeffs = data.table::fread("../Tests/Test_UsingSmallDataFrame/CoefficientsByChromosome.csv")
#
# MyObj = list()
#
# MyObj$significant_SNPs = SignificantSNPs$rs
# MyObj$all_BFs = AllResults
# MyObj$significant_studies = SignificantStudies
# MyObj$all_MRcoeffs = AllCoeffs
# MyObj$log_info = read.csv("../Tests/Test_UsingSmallDataFrame.log")
# colnames(MyObj$all_MRcoeffs)=c("Chrm", "StudyName", "Estimate", "StdError", "TValue", "P")
# colnames(MyObj$significant_studies)=c("StudyName", "Estimate", "StdError", "TValue", "P")
# class(MyObj)="bGWAS"



#' @param obj an object of class bGWAS
#'
#' @return print
#' @export
print.bGWAS <- function(obj) {
  cat("-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n ")

  cat( "bGWAS performed on", format( nrow(obj$all_BFs) , big.mark = "," , scientific = F) ,
                                    "SNPs \n \n" )

  cat("-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n")

  cat(nrow(obj$significant_studies), "studies used to build the prior : \n")
  print(obj$significant_studies[,1:3], row.names=F)

  cat("\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n")

  if(length(obj$significant_SNPs)==1){
    cat(length(obj$significant_SNPs), "significant SNP identified : \n ")
  } else {
    cat(length(obj$significant_SNPs), "significant SNPs identified : \n ")
  }
  cat(obj$significant_SNPs, sep=", ")
  # independant ? which threshold ? get info in log
  cat("\n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ ")

  }

# print(MyObj)


###### Function to create pruned Z-matrix of MR instruments ######



#' Manhattan Plot from bGWAS results
#'
#' Create a Manhattan Plot from bGWAS results (object of class bGWAS obtained
#' when using bGWAS() or bGWASFromprior()),
#'
#'
#' @param obj an object of class bGWAS
#' @param save_file A logical indicating if the graphic should be saved,
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device
#' @param file_name The name of the file saved (is \code{save_file} is \code{TRUE})
#'        \code{default=NULL}, will used NameOfYourAnalysis_ManhattanPlot.png
#' @param threshold The threshold used to draw the significance line,
#'        \code{default=NULL}, will used the threshold used in the analysis
#' @param annotate A logical indicating if the significant SNPs identified in the
#'        analysis should be annotated on the plot, \code{default=TRUE}
#'        If your results are not pruned or if you have a high number of significant SNPs,
#'        \code{annotate=TRUE} might decrease readability of the figure.
#'
#' @return a Manhattan Plot
#' @export

manatthan_plot_bGWAS <- function(obj, save_file=F, file_name=NULL,
                                 threshold=NULL, annotate=T) {

  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(!is.logical(save_file)) stop("save_file : should be logical")
  # if no name, use the one from the analysis (in log file)
  if(save_file && is.null(file_name)){
    file_name = paste0(strsplit(strsplit(obj$log_info[
      grep("The name of your analysis is: ", obj$log_info)],
      "The name of your analysis is: \"", fixed=T)[[1]][2], "\"", fixed=T)[[1]][1],
    "_ManatthanPlot.png")
  }
  if(save_file && !is.character(file_name)) stop("file_name : should be a character")
  if(!is.logical(annotate)) stop("annotate : shoud be logical")
  if(!is.null(threshold) && !is.numeric(threshold)) stop("threshold : shoud be numeric")
  # if no threshold, use the one from analysis (in log file)
  if(is.null(threshold)){
    method = strsplit(strsplit(obj$log_info[
    grep("Significant SNPs will be identified according to ", obj$log_info)],
    "Significant SNPs will be identified according to ")[[1]][2], ".", fixed=T)[[1]][1]
    threshold = as.numeric(strsplit(strsplit(obj$log_info[
      grep("Significant SNPs will be identified according to ", obj$log_info)],
      "The threshold used is :")[[1]][2], ".", fixed=T)[[1]][1])
  } else {
      method="p"
  }



  if(save_file) png(file_name, width = 20, height = 12, units = "cm", res = 500)

  qqman::manhattan(obj$all_BFs ,
                   chr = "chrm" ,
                   bp = "pos" ,
                   p = ifelse( method == "FDR" ,
                               "fdr" ,
                               "BF_p" ) ,
                   snp = "rs",
                   col = c("gray10", "gray60"),
                   suggestiveline = FALSE,
                   genomewideline = -log10(threshold),
                   highlight = NULL,
                   logp = TRUE,
                   annotatePval = NULL,
                   annotateTop = FALSE,
                   ylab=ifelse( method == "FDR" ,
                                expression(-log[10](italic(fdr))) ,
                                expression(-log[10](italic(p))) )
                   )



  if(annotate){ # significant SNPs from the analysis
    # extract them
    value = ifelse(method=="FDR", "fdr", "BF_p")
    SNPs_to_plot = obj$all_BFs[rs %in% obj$significant_SNPs,
                               c("rs", "chrm", "pos", value), with=F]
    all = obj$all_BFs[,c("rs", "chrm", "pos")]
    all <- all[order(all$chrm, all$pos), ]

    # y = -log10(p) ou -log10(fdr)
    SNPs_to_plot$y = -log10(SNPs_to_plot[,..value, with=F])

     # x = have to look at all SNPs chr/pos
    get_posx <- function(snp, all){
      chr = as.numeric(snp[2])
      pos = as.numeric(snp[3])
      if(chr==1){
        posx = pos
      }
      else{
        posx=0
        for(c in 1:(chr-1)){
          posx = posx + tail(subset(all, chrm==c)$pos,1)
        }
        posx = posx + pos
      }
      return(posx)
    }

    SNPs_to_plot$x = apply(SNPs_to_plot, 1, function(x) get_posx(x, all))

    calibrate::textxy(SNPs_to_plot$x, SNPs_to_plot$y,
                      labs=SNPs_to_plot$rs,
                      offset = 0.625, cex=0.45)
  }


  if(save_file) dev.off()
}


# manatthan_plot_bGWAS(MyObj, save_file=T)


#' Ceofficients Plot from bGWAS results
#'
#' Create a Coefficients Plot (causal effect of Prior GWASs) from bGWAS results (object of class bGWAS obtained
#' when using bGWAS())
#'
#'
#' @param obj an object of class bGWAS
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

  cat("\n \n")
  cat(file_name)
  cat("\n \n")

  coeffs = obj$significant_studies
  # add the trait name (if multiple studies for a same trait, add " - X" )
  all_studies =  list_priorGWASs()
  coeffs$Trait = all_studies[match(coeffs$Study, all_studies$File), Trait]
  if(length(unique(coeffs$Trait))!=nrow(coeffs)){
    for(t in unique(coeffs$Trait)){
      if(nrow(coeffs[coeffs$Trait==t]) >1 )
      coeffs$Trait[coeffs$Trait==t] = paste0(t, " (", c(1:sum(coeffs[,coeffs$Trait==t])), ")")
    }
  }


  # add CI
  coeffs$Lower = coeffs$Estimate - 1.96 * coeffs$StdError
  coeffs$Upper = coeffs$Estimate + 1.96 * coeffs$StdError


  # theme
  apatheme=ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
          panel.grid.minor=ggplot2::element_blank(),
          panel.border=ggplot2::element_blank(),
          axis.line=ggplot2::element_line(),
          text=ggplot2::element_text(family='Times'),
          legend.title=ggplot2::element_blank())

 P=  ggplot2::ggplot(data=coeffs, ggplot2::aes(x=Trait, y=Estimate,
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
 # "coefficient estimates")
 # ggplot order the studies by "Trait"
 coeffs = coeffs[order(coeffs$Trait),]
 for(i in 1:nrow(coeffs)){
   t = coeffs$Study[i]
   chrm_estimates = unlist(obj$all_MRcoeffs[obj$all_MRcoeffs$Study==t, "Estimate"])
   for(c in chrm_estimates){
     P = P +
       ggplot2::geom_segment(x=i-0.1, xend=i+0.1, y=as.numeric(c), yend=as.numeric(c),
                              col="darkgrey", lwd = 0.8)
   }
 }
if(save_file) png(file_name, width = 20, height = 12, units = "cm", res = 500)
print(P)
if(save_file) dev.off()
}

# coefficients_plot_bGWAS(MyObj, save_file = T)


#' Exctract Results
#'
#' Extract results from an object of class bGWAS obtained
#' when using bGWAS() or bGWAS_fromPrior()
#'
#'
#' @param obj an object of class bGWAS
#' @param save_file A logical indicating if the graphic should be saved,
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device
#' @param file_name The name of the file saved (is \code{save_file} is \code{TRUE})
#'        \code{default=NULL}, will used NameOfYourAnalysis_CoefficientsPlot.png

#' @return a Coefficients Plot
#' @export

extract_results_bGWAS <- function(obj, save_file=F, file_name=NULL){
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

  cat("\n \n")
  cat(file_name)
  cat("\n \n")

  coeffs = obj$significant_studies
  # add the trait name (if multiple studies for a same trait, add " - X" )
  all_studies =  list_priorGWASs()
  coeffs$Trait = all_studies[match(coeffs$Study, all_studies$File), Trait]
  if(length(unique(coeffs$Trait))!=nrow(coeffs)){
    for(t in unique(coeffs$Trait)){
      if(nrow(coeffs[coeffs$Trait==t]) >1 )
        coeffs$Trait[coeffs$Trait==t] = paste0(t, " (", c(1:sum(coeffs[,coeffs$Trait==t])), ")")
    }
  }


  # add CI
  coeffs$Lower = coeffs$Estimate - 1.96 * coeffs$StdError
  coeffs$Upper = coeffs$Estimate + 1.96 * coeffs$StdError


  # theme
  apatheme=ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   panel.border=ggplot2::element_blank(),
                   axis.line=ggplot2::element_line(),
                   text=ggplot2::element_text(family='Times'),
                   legend.title=ggplot2::element_blank())

  P=  ggplot2::ggplot(data=coeffs, ggplot2::aes(x=Trait, y=Estimate,
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
  # "coefficient estimates")
  # ggplot order the studies by "Trait"
  coeffs = coeffs[order(coeffs$Trait),]
  for(i in 1:nrow(coeffs)){
    t = coeffs$Study[i]
    chrm_estimates = unlist(obj$all_MRcoeffs[obj$all_MRcoeffs$Study==t, "Estimate"])
    for(c in chrm_estimates){
      P = P +
        ggplot2::geom_segment(x=i-0.1, xend=i+0.1, y=as.numeric(c), yend=as.numeric(c),
                              col="darkgrey", lwd = 0.8)
    }
  }
  if(save_file) png(file_name, width = 20, height = 12, units = "cm", res = 500)
  print(P)
  if(save_file) dev.off()
}





### Useful small functions (re-used by different functions)

update_log <- function(log_obj, text, verbose=F){
  log_obj = c(log_obj, text)
  if(verbose) cat(text)
  return(log_obj)
}
