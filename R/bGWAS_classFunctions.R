###### bGWAS class functions ######



#' @param obj an object of class bGWAS
#'
#' @return print
#' @export
print.bGWAS <- function(obj) {
  cat("-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n ")
  cat (paste0(" Analysis : \"",strsplit(strsplit(obj$log_info[
    grep("The name of your analysis is: ", obj$log_info)],
    "The name of your analysis is: \"", fixed=T)[[1]][2], "\"", fixed=T)[[1]][1],"\" \n"))
  cat( "bGWAS performed on", format( nrow(obj$all_BFs) , big.mark = "," , scientific = F) ,
                                    "SNPs \n \n" )

  cat("-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n")

  cat(nrow(obj$significant_studies), "studies used to build the prior : \n")
  print(obj$significant_studies[,1:3], row.names=F)

  cat("\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ \n \n")

  if(length(obj$significant_SNPs)==0){
    cat("No significant SNP identified, because the analysis has been limited to prior estimation")
  } else
    if(length(obj$significant_SNPs)==1){
      cat(length(obj$significant_SNPs), "significant SNP identified : \n ")
      cat(obj$significant_SNPs, sep=", ")
    } else {
      cat(length(obj$significant_SNPs), "significant SNPs identified : \n ")
      cat(obj$significant_SNPs, sep=", ")
    }
  # independant ? which threshold ? get info in log
  cat("\n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ ")

  }

# print(MyObj)



#' @param obj1 an object of class bGWAS
#' @param obj2 an object of class bGWA
#'
#' @return all.equal
#' @export
all.equal.bGWAS <- function(obj1, obj2) {
  if(length(obj1)!= length(obj2)) return(FALSE)
  # log - we don't care
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
#' Create a Manhattan Plot from bGWAS results (object of class bGWAS obtained
#' when using bGWAS() or bGWASfromPrior()),
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
#'
#' @export

manhattan_plot_bGWAS <- function(obj, save_file=F, file_name=NULL,
                                 threshold=NULL, annotate=T) {

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
      method="p-value"
  }



  if(save_file) grDevices::png(file_name, width = 20, height = 12, units = "cm", res = 500)

  qqman::manhattan(obj$all_BFs ,
                   chr = "chrm" ,
                   bp = "pos" ,
                   p = ifelse( method == "FDR" ,
                               "fdr" ,
                               "BF_P" ) ,
                   snp = "rs",
                   col = c("gray10", "gray60"),
                   cex = 0.6,
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
    value = ifelse(method=="FDR", "fdr", "BF_P")
    SNPs_to_plot = subset(obj$all_BFs, rs %in% obj$significant_SNPs,
                               c("rs", "chrm", "pos", value))
    all = obj$all_BFs[,c("rs", "chrm", "pos")]
    all <- all[order(all$chrm, all$pos), ]

    # y = -log10(p) ou -log10(fdr)
    SNPs_to_plot$y = -log10(SNPs_to_plot[,value]) - 0.3 # to make sure all SNPs names are in plotting windows
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


  if(save_file) grDevices::dev.off()
}




#' Coefficients Plot from bGWAS results
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


  coeffs = obj$significant_studies
  # add the trait name (if multiple studies for a same trait, add " - X" )
  Z_matrices = strsplit(obj$log_info[stringr::str_detect(obj$log_info, "The Z-Matrix files are stored in \"")], "\"")[[1]][2]
  all_studies =  list_priorGWASs(Z_matrices = Z_matrices)
  coeffs$Trait = unlist(all_studies[match(coeffs$study, all_studies$File), "Trait"])
  if(length(unique(coeffs$Trait))!=nrow(coeffs)){
    for(t in unique(coeffs$Trait)){
      if(nrow(coeffs[coeffs$Trait==t,]) >1 )
      coeffs$Trait[coeffs$Trait==t] = paste0(t, " (", c(1:nrow(coeffs[coeffs$Trait==t,])), ")")
    }
  }


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

# coefficients_plot_bGWAS(MyObj, save_file = T)


#' Extract SNPs results from bGWAS results
#'
#' Extract SNPs results from an object of class bGWAS obtained
#' when using bGWAS() or bGWAS_fromPrior()
#'
#'
#' @param obj an object of class bGWAS
#' @param SNPs, "all" / "significant"
#'
#' @return a data.frame containing the results for all / significant SNPs
#' @export

extract_results_bGWAS <- function(obj, SNPs="significant"){
  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(!SNPs %in% c("all", "significant")) stop("SNPs : should be \"all\" or \"significant\".")

  if(SNPs=="all"){
    Res = obj$all_BFs
  } else {
    if(length(obj$significant_SNPs)==0){
      stop("You can't extract \"significant\" results , because the analysis has been limited to prior estimation", call. = F)
    }
    Res = subset(obj$all_BFs, rs %in% obj$significant_SNPs)
  }

  return(Res)
}




#' Extract MR coefficients from bGWAS results
#'
#' Extract MR coefficients from an object of class bGWAS obtained
#' when using bGWAS()
#'
#'
#' @param obj an object of class bGWAS created using `bGWAS()`
#'
#' @return a data.frame containing the MR coefficients (1 estimate using
#' all chromosomes + 22 estimates with 1 chromosome masked)
#' @export

extract_MRcoeffs_bGWAS <- function(obj){
  ## check parameters
  if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
  if(is.null(obj$all_MRcoeffs)) stop("The prior has not been created using Prior GWASs, there are no coefficients to extract")

  # For each study :
  # global coeff / chromosomes coeff
  Res=data.frame(obj$significant_studies)
  for(c in 1:22){
    CHRM = subset(obj$all_MRcoeffs, chrm==c)
    Res[,paste0("chrm", c, "_estimate")] = CHRM$estimate[match(Res$study, CHRM$study)]
    Res[,paste0("chrm", c, "_std_error")] = CHRM$std_error[match(Res$study, CHRM$study)]
    Res[,paste0("chrm", c, "_P")] = CHRM$P[match(Res$study, CHRM$study)]
    }
  return(Res)
}



#' Heatmap of SNP effects on prior traits from bGWAS results
#'
#' Create a heatmap of SNP effects on prior traits from bGWAS results
#' (object of class bGWAS obtained when using bGWAS() or bGWASfromPrior()),
#'
#'
#' @param obj an object of class bGWAS
#' @param save_file A logical indicating if the graphic should be saved,
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device
#' @param file_name The name of the file saved (is \code{save_file} is \code{TRUE})
#'        \code{default=NULL}, will used NameOfYourAnalysis_Heatmap.png
#' @return a Heatmap
#'
#' @export





heatmap_bGWAS <- function(obj, save_file=F, file_name=NULL) {


  col=grDevices::colorRampPalette(c("blue", "white", "red"))(16*13)
  ## ALIGNEMENT !!!
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


  if(save_file) grDevices::png(file_name, width = 20, height = 12, units = "cm", res = 500)

  Effects = obj$nonZero_effects
  SNPs = obj$significant_SNPs
  # + keep only top hits !
  Effects = Effects[match(obj$significant_SNPs, Effects$rs),]
  # what if a SNP is a hit but has no significant effect from any risk factor ? Add a line with all 0 ?


  # sign = POS if the risk factor has a positive effect on our trait
  #        NEG if the risk factor has a negative effect on our trait
  Effects_Aligned = Effects[,6:(ncol(obj$nonZero_effects)-1)]

  Sign = ifelse(obj$significant_studies$estimate[match(colnames(Effects_Aligned), obj$significant_studies$study)]>0, "POS", "NEG")
  for(t in 1:ncol(Effects_Aligned)){
    if(Sign[t]=="NEG"){
      tName = colnames(Effects_Aligned)[t]
      Effects_Aligned[, (tName) := -(obj$nonZero_effects[obj$nonZero_effects$rs %in% obj$significant_SNPs,..tName, with=F])]
    }
  }
  # keep 0 = 0 or set to NA ?
  # row distance + clustering
  rd<-dist(Effects_Aligned)
  rc<-hclust(rd)
  # column distance + clustering
  cd<-dist(t(Effects_Aligned))
  cc<-hclust(cd)


  # Add Order ???

  # Heatmap : beta, aligned risk scores
  col=grDevices::colorRampPalette(c("blue", "white", "red"))(16*13)
  SNPsAlleles = paste(Effects$rs, Effects$alt, sep=" - ")

  gplots::heatmap.2(as.matrix(Effects_Aligned), Rowv=as.dendrogram(rc), Colv=as.dendrogram(cc),
            dendrogram="both", col=col, trace="none",
            margins=c(9.5, 7.7),labCol=colnames(Effects_Aligned), labRow=SNPsAlleles , cexCol=0.7, cexRow=0.6,
            #lhei=c(,6),
            lwid=c(1.5, 3.5), keysize=0.75, key.ylab="Count", key.par = list(cex=0.5, cex.main=1, cex.axis=1, cex.lab=1, cex.sub=1), srtCol=45,
            key.xlab="Z-Score")


 # gplots::heatmap.2(as.matrix(Effects_Aligned), Rowv=as.dendrogram(rc), Colv=as.dendrogram(cc),#, reorder(dendo, Order),
 #           dendrogram="col", col=col, trace="none",
 #           margins=c(12, 0),labCol=colnames(Effects_Aligned), labRow=rep("", nrow(Effects_Aligned)) , cexCol=0.7, cex=0.5,
            #lhei=c(2,4),
 #           lwid=c(1.5,3.5), keysize=0.75, key.par = list(cex=0.5)) #keysize=0.9)
 # heatmap(as.matrix(Effects_Aligned), Rowv=as.dendrogram(rc), Colv=as.dendrogram(cc),
 #         col = col, margins=c(14,0),labCol=colnames(Effects_Aligned), labRow=rep("", nrow(Effects_Aligned)))



  if(save_file) grDevices::dev.off()
}


### Useful small function (re-used by other functions)

update_log <- function(log_obj, text, verbose=F){
  log_obj = c(log_obj, text)
  if(verbose) cat(text)
  return(log_obj)
}






# miami_plot_bGWAS() <- function (obj, P1, P2, threshold_min){
#   if(class(obj) != "bGWAS") stop("Function implemented for objets of class \"bGWAS\" only.")
#   if(is.null(obj$all_MRcoeffs)) stop("The prior has not been created using Prior GWASs, there are no coefficients to extract")
#
#   # add corrected ? "Difference"
#   if(!P1 %in% c("Obs", "Prior", "Posterior", "Difference")) stop("P1: should be \"Obs\", \"Prior\", \"Posterior\" or \"Difference\".")
#   if(!P2 %in% c("Obs", "Prior", "Posterior", "Difference")) stop("P2: should be \"Obs\", \"Prior\", \"Posterior\" or \"Difference\".")
#   if(P1==P2) stop("P1 and P2 should be different!.")
#
#   if(!is.numeric(threshold)) stop("threshold: should be numeric.")
#
#   # 1) Identify columns and get P-values
#   myData = obj$all_BFs[,1:3]
#   if(P1=="Obs"){
#     myData$P1 = 2*pnorm(-abs(obj$all_BFs$observed_Z))
#   } else if(P1=="Prior"){
#     ZPrior = obj$all_BFs$prior_estimate/ obj$all_BFs$prior_std_error
#     myData$P1 = 2*pnorm(-abs(ZPrior))
#   } else if(P1=="Posterior"){
#     ZPosterior = obj$all_BFs$posterior_estimate/ obj$all_BFs$posterior_std_error
#     myData$P1 = 2*pnorm(-abs(ZPosterior))
#   } else if(P1=="Difference"){
#     ZDiff = (obj$all_BFs$observed_Z-obj$all_BFs$prior_estimate)/sqrt(1^2+obj$all_BFs$prior_std_error^2)
#     myData$P1 = 2*pnorm(-abs(ZDiff))
#   }
#   if(P2=="Obs"){
#     myData$P2 = 2*pnorm(-abs(obj$all_BFs$observed_Z))
#   } else if(P2=="Prior"){
#     ZPrior = obj$all_BFs$prior_estimate/ obj$all_BFs$prior_std_error
#     myData$P2 = 2*pnorm(-abs(ZPrior))
#   } else if(P2=="Posterior"){
#     ZPosterior = obj$all_BFs$posterior_estimate/ obj$all_BFs$posterior_std_error
#     myData$P2 = 2*pnorm(-abs(ZPosterior))
#   } else if(P2=="Difference"){
#     ZDiff = (obj$all_BFs$observed_Z-obj$all_BFs$prior_estimate)/sqrt(1^2+obj$all_BFs$prior_std_error^2)
#     myData$P2 = 2*pnorm(-abs(ZDiff))
#   }
#
#   # 2) threshold
#   myData = myData[myData$P1 < threshold | myData$P2 < threshold,]
#   myData$P1 = -log10(myData$P1)
#   myData$P2 =  log10(myData$P2)
#   d = data.frame(rs=rep(myData$rs,2),CHR=rep(myData$chrm,2), BP=rep(myData$pos,2),
#                       P=c(myData$P1, myData$P2))
#
#
#   #if (!is.null(x[[snp]]))
#   #  d = transform(d, SNP = x[[snp]])
#   #d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
#   d <- d[order(d$CHR, d$BP), ]
#
#   d$pos = NA
#   d$index = NA
#   ind = 0
#   for (i in unique(d$CHR)) {
#     ind = ind + 1
#     d[d$CHR == i, ]$index = ind
#   }
#   nchr = length(unique(d$CHR))
#   if (nchr == 1) {
#     d$pos = d$BP
#     ticks = floor(length(d$pos))/2 + 1
#     xlabel = paste("Chromosome", unique(d$CHR), "position")
#     labs = ticks
#   }  else {
#     lastbase = 0
#     ticks = NULL
#     for (i in unique(d$index)) {
#       if (i == 1) {
#         d[d$index == i, ]$pos = d[d$index == i, ]$BP
#       }      else {
#         lastbase = lastbase + tail(subset(d, index ==
#                                             i - 1)$BP, 1)
#         d[d$index == i, ]$pos = d[d$index == i, ]$BP +
#           lastbase
#       }
#       ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==
#                                                              i, ]$pos))/2 + 1)
#     }
#     xlabel = "Chromosome"
#     labs <- unique(d$CHR)
#   }
#   xmax = ceiling(max(d$pos) * 1.03)
#   xmin = floor(max(d$pos) * -0.03)
#   def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
#                    las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(ceiling(-max(d$P)),
#                                                                      ceiling(max(d$P))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
#   dotargs <- NULL
#   do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
#                                             names(dotargs)]))
#
#     axis(1, at = ticks, labels = labs)
#
#   col = rep(c("grey", "black"), max(d$CHR))
#   if (nchr == 1) {
#     with(d, points(pos, P, pch = 20, col = col[1], ...))
#   } else {
#     icol = 1
#     for (i in unique(d$index)) {
#       with(d[d$index == unique(d$index)[i], ], points(pos,
#                                                       P, col = col[icol], pch = 20))
#       icol = icol + 1
#     }
#   }
#
#
#   if (suggestiveline)
#     abline(h = suggestiveline, col = "blue")
#   if (genomewideline)
#     abline(h = genomewideline, col = "red")
#   if (!is.null(highlight)) {
#     if (any(!(highlight %in% d$SNP)))
#       warning("You're trying to highlight SNPs that don't exist in your results.")
#     d.highlight = d[which(d$SNP %in% highlight), ]
#     with(d.highlight, points(pos, logp, col = "green3", pch = 20,
#                              ...))
#   }
#   if (!is.null(annotatePval)) {
#     topHits = subset(d, P <= annotatePval)
#     par(xpd = TRUE)
#     if (annotateTop == FALSE) {
#       with(subset(d, P <= annotatePval), textxy(pos, -log10(P),
#                                                 offset = 0.625, labs = topHits$SNP, cex = 0.45),
#            ...)
#     }
#     else {
#       topHits <- topHits[order(topHits$P), ]
#       topSNPs <- NULL
#       for (i in unique(topHits$CHR)) {
#         chrSNPs <- topHits[topHits$CHR == i, ]
#         topSNPs <- rbind(topSNPs, chrSNPs[1, ])
#       }
#       textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625,
#              labs = topSNPs$SNP, cex = 0.5, ...)
#     }
#   }
#   par(xpd = FALSE)
# }

get_RSquared <- function(obj){ # obj should be a bGWAS object
  Line =  obj$log_info[stringr::str_detect(obj$log_info, "Median out-of-sample adjusted R-squared across all chromosomes is ")]
  R2 = as.numeric(stringr::str_split(Line, "is ")[[1]][2])
  return(R2)
}


