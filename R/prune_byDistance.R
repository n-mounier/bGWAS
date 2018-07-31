###### Function to distance-prune SNPs ######



# #' Prune a Z-Matrix
# #'
# #'
# #' @param ZMatrix A subset of the prior GWASs Z-Matrix
# #' @param prune.dist The distance for pruning
# #' @param r2.limit The r2 threshold if LD-pruning, \code{default=1.1}, only distance pruning
# #'        is performed
# #'
# #' @details
# #'
# #'
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":="
# #  Function not exported, no need of extended documentation?


prune_byDistance <- function(data, prune.dist=100, byP=T) {
  # data should be : 1st column rs / 2nd column chr / 3rd column pos / 4th column stat
  # if byP = T : stat = p-value -> min is better
  # if byP = F : stat = Zstat, beta.. -> max is better

     if(byP){
       SNP_order = order(data[,4])
     } else {
       SNP_order =  order(-data[,4])
     }
     data = data[SNP_order,]
     snp=0
     while(T){
       snp=snp+1
       ToRemove=which(data$chr_name==data$chr_name[snp] & abs(data$chr_start - data$chr_start[snp])<prune.dist*1000)
       if(length(ToRemove)>1){
         ToRemove = ToRemove[-1]
         data = data[-ToRemove,]
       }
       if(snp==nrow(data)) break
     }

     return(unlist(data[,1]))

}
