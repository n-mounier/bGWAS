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

     # data should be ordered by chrm/pos
     data <- data[order(data[,2], data[,3])]
     # order the SNPs from the best, to the worst
     if(byP){
       SNP_order = order(data[,4])
     } else {
       SNP_order =  order(-data[,4])
     }

     "%|%" <- function(x,y){
       do.call(y,list(substitute(x)),envir=parent.frame()) # just right
     }

     incrementer = 0
     num.pruned = 0
     #deal with SNPs, from the most significant, to the less significant
     for(i in SNP_order) {
      if(is.na(data[i,4])) { # SNP already pruned
        next
      }
       incrementer = incrementer + 1

       # Best SNP of this region
       rs   = as.character(data[i,1])
       chrm = as.numeric(data[i,2])
       pos  = as.numeric(data[i,3])

       # Find the all the SNPs in the prune.dist-kb window
       suppressWarnings({
         gtools::binsearch( function(x) { # lower
           data[x,2] * 1000000000 + data[x,3]
         } , c(1, nrow(data))
         , target =       chrm * 1000000000 + pos-(prune.dist*1000)
         )$where -> lb
         if(length(lb)==2) {
           lb = lb[2]
         }
         gtools::binsearch( function(x) { # upper
           data[x,2] * 1000000000 + data[x,3]
         } , c(1, nrow(data))
         , target =       chrm * 1000000000 + pos+(prune.dist*1000)
         )$where -> ub
         if(length(ub)==2) {
           ub = ub[1]
         }
       })

       relevant.offsets = c(lb:ub)
       already.pruned = is.na(data[relevant.offsets,4])
       relevant.positions = as.numeric(unlist(data[relevant.offsets,3]))
       relevant.positions[already.pruned] <- NA

       # start looking at these SNPs
       current = lb
       count_them = 0

       while(current <= ub) {
         # check that the SNP is in the window
         stopifnot(chrm == data[current,2])
         stopifnot(as.numeric(data[current,3]) >= pos - (prune.dist*1000))
         stopifnot(as.numeric(data[current,3]) <= pos + (prune.dist*1000))

         # if this is the SNP of interest (i.e. the best one in the window)
         # don't do anything, but check that it's not NA
         if( data[current,1] == rs ) {
           (!is.na(data[current,4])) %|%stopifnot
         }

         else {

         # already pruned, happy days
         if(is.na( data[current,4] )) {
           current = current+1
           next
         }

        # prune it (set stat/p to NA)
         if(!is.na(data[current,4] )) {
           data[current,4] <- NA
           count_them = count_them + 1
           num.pruned = num.pruned + 1
         }

         }

       current = current+1 # go to the next one in the window
     }
     data[i,4] %|%na.fail # make sure I haven't deleted 'myself'
}


     return(unlist(data[!is.na(unlist(data[,4])),1]))

}
