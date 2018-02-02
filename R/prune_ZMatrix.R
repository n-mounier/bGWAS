###### Function to distance pruned SNPs ######



#' Prune a Z-Matrix
#'
#' @inheritParams makeMR_ZMatrix
#' @param ZMatrix A subset of the prior GWASs Z-Matrix
#' @param prune.dist The distance for pruning
#' @param r2.limit The r2 threshold if LD-pruning, \code{default=1.1}, only distance pruning
#'        is performed
#'
#' @details
#' To perform LD-pruning, need binary files !!!
#'
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":="
#  Function not exported, no need of extended documentation?


prune_ZMatrix <- function(ZMatrix, prune.dist=500000, r2.limit=1.1, verbose=F) {

  "%|%" <- function(x,y){
    do.call(y,list(substitute(x)),envir=parent.frame()) # just right
  }


  # remove non-numeric columns
  ZMatrix[,!(colnames(ZMatrix) %in% c('chrm','rs','pos', 'alt', 'ref')),with=F] -> ZMatrix.justZMatrix


  # get the largest Z for each SNP
  apply(ZMatrix.justZMatrix, 1, function(ZMatrix){
    max(abs(ZMatrix))
  }) -> max.ZMatrix

  # convert it to log-p to compare SNPs (calculated directly from z-scores)
  log(2.0)+pnorm(-abs(max.ZMatrix), log=T)    -> logps


  # order the log-p to start pruning by looking at the most significant SNPs
  order(logps) -> order.logps


  incrementer = 0
  num.pruned = 0
# deal with SNPs, from the most significant, to the less significant
  for(i in order.logps) {
    if(is.na(logps[i])) { # why would it happen ?
      next
    }
    incrementer = incrementer + 1
#    PP(i, incrementer, num.pruned, (incrementer+num.pruned)/length(logps))
    rs   = ZMatrix$rs  [i]
    chrm = ZMatrix$chrm[i]
    pos  = ZMatrix$pos [i]
#    PP(rs, chrm, pos, logps[i], max.ZMatrix[i], max(abs(ZMatrix[i,-(1:3)])))
    logps[i] %|%na.fail # this shouldn't happen since we are explicitely "nexting" the NA
    suppressWarnings({
      gtools::binsearch( function(x) {
        ZMatrix$chrm[x] * 1000000000 + ZMatrix$pos[x]
      } , c(1, nrow(ZMatrix))
      , target =       chrm * 1000000000 + pos-prune.dist
      )$where -> lb
      if(length(lb)==2) {
        stopifnot(lb[1]<lb[2])
        lb = lb[2]
      }
      gtools::binsearch( function(x) {
        ZMatrix$chrm[x] * 1000000000 + ZMatrix$pos[x]
      } , c(1, nrow(ZMatrix))
      , target =       chrm * 1000000000 + pos+prune.dist
      )$where -> ub
      if(length(ub)==2) {
        stopifnot(ub[1]<ub[2])
        ub = ub[1]
      }
    })
    stopifnot(length(lb)==1)
    stopifnot(length(ub)==1)

    current = lb
    count_them = 0

    my.seq.int <- function(from, to, by = 1) {
      if(to < from) {
        stopifnot(to+by == from)
        rep.int(0L,0)
      } else {
        seq.int(from=from,to=to,by=by)
      }
    }

    my.seq.int( lb, ub ) -> ld.relevant.offsets
    is.na(logps[ld.relevant.offsets]) -> already.pruned
   # PP       (length(already.pruned) ,  ub-lb)
    stopifnot(length(already.pruned) == ub-lb+1)
    ld.relevant.positions = ZMatrix$pos[ld.relevant.offsets]
    ld.relevant.positions[already.pruned] <- NA
    stopifnot( ld.relevant.positions[ZMatrix$rs[ld.relevant.offsets] == rs] == pos)
    if(r2.limit != 1.1) { # if LD pruning
      cat(chrm)
      cat(pos)
      stop("LD pruning not implemented")
      lookup.LD.inUK10K(chrm, pos, ld.relevant.positions) -> many.rs
    } else { # if distance pruning
      NULL -> many.rs
    }
    while(current <= ub) {
      stopifnot(chrm == ZMatrix$chrm[current])
      stopifnot(ZMatrix$pos [current] >= pos - prune.dist)
      stopifnot(ZMatrix$pos [current] <= pos + prune.dist)
      if( ZMatrix$rs[current] == rs ) {
        (!is.na(logps[current])) %|%stopifnot
        stopifnot(1.0 == many.rs[current-lb+1])
      }
      else {
        if(is.na( logps[current] )) {
          # already pruned. happy days
          current = current+1
          next
        }

        (!is.na(logps[current])) %|%stopifnot

        # might have to apply LD-based pruning
        r2 <- NULL
        if(r2.limit != 1.1) {
          many.rs[current-lb+1] -> r
          na.fail(r)
          r2 <- r*r
        } else { # if distance pruning
          r2 <- -0.1234 # just for pretend
        }

        if(r2.limit == 1.1 || r2 >= r2.limit) {
          if(!is.na( logps[current] )) {
            logps[current] <- NA
            count_them = count_them + 1
            num.pruned = num.pruned + 1
          }
        }
      }
      current = current+1
    }
    logps[i] %|%na.fail # make sure I haven't deleted 'myself'
  }


  stopifnot(length(logps) == nrow(ZMatrix))
  ZMatrix    = ZMatrix   [!is.na(logps),,drop=F]
  logps = logps[!is.na(logps)]

  #if only one row and alt or ref is "T", R will consider the column as a boolean "TRUE" and not a character T, should not happen in real life, but when testing with a really small set of SNPs, this is possible...
  ZMatrix$ref[ZMatrix$ref==TRUE]="T"


  return(ZMatrix)

  }
