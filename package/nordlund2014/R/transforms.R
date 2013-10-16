##' @import predictBase
{}

##' Feature selection based on beta values and kNN-imputation
##' 
##' Dataset preprocessing transform that
##' \tabular{l}{
##'   removes sites with differences in groupwise mean beta values < \code{dbeta} \cr
##'   removes sites with more missing values than \code{na.frac} \cr
##'   applies kNN-imputation
##' }
##' 
##' For more information on how dataset preprocessing transforms work, see
##' \code{\link{pre.trans}}.
##' 
##' @param x Dataset.
##' @param y Class label vector.
##' @param fold Resampling fold.
##' @param dbeta Delta beta cutoff. Only sites with a difference in group means
##'   greater than this value will be included in the model.
##' @param na.frac NA fraction cutoff. Only sites with a fraction of missing
##'   values lower than value will be included in the model.
##' @param k Sent to \code{\link{pre.impute.knn}}. Set to NA to suppress
##'   imputation.
##' @param feat An initial feature filtering vector.
##' @param cluster A snow cluster to be used for parallelization.
##' @param distmat Sent to \code{\link{pre.impute.knn}}.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
pre.trans.450k <- function(x, y, fold, dbeta, na.frac, k=.05, distmat, feat=rep(TRUE, ncol(x)), cluster){
    if(!missing(dbeta)){
        y.int <- as.integer(y)
        y.int[na.fill(fold, TRUE)] <- 0L
        feat <- feat & column.max.delta.beta(x, y.int, feat) > dbeta
    }
    if(!missing(na.frac)){
        feat <- feat & column.na.frac(x, na.frac, na.fill(!fold, FALSE), feat)
    }
    if(!any(feat)) stop("No sites passed the filters.")

    #if(is.character(distmat) && distmat == "auto") x <- x[,feat,drop=FALSE]
    c(pre.impute.knn(x=x[,feat,drop=FALSE], y=y, fold=fold, k=k, distmat=distmat),
      list(features=feat))
}

