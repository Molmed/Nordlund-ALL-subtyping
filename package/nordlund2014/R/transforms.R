##' @import predict
{}

##' Feature selection based on beta values and imputation
##' 
##' Dataset preprocessing transform that
##' \tabular{l}{
##'   removes sites with differences in groupwise mean beta values < \code{dbeta} \cr
##'   removes sites with more missing values than \code{na.frac} \cr
##'   applies median- or kNN-imputation (depending on if \code{distmat})
##' }
##' 
##' For more information on how dataset preprocessing transforms work, see
##' \code{\link{pre.trans}}.
##' 
##' @param x Dataset.
##' @param y Class label vector.
##' @param fold Resampling fold.
##' @param feat An initial feature filtering vector.
##' @param dbeta Delta beta cutoff. Only sites with a difference in group means
##'   greater than this value will be included in the model.
##' @param na.frac NA fraction cutoff. Only sites with a fraction of missing
##'   values lower than value will be included in the model.
##' @param k Sent to \code{\link{pre.impute.knn}}. Set to NA to suppress
##'   imputation.
##' @param distmat Sent to \code{\link{pre.impute.knn}}. If missing, median
##'   imputation is used.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
pre.trans.450k <- function(x, y, fold, feat=rep(TRUE, ncol(x)), dbeta, na.frac, k=.05, distmat){
    if(!missing(dbeta)){
        y.int <- as.integer(y)
        y.int[na.fill(fold, TRUE)] <- 0L
        feat <- feat & column.max.delta.beta(x, y.int, feat) > dbeta
    }
    if(!missing(na.frac)){
        feat <- feat & column.na.frac(x, na.frac, na.fill(!fold, FALSE), feat)
    }
    if(!any(feat)) stop("No sites passed the filters.")

    if(missing(distmat)){
        return(c(pre.impute.median(x=x[,feat,drop=FALSE], y=y, fold=fold),
                 list(features=feat)))
    } else {
        return(c(pre.impute.knn(x=x[,feat,drop=FALSE], y=y, fold=fold, k=k, distmat=distmat),
                 list(features=feat)))
    }
}

