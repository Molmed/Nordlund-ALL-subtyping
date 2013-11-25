##' Machine level CpG-site filters
##' 
##' @param x Dataset. Must be numeric matrix.
##' @param groups Vector of group memberships. Only two groups labeled \code{1}
##'   and \code{2} are allowed. Use \code{NA} to exclude from calculation.
##' @param feats Feature subset, logical vector. Excluded features will be
##'   marked with \code{NA} in the results.
##' @return A vector of delta beta values.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
column.max.delta.beta <- function(x, groups, feats){
    if(is.logical(groups)) groups <- groups + 1
    groups <- na.fill(as.integer(groups), 0)
    if(length(groups) != nrow(x))
        stop("Group vector does not match number of rows in x.")
    if(missing(feats)) feats <- rep(TRUE, ncol(x))
    if(!is.logical(feats) || any(is.na(feats)))
        stop("Feature subset must be a logical vector without NAs.")
    if(!is.numeric(x)) stop("x must be numeric.")
    .Call("column_max_delta_beta", PACKAGE="nordlund2013", x, groups, feats)
}

##' @param thres NA frequency threshold.
##' @param subset Subset of rows to inspect.
##' @return A logical vector describing which columns have fewer \code{NA} than
##'   \code{thres}.
##' @rdname column.max.delta.beta
##' @export
column.na.frac <- function(x, thres, subset, feats){
    if(!is.numeric(x)) stop("x must be numeric.")
    if(!is.numeric(thres) || thres < 0 || thres >= 1)
        stop("`thres` must be a number in [0, 1).")
    if(missing(subset))
        subset <- rep(TRUE, ncol(x))
    if(!is.logical(subset) || any(is.na(subset)))
        stop("Subset must be a logical vector without NAs.")
    if(missing(feats)) feats <- rep(TRUE, ncol(x))
    if(!is.logical(feats) || any(is.na(feats)))
        stop("Feature subset must be a logical vector without NAs.")
    .Call("column_na_frac", PACKAGE="nordlund2013", x, thres, subset, feats)
}
