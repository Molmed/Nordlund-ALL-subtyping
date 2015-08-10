#!/usr/bin/Rscript

#===============================================================================
#  The final subtype classifier
#
#    To use it source this script, run the setup function and call the returned
#    function with your own methylation data.
#
#  Usage:
#
#    source("classifier.R")
#    f <- get.consensus.classifier()
#    f(met.data, met.annot$TargetID)
#
#  Arguments (to `f` above):
#
#             x: Methylation data.
#
#      TargetID: Illumina probe IDs defining the CpG-sites of `x`.
#
#    samples.as: Whether `x` contains samples as rows or columns.
#
#-------------------------------------------------------------------------------

get.consensus.classifier <- function(){
    if(!require("pamr")){
        cat("The `pamr` package is required for running the consensus classifier. Would you like to install it? [Y/n]\n")
        if(!tolower(readline()) %in% c("", "y"))
            stop("Cannot run classifier without the `pamr` package.")
        install.packages("pamr")
        require("pamr")
    }

    # Load the `cons` object
    load("classifier.Rdata")

    # Load the consensus sites
    cons.sites <- read.csv("consensus_sites.csv", colClasses=c("factor", "character", "factor", "integer"))

    # The classifier to be returned to the user
    function(x, TargetID, samples.as=c("rows", "columns")){
        ind <- match(cons.sites$TargetID, TargetID)
        if(!length(TargetID) %in% dim(x))
            stop("`x` and `TargetID` does not match.")

        # Orient dataset
        samples.as <- match.arg(samples.as)
        if(samples.as == "rows"){
            x <- t(x[,ind])
        } else {
            x <- x[ind,]
        }
        na.ind <- which(is.na(x), arr.ind=TRUE)
    
        probability <- vapply(cons, function(subtype.classifier){
            # Remove sites on the X-chromosome if we are not classifying the sex
            if(length(subtype.classifier$centroid.overall) < nrow(x)){
                x <- x[!cons.sites$Subtype %in% "sex",]
                na.ind <- na.ind[!cons.sites$Subtype[na.ind[,"row"]] %in% "sex",]
            }

            # Impute missing values by the overall centroid,
            # causing them to not favour any class
            if(nrow(na.ind) > 0)
                x[na.ind] <- subtype.classifier$centroid.overall[na.ind[,"row"]]
            pamr.predict(subtype.classifier, x, type = "posterior",
                         threshold = subtype.classifier$threshold)[,1]
        }, numeric(ncol(x)))
        colnames(probability)[10] <- "female"

        i <- !colnames(probability) %in% c("reference", "sex")
        subtype.hits <- probability[,i] >= .5
        n.hits <- apply(subtype.hits, 1, sum)
        data.frame(subtype = factor(
            ifelse(probability[,"reference"] >= .5, "non-leukemic",
            ifelse(n.hits == 0, "unknown",
            ifelse(n.hits > 1, "multiple",
                colnames(probability)[i][apply(probability[,i], 1, which.max)]))),
            levels=c(colnames(probability)[i], "multiple", "unknown", "non-leukemic")),
            sex = factor(ifelse(probability[,"female"] >= .5, "female", "male")),
            probability
        )
    }
}

