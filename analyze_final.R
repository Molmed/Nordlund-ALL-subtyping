#===============================================================================
#
#   This script performs the analyses that the paper is based on.
#   
#   Since it takes quite some time to carry out all of it, it is written in a
#   way that allows halting and resuming without substantial loss of results.
#
#-------------------------------------------------------------------------------

source("analyze_init.R")

# TODO: remove tryCatches

# TEMP
site.ind <- 1:1000
met.data <- t(met.data[site.ind, sample.idx])
met.annot <- met.annot[site.ind,]


#===============================================================================
#   Perform classification
#
#   This section took about 8 h to run, so we decided not to bother with
#   parallelization.
#-------------------------------------------------------------------------------


#-------------------------------o
#   Select sites for each class

for(my.class in names(y)[sapply(feat.sel, is.null)]){
    tryCatch({
        cat("\nSelecting features for", my.class, "\n")
        feat.sel[[my.class]] <- design("feature_selection", met.data, y[[my.class]],
                        chr = if(my.class == "sex") c(1:22, "X") else 1:22, cv = cv)
        save.workspace()
    }, error=function(err){
        print(err)
    })
}
if(any(sapply(feat.sel, is.null)))
    stop("Encountered problems during feature selection. Please fix and complete before continuing.")


#-------------------------------o
#   Output a summary of the selected sites

cons.sites <- do.call(rbind, lapply(names(y), function(my.class){
    data.frame(Subtype = my.class,
               TargetID = met.annot$TargetID[feat.sel[[my.class]] >= 10],
               stringsAsFactors=FALSE)
}))
write.table(cons.sites, "results/consensus_sites.csv", quote=FALSE, sep="\t",
            row.names=FALSE)


#-------------------------------o
#   Train consensus classifiers

cons.met <- impute.knn(met.data[,met.annot$TargetID %in% cons.sites$TargetID],
                       distmat="auto")
for(my.class in names(y)[sapply(cons, is.null)]){
    cat("Making final classifier for", my.class, "\n")
    idx <- !is.na(y[[my.class]])
    cons[[my.class]] <- design("nsc", cons.met[idx,], y[[my.class]][idx])
}


#-------------------------------o
#   Predict class probabilities of all samples based on the selected sites,
#   including the samples used for training.

pred <- lapply(cons, predict, cons.met)
save.workspace()
cons.pred <- data.frame(met.pheno,
    lapply(pred, function(p) p$prob[,1]))
names(cons.pred)[ncol(cons.pred)] <- "sex.female"
write.table(cons.pred, "results/consensus_predictions.csv",
    quote=FALSE, sep="\t", row.names=FALSE)


