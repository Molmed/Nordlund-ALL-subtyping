#!/opt/apps/R/3.0.1/bin/Rscript

#===============================================================================
#
#   This script performs the analyses that the paper is based on.
#   
#   Since it takes quite some time to carry out all of it, it is written in a
#   way that allows halting and resuming without substantial loss of results.
#
#-------------------------------------------------------------------------------

source("analyze_init.R")

library(doMC)
registerDoMC(3)


#===============================================================================
#   Assemble results from parameter tuning
#-------------------------------------------------------------------------------
#
#   Two alternative feature selection methods were considered.
#
#    1. Only let the subtype classifiers use the sites that has been selected
#       for its particular subtype. Variables belonging to this approach are
#       named `sep.###` (for "separate").
#    2. Let all classifiers use all sites selected for any subtype. Varibles
#       belonging to this approach are named `comb.###` (for combined).
#
#   On top this distinction, two modelparameters are also tuned.
#
#    1. The number of folds a variable must be selected in to qualify for the
#       final model (`times.selected`).
#    2. The threshold to use for calling a subtype (`thres`).
#
#-------------------------------------------------------------------------------

out.file <- "results/tuning.Rdata"
if(file.exists(out.file)){
    load(out.file)
} else {
    thres.axs <- 0:100/100
    sep.probs <- sep.errors <- comb.probs <- comb.errors <- vector("list", length(cv))
}
cv.feat.sel <- vector("list", length(cv))
save.assembly <- function()
    save(sep.probs, sep.errors, comb.probs, comb.errors, thres.axs,
         file=out.file)

trace.msg(1, "Confirming that the tuning has completed", linebreak=FALSE)
for(i in seq_along(cv)){
    tryCatch({
        load(sprintf("tuning/fold_%i.Rdata", i))
        stopifnot(!any(sapply(fold.feat.sel, is.blank)))
        cv.feat.sel[[i]] <- fold.feat.sel
        rm(fold.feat.sel)
    }, error=function(err){
        stop("Model tuning has not completed successfully.")
    })
    cat(".")
}
cat("\n")


trace.msg(1, "Assembling tuning results")
for(i in seq_along(cv)){
    trace.msg(2, "Fold %i", i)

    # Train classifiers and calculate test set probabilities using the separate approach
    if(is.blank(sep.probs[[i]])){
        trace.msg(3, "Classifying based on separate datasets ", linebreak=FALSE)
        sep.probs[[i]] <- foreach(times.chosen=1:25) %dopar% {
            my.pred <- lapply(seq_along(y), function(j){
                idx <- cv.feat.sel[[i]][[j]] >= times.chosen
                if(!any(idx)) return(NA)
                my.met <- met.data[, idx, drop=FALSE]
                my.y <- y[[j]]
                my.y[is.na(my.y) & cv[[i]]] <- levels(my.y)[2]
                batch.predict(my.met, my.y, "nsc", cv[i],
                    pre.trans = pre.impute.median)
            })
            names(my.pred) <- names(y)
            cat(".")
            if(any(sapply(my.pred, is.null))) return(NA)
            sapply(my.pred, function(p) p$cv[[1]]$nsc$prob[,1])
        }
        cat("\n")
    }
    if(is.blank(sep.errors[[i]])){
        trace.msg(3, "Calculating errors ", linebreak=FALSE)
        sep.errors[[i]] <- foreach(j = 1:25, .combine=cbind) %dopar% {
            if(is.blank(sep.probs[[i]][[j]])) return(NA)
            tmp <- foreach(thres = thres.axs, .combine=c) %do% {
                preds <- lapply(apply(sep.probs[[i]][[j]], 1, function(x)
                                as.list(which(x >= thres))), names)
                correct <- mapply(function(s, p){
                    if(s == "reference"){
                        "reference" %in% p
                    } else if(s %in% class.types){
                        length(intersect(p, class.types)) == 1 & s %in% p
                    } else {
                        NA
                    }
                }, met.pheno$subtype[cv[[i]]], preds)
                mean(!correct, na.rm=TRUE)
            }
            cat(".")
            tmp
        }
        cat("\n")
    }

    # And again using the combined approach
    if(is.blank(comb.probs[[i]])){
        trace.msg(3, "Classifying based on combined dataset ", linebreak=FALSE)
        v <- Reduce("pmax", cv.feat.sel[[i]])
        comb.probs[[i]] <- foreach(times.chosen=1:25) %dopar% {
            my.met <- met.data[, v >= times.chosen]
            my.pred <- lapply(y, function(my.y){
                my.y[is.na(my.y) & cv[[i]]] <- levels(my.y)[2]
                batch.predict(my.met, my.y, "nsc", cv[i],
                    pre.trans = pre.impute.median)
            })
            cat(".")
            if(any(sapply(my.pred, is.null))) return(NULL)
            sapply(my.pred, function(p) p$cv[[1]]$nsc$prob[,1])
        }
        cat("\n")
    }
    # Calculate errors
    if(is.blank(comb.errors[[i]])){
        trace.msg(3, "Calculating errors ", linebreak=FALSE)
        comb.errors[[i]] <- foreach(j = 1:25, .combine=cbind) %dopar% {
            if(is.blank(comb.probs[[i]][[j]])) return(NA)
            tmp <- foreach(thres = thres.axs, .combine=c) %do% {
                if(is.null(comb.probs[[i]][[j]])) return(NA)
                preds <- lapply(apply(comb.probs[[i]][[j]], 1, function(x)
                                as.list(which(x >= thres))), names)
                correct <- mapply(function(s, p){
                    if(s == "reference"){
                        "reference" %in% p
                    } else if(s %in% class.types){
                        length(intersect(p, class.types)) == 1 & s %in% p
                    } else {
                        NA
                    }
                }, met.pheno$subtype[cv[[i]]], preds)
                mean(!correct, na.rm=TRUE)
            }
            cat(".")
            tmp
        }
        cat("\n")
    }

    save.assembly()
}


#   At this point we found that there was only marginal difference in
#   performance between the different models and choices of parameters (unless
#   `thres` is far from 0.5). The best model turned out to be feature selection
#   based on a 14-fold overlap and a classification threshold of 0.47.


#===============================================================================
#   Perform final classification
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
               TargetID = met.annot$TargetID[feat.sel[[my.class]] >= 14],
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


#-------------------------------o
#   Predict class of the blinded validation samples that has not been part of
#   any previous stage of the analysis

library(analyse450k)
load.450k.data("subtype_validation", complete=TRUE)

val.cons <- impute.knn(t(val.met[met.annot$TargetID %in% cons.sites$TargetID,]),
                       distmat="auto")
val.pred <- data.frame(val.pheno,
    lapply(cons, function(fit) predict(fit, val.cons)$prob[,1]))

write.table(val.pred, "results/validation_predictions.csv",
    quote=FALSE, sep="\t", row.names=FALSE)


