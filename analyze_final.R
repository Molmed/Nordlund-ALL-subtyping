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
number.of.cores <- 3

# If we run on the UPPMAX cluster, maximize the number of processes
if(grepl("^[qm]\\d+\\.uppmax\\.uu\\.se$", Sys.info()["nodename"])){
    max.mem <- as.integer(sub("^MemTotal:\\s+(\\d+) kB$", "\\1",
        system("head /proc/meminfo -n 1", intern=TRUE)))
    number.of.cores <- floor(max.mem/12e6)
}

library(parallel)
options(mc.cores = number.of.cores)


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
#   On top this distinction, a model parameter `times.chosen` is also tuned,
#   referred to as `F` in the paper. In controls the number of folds a variable
#   must be selected in to qualify for the final model.
#
#-------------------------------------------------------------------------------

out.file <- "results/tuning.Rdata"
if(file.exists(out.file)){
    load(out.file)
} else {
    n.sites <- array(NA, c(length(y), length(cv), length(cv)))
    probs <- list(separate=vector("list", length(cv)),
                  combined=vector("list", length(cv)))
    error <- error.sex <-
        list(separate = matrix(NA, length(cv), length(cv)),
             combined = matrix(NA, length(cv), length(cv)))
    n.error <- list(separate = array(NA, c(2, length(cv), length(cv))),
                    combined = array(NA, c(2, length(cv), length(cv))))
    conf.tab <- list(separate = array(NA, c(2, 2, length(y), length(cv), length(cv))),
                     combined = array(NA, c(2, 2, length(y), length(cv), length(cv))))
    #sens.spec <- true.call.frac <-
    #    list(separate = array(NA, c(2, length(y), length(cv), length(cv))),
    #         combined = array(NA, c(2, length(y), length(cv), length(cv))))
}
cv.feat.sel <- vector("list", length(cv))
save.assembly <- function()
    save(n.sites, probs, error, error.sex, n.error, conf.tab, file=out.file)

trace.msg(1, "Confirming that the tuning has completed", linebreak=FALSE)
for(i in seq_along(cv)){
    tryCatch({
        load(sprintf("tuning/fold_%i.Rdata", i))
        fold.idx <- sapply(fold.feat.sel, is.null)
        stopifnot(!any(fold.idx))
        cv.feat.sel[[i]] <- fold.feat.sel
        rm(fold.feat.sel)
    }, error=function(err){
        stop(sprintf("Model tuning has not completed successfully, please rerun fold %i %s.",
                     i, paste(names(y)[fold.idx], collapse=", ")))
    })
    n.sites[,i,] <- t(sapply(cv.feat.sel[[i]], function(x)
            rev(cumsum(rev(unname(table(cut(x, 0:25))))))))
    cat(".")
}
cat("\n")
save.assembly()
#save(cv.feat.sel, n.sites, file="results/cv_feat_sel.Rdata")
#load("results/cv_feat_sel.Rdata")


trace.msg(1, "Fitting classifiers and predicting test sets classes")
# Overlayer the default predict function with one that discards as few sites
# as possible, since we already tune complexity with `times.chosen`.
predict.nsc <- function(...)
    predict:::predict.nsc(..., thres=min)
for(i in seq_along(cv)){
    trace.msg(2, "Fold %i", i)
    need.save <- FALSE

    # Train classifiers and calculate test set probabilities
    for(method in names(probs)){
        if(is.blank(probs[[method]][[i]])){
            trace.msg(3, "Classifying based on %s feature sets ", method, linebreak=FALSE)
            
            if(method == "combined")
                v <- Reduce("pmax", cv.feat.sel[[i]])
            probs[[method]][[i]] <- vector("list", 25)
            for(times.chosen in 1:25){
                my.pred <- mclapply(names(y), function(my.class){
                    my.site.idx <- switch(method,
                        separate = cv.feat.sel[[i]][[my.class]] >= times.chosen,
                        combined = v >= times.chosen)
                    if(!any(my.site.idx) || sum(my.site.idx) > 10000) return(NULL)
                    my.sample.idx <- cv[[i]] | !is.na(y[[my.class]])

                    my.met <- met.data[my.sample.idx, my.site.idx, drop=FALSE]
                    my.y <- y[[my.class]][my.sample.idx]
                    my.y <- na.fill(my.y, levels(my.y)[2])
                    batch.predict(my.met, my.y,
                        models = list(nsc=list(cv=list(list(nrep=25, nfold=8)))),
                        # 8 comes from that there are only 8 confirmed iAMP samples
                        test.subset = cv[my.sample.idx,][i],
                        pre.trans = pre.impute.median)
                })
                names(my.pred) <- names(y)
                cat(".")
                probs[[method]][[i]][[times.chosen]] <- if(any(sapply(my.pred, is.null))){
                    NA
                } else {
                    sapply(my.pred, function(p) p$cv[[1]]$nsc$prob[,1])
                }
            }
            cat("\n")
            need.save <- TRUE
        }
    }
    if(need.save) save.assembly()
}
rm(predict.nsc)


trace.msg(1, "Evaluating performance and tuning parameters")
truth <- sapply(y, function(my.y) as.integer(my.y) == 1)
truth[y$reference %in% "reference", 2:9] <- FALSE

for(method in names(error)){
    trace.msg(2, "Evaluating performance of the %s method ", method, linebreak=FALSE)
    for(i in seq_along(cv)){
        for(times.chosen in seq_along(cv)){
            if(!is.blank(probs[[method]][[i]][[times.chosen]])){

                # Call classes from probability estimates
                p <- probs[[method]][[i]][[times.chosen]] >= .5
                # Subtype classification is irrelevant when classed as reference
                p[p[,1],2:9] <- FALSE
                # Samples of unknown classes can not be part of subtype performance evaluation
                p[is.na(truth[cv[[i]],1]), 1:9] <- NA
                # Sex does not apply for purified reference samples, e.g. CD19+
                p[is.na(truth[cv[[i]],"sex"]), "sex"] <- NA 
                correct <- p == truth[cv[[i]],]

                n.error[[method]][,i,times.chosen] <- table(factor(
                    !apply(correct[,1:9], 1, all, na.rm=TRUE)[!apply(is.na(correct[,1:9]), 1, all)],
                    levels=c(FALSE, TRUE)))
                error[[method]][i, times.chosen] <- prop.table(n.error[[method]][,i,times.chosen])[2]
                error.sex[[method]][i, times.chosen] <- mean(!correct[,10], na.rm=TRUE)

                conf.tab[[method]][,,,i,times.chosen] <- 
                    mapply(function(yt, yp){
                        table(yt, factor(yp, levels=c(TRUE, FALSE)))
                    }, y[cv[[i]],], as.data.frame(p))
            }
        }
        cat(".")
    }
    cat("\n")
}
save.assembly()


#===============================================================================
#   Perform final classification
#
#   This section took about 8 h to run, so we decided not to bother with
#   parallelization.
#-------------------------------------------------------------------------------

#merr <- sapply(error, apply, 2, mean)
#which(merr == min(merr, na.rm=TRUE), arr.ind=TRUE)

times.chosen <- 17


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
    data.frame(Subtype = factor(my.class, levels=names(y)),
               TargetID = met.annot$TargetID[feat.sel[[my.class]] >= times.chosen],
               stringsAsFactors=FALSE)
}))
write.table(cons.sites, "results/consensus_sites.csv", quote=FALSE, sep=",",
            row.names=FALSE)


#-------------------------------o
#   Train consensus classifiers

#cons.met <- list(
#    impute.knn(met.data[, met.annot$TargetID %in% cons.sites$TargetID &
#                          met.annot$CHR %in% 1:22], distmat="auto"),
#    impute.knn(met.data[,met.annot$TargetID %in% cons.sites$TargetID],
#               distmat="auto"))
#for(my.class in names(y)){
#    cat("Making final classifier for", my.class, "\n")
#    idx <- !is.na(y[[my.class]])
#    cons[[my.class]] <- design("nsc",
#        cons.met[[if(my.class == "sex") 2 else 1]][idx,], y[[my.class]][idx],
#        cv = list(nrep=25, nfold=8),
#        slim.fit = TRUE)
#}
cons.met <- lapply(with(cons.sites, split(TargetID, Subtype)),
    function(ids) impute.knn(met.data[,met.annot$TargetID %in% ids], distmat="auto"))
for(my.class in names(y)){
    cat("Making final classifier for", my.class, "\n")
    idx <- !is.na(y[[my.class]])
    cons[[my.class]] <- design("nsc",
        cons.met[[my.class]][idx,], y[[my.class]][idx],
        thres = 0, cv = list(nrep=25, nfold=8),
        slim.fit = TRUE)
}


#-------------------------------o
#   Predict class probabilities of all samples based on the selected sites,
#   including the samples used for training.

#pred <- mapply(predict, cons, cons.met[rep(1:2, c(9,1))],
#               MoreArgs=list(threshold=min), SIMPLIFY=FALSE)
pred <- mapply(predict, cons, cons.met,
               MoreArgs=list(threshold=min), SIMPLIFY=FALSE)
save.workspace()

cons.pred <- data.frame(met.pheno,
    lapply(pred, function(p) p$prob[,1]))
names(cons.pred)[ncol(cons.pred)] <- "sex.female"
write.table(cons.pred, "results/consensus_predictions.csv",
    quote=FALSE, sep=",", row.names=FALSE)


#-------------------------------o
#   Predict class of the blinded validation samples that has not been part of
#   any previous stage of the analysis

library(analyse450k)
load.450k.data("subtype_validation", complete=TRUE)
val.met <- t(val.met[site.idx,])

#val.cons <- list(
#    impute.knn(val.met[, met.annot$TargetID %in% cons.sites$TargetID &
#                         met.annot$CHR %in% 1:22], distmat="auto"),
#    impute.knn(val.met[, met.annot$TargetID %in% cons.sites$TargetID],
#               distmat="auto"))
val.cons <- lapply(with(cons.sites, split(TargetID, Subtype)),
    function(ids) impute.knn(val.met[,met.annot$TargetID %in% ids], distmat="auto"))

val.pred <- data.frame(val.pheno,
    mapply(function(fit, dat) predict(fit, dat, threshold=min)$prob[,1],
           cons, val.cons[rep(1:2, c(9,1))], SIMPLIFY=FALSE),
    check.names=FALSE)

write.table(val.pred, "results/validation_predictions.csv",
    quote=FALSE, sep=",", row.names=FALSE)

save.workspace()


#===============================================================================
#   Make a summary table
#-------------------------------------------------------------------------------

sens <- apply(conf.tab$combined[1,,,,times.chosen], 2:3, prop.table)[1,,]
spec <- apply(conf.tab$combined[2,,,,times.chosen], 2:3, prop.table)[2,,]

write.table(
    data.frame(
        subtype = names(y),
        sens.mean = apply(sens, 1, mean),
        sens.sd   = apply(sens, 1, sd),
        spec.mean = apply(spec, 1, mean),
        spec.sd   = apply(spec, 1, sd),
        n.min     = apply(n.sites[,,times.chosen], 1, min),
        n.max     = apply(n.sites[,,times.chosen], 1, max),
        n.cc      = sapply(feat.sel, function(x) sum(x >= times.chosen))),
    file="results/sens_spec_table.csv", sep=",", row.names=FALSE)

