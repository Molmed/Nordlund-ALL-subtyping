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
dir.create("results", showWarnings=FALSE)


#===============================================================================
#   Assemble results from parameter tuning
#-------------------------------------------------------------------------------

tryCatch({
    for(i in seq_along(cv)){
        load(sprintf("tuning/fold_%i.Rdata", i))
        v <- Reduce("pmax", fold.feat.sel)

        prob <- foreach(times.chosen=1:25) %dopar% {
            #cons.met <- pre.impute.knn(met.data[, v >= times.chosen], fold=cv[[i]], distmat="auto")
            my.pred <- lapply(y, function(my.y){
                my.cv <- cv
                my.cv[is.na(my.y[!cv[[i]]]),] <- NA

                batch.predict(met.data[, v >= times.chosen], my.y
                fit <- design("nsc", cons.met[!test.idx,], my.y[!test.idx])
                predict(fit, cons.met)
            })
            my.prob <- sapply(my.pred, function(p) p$prob[,1])
            thres <- .5
            preds <- lapply(apply(my.prob, 1, function(x) as.list(which(x >= thres))), names)
            correct <- mapply(function(s, p){
                if(s == "reference"){
                    "reference" %in% p
                } else if(s %in% class.types){
                    length(intersect(p, class.types)) == 1 & s %in% p
                } else {
                    NA
                }
            }, met.pheno$subtype, preds)
            #for(k in head(which(!correct))){
            #    print(barchart(~my.prob[k,], main=as.character(met.pheno$subtype[k]), xlim=0:1))
            #    invisible(readline())
            #}
                 

            errors <- ifelse(y$reference == "reference",
                             my.prob[,1] < thres,
                             ifelse(my.prob[,1] >= thres,
                                    TRUE,
                                    predsmy.prob[,class.types] >= thres

        }
        reference.prob <- sapply(prob, "[[", "reference")

        sex.prob <- sapply(prob, "[[", "sex")



        my.cons.sites <- do.call(rbind, lapply(names(y), function(my.class){
            data.frame(Subtype = my.class,
                       TargetID = met.annot$TargetID[feat.sel[[my.class]] >= 10],
                       stringsAsFactors=FALSE)
        }))

        # Train consensus classifiers
        my.cons <- structure(vector("list", ncol(y)), names=names(y))
        my.cons.met <- impute.knn(
            met.data[,met.annot$TargetID %in% cons.sites$TargetID], distmat="auto")
        for(my.class in names(y)[sapply(my.cons, is.null)]){
            cat("Making final classifier for", my.class, "\n")
            idx <- !is.na(y[[my.class]]) & na.fill(!cv[[i]], FALSE)
            my.cons[[my.class]] <- design("nsc", cons.met[idx,], y[[my.class]][idx])
        }

        # Predict class probabilities of held out samples
        cv.pred[[i]] <- lapply(my.cons, predict, my.cons.met)
    }
}, error=function(err){
    print(err)
    stop("Model tuning has not completed successfully.")
})


save.workspace()



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


