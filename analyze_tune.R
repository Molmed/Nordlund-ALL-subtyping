
#===============================================================================
#
#   Parameter tuning and performance estimation
#
#   It would be criminally inefficient not to parallelize the following section.
#   However, since each worker requires about 15 GB of RAM be sure to set the
#   number of workers carefully. We used 3 nodes on Uppsala University's UPPMAX
#   cluster with only one worker each.
#
#-------------------------------------------------------------------------------

number.of.workers <- 3

source("analyze_init.R")


#-------------------------------o
#   Setup parallelization

library(doSNOW)
cl <- makeCluster(number.of.workers)
registerDoSNOW(cl)
clusterExport(cl, c("met.pheno", "met.data", "met.annot", "y", "inner.cv", "design.feature_selection"))
clusterEvalQ(cl, library(nordlund2014))
dir.create("tuning", showWarnings=FALSE)


#===============================================================================
#   Calculate
#-------------------------------------------------------------------------------

foreach(i=seq_along(cv)) %dopar% {
    outfile <- sprintf("tuning/fold_%i.Rdata", i)
    if(file.exists(outfile)){
        load(outfile)
        if(!any(sapply(fold.feat.sel, is.blank)))
            return(NULL)
    } else {
        fold.feat.sel <- structure(vector("list", ncol(y)), names=names(y))
    }
    save.workspace <- function()
        save(fold.feat.sel, file=outfile)
 
    sink(sprintf("tuning/fold_%i.log", i))
    for(my.class in names(y)[sapply(fold.feat.sel, is.null)]){
        tryCatch({
            cat("\nSelecting features for", my.class, "\n")
            fold.feat.sel[[my.class]] <- design("feature_selection", met.data, y[[my.class]],
                    chr = if(my.class == "sex") c(1:22, "X") else 1:22, cv = inner.cv[[i]])
            save.workspace()
        }, error=function(err){
            print(err)
        })
    }
    sink()

    return(NULL)
}
stopCluster(cl)


#===============================================================================
#   Assumble the results
#-------------------------------------------------------------------------------

for(i in seq_along(cv)){


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
save.workspace()

