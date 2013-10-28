#!/opt/apps/R/3.0.1/bin/Rscript

#===============================================================================
#   Parameter tuning and performance estimation
#-------------------------------------------------------------------------------
#
#   This script supports two forms of parallelizations:
#
#       - Multithreading within a single CPU with multiple cores
#       - Manual batch execution in parallel on multiple machines
#
#   Multithreading peak memory requirement of about 13 GB of RAM/core, so set
#   the variable `number.of.cores` accordingly.
#
#   To use only one machine, simply source this script into an R session, or
#   batch execute it in any of the following ways:
#
#       $ R -f analyze_tune.R
#       $ ./analyze_tune.R --no-save --no-restore
#       $ R CMD BATCH analyze_tune.R --no-save --no-restore
#   
#   To use more than one machine, divide the folds between them manually, and
#   execute the script like this (on the separate machines of course):
#
#       $ R -f analyze_tune.R --args 1 2 3 4 5 6 7 8 9
#       $ R -f analyze_tune.R --args 10 11 12 13 14 15 16 17
#       $ R -f analyze_tune.R --args 18 19 20 21 22 23 24 25
#
#   The above will create three separate processes that all multithread
#   internally. If the number of folds given exceeds the designated number of
#   cores they will be queued.
#
#   Should the script be terminated unexpectedly, it will resume from where it
#   stopped when restarted.
#
#-------------------------------------------------------------------------------

number.of.cores <- 3

# If we run on the UPPMAX cluster, maximize the number of processes
if(grepl("^q\\d+\\.uppmax\\.uu\\.se$", Sys.info()["nodename"])){
    max.mem <- as.integer(sub("^MemTotal:\\s+(\\d+) kB$", "\\1",
        system("head /proc/meminfo -n 1", intern=TRUE)))
    number.of.cores <- floor(max.mem/13e6)
}

a <- commandArgs()
if("--args" %in% a){
    my.folds <- as.integer(a[(which(a == "--args") + 1):length(a)])
} else {
    my.folds <- seq_along(cv)
}
number.of.cores <- min(number.of.cores, length(my.folds))

source("analyze_init.R")
dir.create("tuning", showWarnings=FALSE)

library(doMC)
registerDoMC(number.of.cores)


#-------------------------------o
#   Calculate

foreach(i=my.folds) %dopar% {
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

