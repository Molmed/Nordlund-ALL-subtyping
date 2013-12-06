#!/opt/apps/R/3.0.1/bin/Rscript

#===============================================================================
#
#   ooooo      ooo                 .               .o. 
#   `888b.     `8'               .o8               888 
#    8 `88b.    8    .ooooo.   .o888oo   .ooooo.   888 
#    8   `88b.  8   d88' `88b    888    d88' `88b  Y8P 
#    8     `88b.8   888   888    888    888ooo888  `8' 
#    8       `888   888   888    888 .  888    .o  .o. 
#   o8o        `8   `Y8bod8P'    "888"  `Y8bod8P'  Y8P 
#
#
#   This file is sourced by the other analysis scripts and should not be run
#   directly by the user.
#
#   The files should be run in the order:
#    - setup.R
#    - analyze_tune.R (preferably on a computation cluster or multicore machine)
#    - analyze_final.R
#
#===============================================================================
#   Load data and packages
#-------------------------------------------------------------------------------

tryCatch({
    library(pamr)
    library(predict)
    library(nordlund2013)

    load("data/methylation.Rdata")
    load("data/annotations.Rdata")
    load("data/phenotypes.Rdata")

    # Remove sample not used in the study
    # and sites that did not pass quality control
    sample.idx <- !is.na(met.pheno$subtype) & !grepl("rep[2-9]$", met.pheno$id)
    site.idx <- met.annot$analyzed
        # ^^^ needs to be saved as a separate variable since we need it later
    met.pheno <- met.pheno[sample.idx,]
    met.data <- t(met.data[site.idx, sample.idx])
    met.annot <- met.annot[site.idx,]
}, error=function(...)
    cat("Could not setup analysis environment. Did you run `setup.R` first?\n"))


#===============================================================================
#   Prepare responses and output variables
#-------------------------------------------------------------------------------

dir.create("results", showWarnings=FALSE)
if(file.exists("results/pred.Rdata")){
    load("results/pred.Rdata")
} else {
    attach(met.pheno)


    # Make the response vectors

    class.types <- c("T-ALL", "HeH", "t(12;21)", "11q23/MLL", "t(1;19)", "dic(9;20)", "t(9;22)", "iAMP21")
    known.types <- c(class.types, "<45chr", ">67chr", "biclone")
    unknown.types <- c("other", "normal", "no result")
    y <- data.frame(
        reference = factor(ifelse(subtype == "reference", 1,
                           ifelse(subtype %in% known.types, 2, NA)),
                     labels=c("reference", "ALL")),
        lapply(class.types, function(lev){
            factor(ifelse(subtype %in% known.types, ifelse(subtype == lev, 1, 2), NA),
                   labels=c(lev, paste("not", lev)))
        }),
        sex = factor(ifelse(sample.type %in% c("diagnosis", "remission"),
                            sex, NA),
                     levels=1:2, labels=levels(sex))
    )
    names(y)[1+seq(class.types)] <- class.types


    # Setup cross validation

    set.seed(123)
    tmp <- factor(paste(subtype, sex))
    cv <- resample.crossval(tmp, nfold=5, nrep=5)
    inner.cv <- lapply(cv, function(idx) resample.crossval(tmp, 5, 5, subset=!idx))
    rm(tmp)
    feat.sel <- cons <- structure(vector("list", ncol(y)), names=names(y))
    
    save.workspace <- function(){
        obj <- c("sample.idx", "y", "cv", "inner.cv",
            "class.types", "known.types", "unknown.types",
            "feat.sel", "cons.sites", "cons", "pred", "cons.pred",
            "val.pred", "save.workspace")
        save(list=obj[sapply(obj, exists)], file="results/pred.Rdata")
    }
    save.workspace()
    detach(met.pheno)
}


#===============================================================================
#   Feature selection plugin for predict::batch.predict()
#-------------------------------------------------------------------------------

design.feature_selection <- function(x, y, chr, cv){
    if(missing(cv)){
        cv <- resample.crossval(y, 5, 5)
    } else {
        cv[is.na(y),] <- NA
    }
    my.pre.trans <- function(...){
        sets <- pre.trans.450k(..., feat=met.annot$CHR %in% chr, dbeta=.2, na.frac=.1)
        predict::trace.msg(3, "%i features passed filter.", ncol(sets$design), time=FALSE)
        sets
    }

    my.pred <- batch.predict(x, y, "nsc", test.subset=cv,
        pre.trans=my.pre.trans, save.fits=TRUE, save.vimp=TRUE, .verbose=TRUE)
    # Return an integer vector of how many times each feature had a coefficinent != 0
    apply(sapply(subtree(my.pred$cv, T, "nsc", "vimp", flatten=2),
                         function(x) apply(abs(x), 1, sum) > 1e-12), 
          1, sum, na.rm=TRUE)
    # 1e-12 for avoiding comparison with exactly 0.
    # Might not be an issue, but doesn't harm.
}

