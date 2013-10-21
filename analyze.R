#
#   This script performs the analyses that the paper is based on.
#   
#   Since it takes quite some time to carry out all of it, it is written in a
#   way that allows halting and resuming without substantial loss of results.
#
#   
#===============================================================================
#   Load data and packages
#-------------------------------------------------------------------------------

tryCatch({
    setwd("data")
    if(!exists("met.data"))  load("methylation.Rdata")
    if(!exists("met.annot")) load("annotations.Rdata")
    if(!exists("met.pheno")) load("phenotypes.Rdata")
    setwd("../results")

    library(pamr)
    library(predict)
    library(nordlund2014)
}, error=function(...)
    cat("Could not setup analysis environment. Did you run `setup.R` first?\n"))

# Remove sample not used in the study
sample.idx <- !is.na(met.pheno$subtype) & !grepl("rep[2-9]$", met.pheno$id)
met.pheno <- met.pheno[sample.idx,]
#met.data <- t(met.data[,sample.idx])

   # TEMP
   site.ind <- 1:1000
   met.data <- t(met.data[site.ind, sample.idx])
   met.annot <- met.annot[site.ind,]


#===============================================================================
#   Initialize
#-------------------------------------------------------------------------------

if(file.exists("pred.Rdata")){
    load("pred.Rdata")
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

    feat.sel <- cons <- structure(vector("list", ncol(y)), names=names(y))
    pred <- NULL
    
    save.workspace <- function(){
        save(sample.idx, y, cv, class.types, known.types, unknown.types,
             feat.sel, cons, pred, file="pred.Rdata")
    }
    save.workspace()
    detach(met.pheno)
}


#===============================================================================
#   Perform classification
#-------------------------------------------------------------------------------

design.feature_selection <- function(x, y, chr, cv){
    if(missing(cv)){
        cv <- resample.crossval(y, 5, 5)
    } else {
        cv[is.na(y),] <- NA
    }
    my.pre.trans <- function(...)
        pre.trans.450k(..., feat=met.annot$CHR %in% chr, dbeta=.2, na.frac=.1)
    my.pred <- batch.predict(x, y, "nsc", test.subset=cv, error.fun=error.rate,
        pre.trans=my.pre.trans, save.fits=FALSE, save.vimp=TRUE, .verbose=TRUE)
    # Return an integer vector of how many times each feature had a coefficinent != 0
    apply(sapply(subtree(my.pred$cv, T, "nsc", "vimp", flatten=2),
                         function(x) apply(abs(x), 1, sum) > 1e-12), 
          1, sum, na.rm=TRUE)
    # 1e-12 for avoiding comparison with exactly 0.
    # Might not be an issue, but doesn't harm.
}


# Select sites
for(my.class in names(y)[sapply(feat.sel, is.null)]){
    cat("\nSelecting features for", my.class, "\n")
    feat.sel[[my.class]] <- design("feature_selection", met.data, y[[my.class]],
                    chr = if(my.class == "sex") c(1:22, "X") else 1:22, cv = cv)
    save.workspace()
}
cons.sites <- do.call(rbind, lapply(names(y), function(my.class){
    data.frame(Subtype = my.class,
               TargetID = met.annot$TargetID[feat.sel[[my.class]] >= 10],
               stringsAsFactors=FALSE)
}))
write.table(cons.sites, "consensus_sites.csv", quote=FALSE, sep="\t",
            row.names=FALSE)


# Train consensus classifiers
cons.met <- impute.knn(met.data[,met.annot$TargetID %in% cons.sites$TargetID],
                       distmat="auto")
for(my.class in names(y)[sapply(cons, is.null)]){
    cat("Making final classifier for", my.class, "\n")
    idx <- !is.na(y[[my.class]])
    cons[[my.class]] <- design("nsc", cons.met[idx,], y[[my.class]][idx])
}


# Predict class probabilities of all samples
# (including the ones used for training)
pred <- lapply(cons, predict, cons.met)
save.workspace()
cons.pred <- data.frame(met.pheno,
    lapply(pred, function(p) p$prob[,1]))
names(cons.pred)[ncol(cons.pred)] <- "sex.female"
write.table(cons.pred, "consensus_predictions.csv",
    quote=FALSE, sep="\t", row.names=FALSE)


#===============================================================================
#   Estimate performance and tune with an additional layer of CV
#-------------------------------------------------------------------------------



