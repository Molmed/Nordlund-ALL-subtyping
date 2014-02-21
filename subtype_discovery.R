library(predict)
library(survival)
library(parallel)
workers <- makeCluster(16)
clusterEvalQ(workers, {
    library(survival)
    library(predict)
})
load("data/annotations.Rdata")
attach(met.annot)

load("/proj/b2010028/private/nobackup/jess/subtyping_paper/reclassified_type_2014-01-09.Rdata")
# We need our dataset to get the survival data
library(analyse450k)
load.450k.data("all")
sample.idx <- with(all.pheno,
    !is.na(months.to.primary.event) &
    primary.event %in% c("CR1", "relapse") &
    public.id %in% sub.pheno[
        sub.pheno[,"sub.pheno"] %in% c("Multi Class", "No Class"),
        "public.id"
    ]
)

all.met <- t(unname(all.met[,sample.idx]))
attach(all.pheno[sample.idx,])
yo <- outcome(months.to.primary.event, primary.event, "relapse", "CR1")
y <- as.Surv(yo)

p.value.coxph <- function(x, log=FALSE, test=c("logrank", "wald", "likelihood"), ...){
    test <- match.arg(test)
    df <- sum(!is.na(x$coefficients))
    switch(test,
        logrank = pchisq(x$score, df, lower.tail=FALSE, log=log),
        wald = pchisq(x$score, df, lower.tail=FALSE, log=log),
        likelihood = pchisq(x$logtest, df, lower.tail=FALSE, log=log))
}

clusterExport(workers, "p.value.coxph")
my.y <- "Grrrrr!"
my.strat <- "iwiwiwi"
fun <- function(idx, strat=FALSE){
    cat(".")
    my.y <<- y[idx]            # Hack to get around the fact that
    if(strat){
        my.strat <<- as.integer(risk.group[idx]) - 2
        clusterExport(workers, c("my.y", "my.strat"))   # clusterExport only looks in globalenv
        parCapply(workers, all.met[idx,], function(x) p.value(coxph(my.y ~ x + strata(my.strat))))
    } else {
        clusterExport(workers, "my.y")
        parCapply(workers, all.met[idx,], function(x) p.value(coxph(my.y ~ x)))
    }
}
pval <- data.frame(all = fun(TRUE, TRUE),
                   IR = fun(risk.group %in% "IR"),
                   HR = fun(risk.group %in% "HR"))
save(pval, sample.idx, file="pval.Rdata")


# Summary
site.ind <- lapply(pval, function(p){
    x <- which(p.adjust(p, "fdr") < .05)
    lapply(split(x, CHR[x]), function(x) x[order(MAPINFO[x])])
})
min.dist <- sapply(site.ind, sapply, function(x)
    if(length(x) > 1) min(dist(matrix(x))) else NA)

bh <- lapply(pval, p.adjust)
lapply(bh, function(x) table(x < .01))
{
pdf("discovery_plots/clustering.pdf", 8, 12)
par(mar=c(.2, .2, 1, 5), oma=c(4,4,2,0))
layout(matrix(1:2), heights=c(1,3))
for(i in seq_along(site.ind)){
    surv.met <- all.met[,unlist(site.ind[[i]])]
    cl <- hclust(dist(surv.met))
    clust.plot(cl)
    met.image(t(surv.met[cl$order,]))
    color.cluster.leaves(factor.events(yo)[cl$order])
    #mtext(sprintf("%i: %i", CHR[site.ind], MAPINFO[site.ind]),
    #      4, at=seq_along(site.ind), las=1, cex=.7)
}
dev.off()
}


# QQ-plot
pval.o <- apply(pval, 2, function(x) -sort(x))
pval.e <- -log(1:ncol(all.met)/ncol(all.met))
ind <- sort(sample(length(pval.e), 1000))
pdf("discovery_plots/qq_plot.pdf")
matplot(pval.e[ind], pval.o[ind,], type="l", col="#0000aa33", lty=1)
segments(0,0,100,100)
dev.off()

    
# Mean SD histogram
plot.data <- cbind(mean = parRapply(cl, pval, mean),
                   sd = parRapply(cl, pval, sd))
h <- hist_2d(plot.data, bins=c(100,100))
pdf("discovery_plots/mean_sd_histogram.pdf")
image(h, log=10)
dev.off()


# Clustering

thres <- sort(plot.data[,"mean"])[1000]
ranks <- rank(plot.data[,"mean"])
site.idx <- 10000 < ranks & ranks < 11000

surv.met <- all.met[,site.idx]

cl <- list(hclust(dist(surv.met)), hclust(dist(t(surv.met))))
pdf("discovery_plots/clustering.pdf")
met.image(t(surv.met[cl[[1]]$order, cl[[2]]$order]))
dev.off()

table(cutree(cl[[1]], 15))
eval.survival(as.outcome(y), cutree(cl[[1]], 15) > 1)



eval.survival <- function (y, groups, time.point = 60) {
    library(cmprsk)
    subset <- !is.na(y) & !is.na(groups)
    sf <- survival::survfit(as.Surv(y) ~ groups)
    ssf <- summary(sf)
    surv <- sapply(split(data.frame(time = ssf$time, surv = ssf$surv), 
        ssf$strata), function(dsf) {
        if (nrow(dsf) == 0 || dsf$time[1] > time.point) 
            return(1)
        dsf$surv[max(which(dsf$time <= time.point))]
    })
    n <- table(groups[subset])
    n.relapse <- sapply(split(integer.events(y)[subset] == 1, 
        groups[subset]), sum, na.rm = TRUE)
    pval <- p.value(cuminc(y$time, integer.events(y), groups))
    out <- list(fit = sf, surv = surv, n = n, n.relapse = n.relapse, 
        p.value = pval)
    class(out) <- "eval_survival"
    return(out)
}


#----------------------------------------------------------------------------------

cv <- resample.crossval(factor.events(yo), 2, 10, subset=risk.group %in% "HR")

my.y <- "Grrrrr!"
fun <- function(fold){
    cat(".")
    idx <- na.fill(!fold, FALSE)
    my.y <<- y[idx]           # Hack to get around the fact that
    clusterExport(workers, "my.y")   # clusterExport only looks in globalenv
    parCapply(workers, all.met[idx,], function(x) p.value(coxph(my.y ~ x), log=TRUE))
}
cv.pval <- vapply(cv, fun, numeric(ncol(all.met)))
save(pval, sample.idx, cv, file="cv_pval.Rdata")


# Plot
h <- vector("list", attr(cv, "nrep"))
pp <- matrix(NA, ncol(all.met), attr(cv, "nrep"))
for(i in seq_len(attr(cv, "nrep"))){
    pp[,i] <- pval[,i*2-1] * pval[,i*2]
    h[[i]] <- hist_2d(pval[,i*2 - 1:0], bins=c(100, 100))
}
met.range <- apply(all.met, 2, function(x) diff(range(x, na.rm=TRUE)))
pp[met.range < .2,] <- NA

o <- order(apply(pp, 1, min))
site.ind <- o[1:10]
site.ind <- order(bh[,1])[1:30]
surv.met <- all.met[,site.ind]

cl <- list(hclust(dist(surv.met)), hclust(dist(t(surv.met))))
pdf("discovery_plots/clustering.pdf")
par(mar=c(.2, .2, 1, 5), oma=c(4,4,2,0))
layout(matrix(1:2), heights=1:2)
clust.plot(cl[[1]])
met.image(t(surv.met[cl[[1]]$order, cl[[2]]$order]))
color.cluster.leaves(factor.events(y)[cl[[1]]$order])
mtext(sprintf("%i: %i", CHR[site.ind][cl[[2]]$order],
              MAPINFO[site.ind][cl[[2]]$order]),
      4, at=seq_along(site.ind), las=1, cex=.8)
dev.off()

eval.survival(y, cutree(cl[[1]], 3))

