
# This file will not accompany the paper!

library(predict)
library(nordlund2013)

load("data/methylation.Rdata")
load("data/annotations.Rdata")
load("data/phenotypes.Rdata")
load("results/pred.Rdata")

attach(met.pheno)
attach(met.annot)
sample.idx <- met.pheno$sample.type %in% c("relapse 1", "relapse 2")
#rel.met <- list(
#    impute.knn(t(met.data[TargetID %in% cons.sites$TargetID & CHR %in% 1:22, sample.idx]),
#               distmat="auto"),
#    impute.knn(t(met.data[TargetID %in% cons.sites$TargetID, sample.idx]),
#               distmat="auto"))
rel.met <- lapply(with(cons.sites, split(TargetID, Subtype)), function(ids)
    impute.knn(t(met.data[met.annot$TargetID %in% ids, sample.idx]), distmat="auto"))
                  
#rel.pred <- mapply(predict, cons, rel.met[c(rep(1,9),2)], SIMPLIFY=FALSE)
rel.pred <- mapply(predict, cons, rel.met, SIMPLIFY=FALSE)
rel.prob <- data.frame(lapply(rel.pred, function(x) x$prob[,1]), check.names=FALSE)
rel2diag <- match(sub("r\\d$", "", id[sample.idx]), id)
rel.tab <- data.frame(
    diagnosis.id = id[rel2diag],
    relapse.id = id[sample.idx],
    subtype = subtype[rel2diag],
    predicted = apply(rel.prob[1:9] >= .5, 1, function(x) paste(names(which(x)), collapse="; ")),
    rel.prob,
    check.names=FALSE)

write.table(rel.tab, "results/relapse_predictions.csv",
            quote=FALSE, sep=",", row.names=FALSE)

