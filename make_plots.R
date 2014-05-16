#!/opt/apps/R/3.0.1/bin/Rscript

#===============================================================================
#   Supplementary figure 3
#-------------------------------------------------------------------------------

library(nordlund2013)
library(zoo)

is.complete <- function(x) all(!is.na(x))
qi <- seq(0, 1, length.out=500)
nk <- function(k, x=c(.025, .975)){
    if(any(is.na(k))){
        rep(NA, length(x))
    } else {
        p <- cumsum(qi^k[2]*(1-qi)^k[1])
        p <- p/tail(p,1)
        approx(p, qi, x)$y
    }
}
class.labels <- sub("^ref", "Ref", sub("^sex", "Sex",
        colnames(probs$separate[[11]][[11]])))

load("results/tuning.Rdata")
error.mean <- sapply(error, apply, 2, mean, na.rm=TRUE)
#error.sd <- sapply(error, apply, 2, sd, na.rm=TRUE)
#appropriate.site.range <- which(diff(apply(error[[1]], 2, is.complete)) != 0) + .5
error.mean.min <- min(error.mean, na.rm=TRUE)
ci <- lapply(n.error, apply, 2:3, nk)
ci.mean <- lapply(ci, apply, c(1,3), mean, na.rm=TRUE)

sens <- lapply(conf.tab, function(x)
    apply(x[1,,,,], 2:4, function(xx) prop.table(xx)[1]))
spec <- lapply(conf.tab, function(x)
    apply(x[2,,,,], 2:4, function(xx) prop.table(xx)[2]))
msens <- lapply(sens, apply, c(1,3), mean, na.rm=TRUE)
mspec <- lapply(spec, apply, c(1,3), mean, na.rm=TRUE)
sens.lim <- range(unlist(msens))
spec.lim <- range(unlist(mspec))


true.call.frac <- lapply(conf.tab, function(x)
    apply(x[,1,,,], 2:4, function(xx) prop.table(xx)[1]))
true.call.frac.mean <- lapply(true.call.frac, apply, c(1,3), mean, na.rm=TRUE)
call.lim <- range(unlist(true.call.frac.mean))


X11(,14/cm(1),18/cm(1))
tryCatch({
    pdf("results/S3_tuning_spec.pdf", 14/cm(1), 18/cm(1))
    m <- .5
    pars <- list(ps=8, tcl=-.3, mar=c(3,3,m,m), cex=1)
    pal <- c("blue", "red")
    xd <- .7
    # The lower panels' relative widths
    xf <- c(0, 1.5, 1,1,1,1)
    # The lower panels' relative left/right coordinates in the figure
    xf <- cumsum(xf)
    xf <- xf/tail(xf, 1)
    # The lower panels' relative heights
    yf <- c(0, 1.5, 1, 1.5, 1)
    # All panels' relative bottom/top coordinates in the figure
    yf <- cumsum(yf)
    yf <- c(yf/tail(yf, 1)*2/3, 1)
    
    screens <- split.screen(cbind(
         left = c(0, xd, xf[rep(1:5, 4)]),
         right = c(xd, 1, xf[rep(2:6, 4)]),
         bottom = rep(yf[5:1], c(2,5,5,5,5)),
         top = rep(yf[6:2], c(2,5,5,5,5))))

    screen(screens[1])
    do.call(par, pars)
    plot(c(0, 100), rep(error.mean.min, 2), type="l", col="#e6e6e6", axes=FALSE, ann=FALSE, bty="n",
         xlim=c(1,25), ylim=range(c(unlist(error.mean), unlist(ci.mean)), na.rm=TRUE))
    #vlines(appropriate.site.range, col="#e6e6e6", lty=2)
    matplot(error.mean, type="l", lty=1, col=pal, add=TRUE)
    for(i in 1:2) matplot(t(ci.mean[[i]]), type="l", lty=2, col=pal[i], add=TRUE)

    #plot(c(0, 100), rep(error.mean.min, 2), type="l", col="#e6e6e6", axes=FALSE, ann=FALSE, bty="n",
    #     xlim=c(1,25), ylim=range(c(unlist(error.mean), unlist(error.quantile)), na.rm=TRUE))
    #matplot(error.mean, type="l", lty=1, col=pal, add=TRUE)
    #for(i in 1:2) matplot(t(error.quantile[[i]]), type="l", lty=2, col=pal[i], add=TRUE)
    nice.axis(1, at=c(1,1:5*5), mgp=c(2,.4,0))
    nice.axis(1, at=setdiff(2:24, 1:4*5), labels=rep("", 19), mgp=c(2,.4,0), tck=-0.02)
    nice.axis(2, mgp=c(2,.6,0))
    nice.box()
    mtext(expression(tau), 1, 1.3)
    mtext("Error rate", 2, 2.2)

    screen(screens[2])
    par(mar=c(m,m,m,m))
    blank.plot()
    legend(-1.4, 1.4, c(
        "Consensus CpG sites\nused together\n(combined approach)\n",
        "Subtype CpG sites\nused separately\n(separate approach)\n"),
        lty=c(1,1), col=rev(pal), xpd=TRUE, bty="n", adj=c(0, .94))
    legend(-1.4, 0, c("Mean\n", "95% credibile\ninterval"),
        lty=1:2, col=c("black", "black"), xpd=TRUE, bty="n", adj=c(0, .75))

    pars$mar[3] <- 1
    # Sensitivity
    for(i in 1:10){
        screen(screens[i+2])
        pars$mar[1] <- if(i < 6) m else 2.5
        pars$mar[2] <- if(i %% 5 == 1) 3 else m
        do.call(par, pars)
        blank.plot(c(1,25), c(0,0), ylim=sens.lim)
        hlines(2:5/5, col="#e6e6e6")
        matplot(sapply(msens, "[", i, T), type="l", col=pal, lty=1, add=TRUE)
        if(i > 5){
            nice.axis(1, at=c(1, 1:5*5), mgp=c(2,.4,0))
            mtext(expression(tau), 1, 1.2)
        }
        if((i-1) %% 5 == 0){
            nice.axis(2, mgp=c(2,.6,0), at=2:5/5)
            if(i == 1) mtext("Sensitivity", 2, 2.2, at=fig.usr()[1], xpd=TRUE)
        }
        nice.box()
        mtext(class.labels[i], 3, .1)
    }
    # Specificity
    for(i in 1:10){
        screen(screens[i+12])
        pars$mar[1] <- if(i < 6) m else 2.5
        pars$mar[2] <- if(i %% 5 == 1) 3 else m
        do.call(par, pars)
        blank.plot(c(1,25), c(0,0), ylim=spec.lim)
        hlines(46:50/50, col="#e6e6e6")
        matplot(sapply(mspec, "[", i, T), type="l", col=pal, lty=1, add=TRUE)
        if(i > 5){
            nice.axis(1, at=c(1, 1:5*5), mgp=c(2,.4,0))
            mtext(expression(tau), 1, 1.2)
        }
        if((i-1) %% 5 == 0){
            nice.axis(2, mgp=c(2,.6,0), at=46:50/50)
            if(i == 1) mtext("Specificity", 2, 2.2, at=fig.usr()[1], xpd=TRUE)
        }
        nice.box()
        mtext(class.labels[i], 3, .1)
    }
    close.screen(all=TRUE)
},
error = function(...) cat("Oh no!"),
finally = {
    close.screen(all=TRUE)
    dev.off()
})


#===============================================================================
#   Subtype-specific error rates
#-------------------------------------------------------------------------------

sub.mean <- lapply(sub.error, apply, c(1,3), mean, na.rm=TRUE)
null.err <- sapply(y, function(x) min(prop.table(table(x))))

pdf("results/subtype_specific_error.pdf", 14/cm(1), 6/cm(1))
layout(matrix(1:10, 2, byrow=TRUE))
m <- .5
par(ps=8, tcl=-.3, mar=c(m,m,1,m), oma=c(2,2.5,0,0), cex=1)
for(i in 1:10){
    blank.plot(c(1,25), c(0,0), ylim=c(0, null.err[i]*1.1))
    hlines(null.err[i], col="#cccccc", lty=1)
    my.err <- sapply(sub.mean, "[", i, T)
    matplot(my.err, type="l", col=pal, lty=1, add=TRUE)
    points(apply(my.err, 2, which.min), apply(my.err, 2, min), col=pal)
    if(i > 5){
        nice.axis(1, at=c(1, 1:5*5), mgp=c(2,.4,0))
        mtext("F", 1, 1.5)
    }
    if((i-1) %% 5 == 0){
        nice.axis(2, mgp=c(2,.6,0))
        mtext("Error rate", 2, 2.2)
    }
    nice.box()
    mtext(class.labels[i], 3, .1)
}
dev.off()


#===============================================================================
#   Plot validation set classifications
#-------------------------------------------------------------------------------

library(lattice)
library(gtools)
load("results/pred.Rdata")

pal <- rev(c(Reference="black", `T-ALL`="#d70000", HeH="#00ad23",
             `t(12;21)`="#0073ff", `11q23/MLL`="#85624f", `t(1;19)`="#52e0e0",
             `dic(9;20)`="#ff7182", `t(9;22)`="#003aae", iAMP21="#810086"))
plot.data <- data.frame(
    Sample=factor(val.pred$ID, levels=mixedsort(val.pred$ID)),
    stack(val.pred[9:18]))
names(plot.data)[2:3] <- c("Probability", "Class")
plot.data$Class <- factor(as.character(plot.data$Class), levels=rev(names(y)),
                        labels=c("Sex", names(pal)))
plot.data <- plot.data[order(plot.data$Sample),]
sexes <- with(plot.data, Probability[Class == "Sex"])
sexes <- sprintf("%s %3.0f%%", ifelse(sexes > .5, "Female", "Male"),
                 100*ifelse(sexes > .5, sexes, 1 - sexes))

pdf("results/validation.pdf", 8, 9)
trellis.par.set(list(fontsize = list(text = 8)))
counter <- 0
print(barchart(Class ~ Probability | Sample, plot.data, subset = Class != "Sex",
    col=rev(pal), as.table=TRUE, scales=list(tck=.5),
    panel=function(...){
        counter <<- counter + 1
        panel.segments(.5, -100, .5, 1000, col="#aaaaaa")
        panel.barchart(...)
        panel.text(1.04, 1, sexes[counter], adj=c(1,.5))
    }))
dev.off()


#===============================================================================
#   Plot auto-predictions of the final classifier
#-------------------------------------------------------------------------------

x <- cons.pred[11:20]
x[cons.pred$subtype %in% "reference", 2:9] <- NA

pdf("results/auto_predictions.pdf", 8/cm(1), 8/cm(1))
par(mar=c(2, 3.5, .5, .5), ps=8)
blank.plot(ylim=c(.5,10.5), xlim=0:1)
segments(.5, -10, .5, 100, col="#cccccc")
points(unlist(x),
       jitter(rep(10:1, each=nrow(cons.pred)), amount=.25),
       cex=.75, col=c("#ff000033", "#00000033")[sapply(y, as.integer)])
mtext(c("Reference", names(y)[2:9], "Female"), 2, .3, at=10:1, las=1)
nice.axis(1)
nice.box()
dev.off()

