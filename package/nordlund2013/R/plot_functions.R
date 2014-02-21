##' Initiate a plot without drawing anything
##' 
##' @param x Mainly for not having to set \code{xlim} explicitly.
##' @param y Mainly for not having to set \code{ylim} explicitly.
##' @param ... Sent to \code{\link{plot}}.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
blank.plot <- function(x=0, y=0, ...)
    plot(x, y, ..., type="n", ann=FALSE, axes=FALSE, xaxs="i", yaxs="i")


##' Add vertical or horizontal lines to a plot
##'
##' @param x Coordinates of vertical lines.
##' @param lend Line ending style, see \code{\link{par}}.
##' @param ... Sent to \code{\link{segments}}.
##' @examples
##' plot(0:10, 0:10, type="n")
##' hlines(0:4*2.5, col="#dddddd")
##' points(0:10, 0:10)
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
vlines <- function(x, lend=1, ...)
    segments(x, par("usr")[3], x, par("usr")[4], lend=lend, ...)
##' @param y Coordinates of horizontal lines.
##' @rdname vlines
##' @export
hlines <- function(y, lend=1, ...)
    segments(par("usr")[1], y, par("usr")[2], y, lend=lend, ...)


##' Plots an axis the way an axis should be plotted.
##'
##' @param ... Sent to \code{\link{axis}}.
##' @param las Rotation of axis labels. Always horizontal by default.
##' @param lwd Width of the line drawn along the plot area. Omitted by default
##'   since it overlaps with \code{\link{box}} and causes it to look thicker
##'   where the axis is.
##' @param lwd.ticks Width of the tick lines. These are kept by default.
##' @param lend Line endings, see \code{\link{par}}.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
nice.axis <- function(..., las=1, lwd=0, lwd.ticks=par("lwd"), lend=2){
    axis(..., las=las, lwd=0, lwd.ticks=lwd.ticks, lend=1)
    if(lwd){
        # Plot the line along the axis
        args <- c(list(...), list(lwd=lwd, lwd.ticks=0, lend=lend))
        args$labels <- FALSE
        do.call(axis, args)
    }
}
##' Plots a box around a plot
##'
##' @param lend Line ending style, see \code{\link{par}}. Defaults to square.
##' @param ljoin Line joint style, see \code{\link{par}}. Defaults to mitre,
##'   i.e. 90 degree corners in this case.
##' @param ... Sent to \code{\link{box}}.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
nice.box <- function(lend=2, ljoin=1, ...) box(lend=lend, ljoin=ljoin, ...)


