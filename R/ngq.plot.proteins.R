#' Cumulative density plot of protein logratios.
#' 
#' This provides plotting functions
#' for quantitative protein-level NGQ measurements. 
#' 
#' @export
#' @param proteins A \code{proteins} data frame.
#' @param type The plot type: "cdf" (the default) or "volcano". See details.
#' @param colors A color vector; one color for every plex.
#' @param pch The point style. Defaults to \code{NA}.
#' @param bty The box style. Who wants complete boxes anyway?
#' @param median Boolean indicator if the median fold change should be
#'               marked by a vertical line in the CDF plot. Default is
#'               \code{FALSE}.
#' @param annotation.cex Relative size of annotation text 
#' @param ... Additional graphics parameters that are passed on to 
#'            \code{plot}/\code{points}.
#' 
#' @section CDF plots: This is a simple empirical cumulative density plot. In
#' order to center the plot, the user must provide suitable \code{xlim}
#' parameters (see \code{?par}). In many cases \code{lwd=2} seems to be a
#' good choice.
#' 
#' @section Volcano plot: The volcano plot illustrates the results of
#' the p-value estimation carried out in \code{ngq.significance}. The
#' x-axis shows (log2-based) logratio changes, the y-axis corresponds to
#' \code{-log10(p)}. The median of all logratio changes is set to zero
#' (this corresponds to the test assumptions).
#' 
#' @examples \dontrun{
#'     data(xilac_peptides)
#'     p <- ngq.peptides(xilac_peptides)
#'     q <- ngq.proteins(p)
#'     s <- ngq.significance(p, q, p.adjust.method="BH")
#'     par(mfrow=c(1,2))
#'     # plot CDF
#'     ngq.plot.proteins(s, type="cdf", median=TRUE)
#'     # and a volcano plot
#'     ngq.plot.proteins(s, type="volcano")  
#' }
#' 
ngq.plot.proteins <- function(proteins, type="cdf", colors=NULL,
        pch=NA, bty="l", median=FALSE, annotation.cex=NA, ...)
{
    # check the plot type
    if (type == "cdf") {
        cols <- proteins$meta$cols.logratio
        n <- length(cols)
        if (is.null(colors)) {
            colors <- brewer.pal(n, "Set1")
        }
        add <- F
        for (i in 1:n) {
            if (i > 1) {
                add <- T
            }
            d <- proteins$proteins[[cols[i]]]
            plot(ecdf(d), add=add, pch=pch, bty=bty, col=colors[i], ...)
            abline(v=0, col=colors[i], lty=3)
            if (median) {
                abline(v=median(d), col=colors[i], lty=2)
            }
        }
    }
    if (type == "volcano") {
        logratio.cols <- proteins$meta$cols.logratio
        qvalue.cols <- proteins$meta$cols.qvalue
        n <- length(logratio.cols)
        if (is.null(colors)) {
            colors <- brewer.pal(n, "Set1")
        }
        fun <- plot
        for (i in 1:n) {
            if (i > 1) {
                fun <- points
            }
            x <- proteins$proteins[[logratio.cols[i]]]
            x <- x - median(x)
            p <- -log10(proteins$proteins[[qvalue.cols[i]]])
            if (is.na(pch)) {
                pch <- 20
            }
            fun(x, p, pch=pch, bty=bty, col=colors[i], ...)
            abline(v=0, col=colors[i], lty=3)
            # text
            ix <- ((!(abs(x) < 1)) & p > -log10(0.01))
            if (!all(!ix)) {
                text(x[ix], p[ix], labels=proteins$proteins$ProteinGroup[ix],
                        cex=annotation.cex, pos=4)
            }
        }
        abline(v=c(-1,1), col="lightgray", lty=2)
        abline(h=-log10(0.05), col="lightgray", lty=2)
        abline(h=-log10(0.01), col="lightgray", lty=2)
    }
}