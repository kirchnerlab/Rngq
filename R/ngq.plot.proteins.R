#' Cumulative density plot of protein logratios.
#' 
#' This provides a simple an straightforward plotting function
#' for the quantitative protein-level NGQ measurements. 
#' 
ngq.plot.proteins <- function(proteins, type="cdf", colors=NULL,
        pch=NA, bty="l", median=FALSE, ...)
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
            text(x[ix], p[ix], labels=proteins$proteins$GeneName[ix])
        }
        abline(v=c(-1,1), col="lightgray", lty=2)
        abline(h=-log10(0.05), col="lightgray", lty=2)
        abline(h=-log10(0.01), col="lightgray", lty=2)
    }
}