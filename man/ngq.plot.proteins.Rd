\name{ngq.plot.proteins}
\alias{ngq.plot.proteins}
\title{Cumulative density plot of protein logratios.}
\usage{
  ngq.plot.proteins(proteins, type = "cdf", colors = NULL,
    pch = NA, bty = "l", median = FALSE,
    annotation.cex = NA, ...)
}
\arguments{
  \item{proteins}{A \code{proteins} data frame.}

  \item{type}{The plot type: "cdf" (the default) or
  "volcano". See details.}

  \item{colors}{A color vector; one color for every plex.}

  \item{pch}{The point style. Defaults to \code{NA}.}

  \item{bty}{The box style. Who wants complete boxes
  anyway?}

  \item{median}{Boolean indicator if the median fold change
  should be marked by a vertical line in the CDF plot.
  Default is \code{FALSE}.}

  \item{annotation.cex}{Relative size of annotation text}

  \item{...}{Additional graphics parameters that are passed
  on to \code{plot}/\code{points}.}
}
\description{
  This provides plotting functions for quantitative
  protein-level NGQ measurements.
}
\section{CDF plots}{
  This is a simple empirical cumulative density plot. In
  order to center the plot, the user must provide suitable
  \code{xlim} parameters (see \code{?par}). In many cases
  \code{lwd=2} seems to be a good choice.
}

\section{Volcano plot}{
  The volcano plot illustrates the results of the p-value
  estimation carried out in \code{ngq.significance}. The
  x-axis shows (log2-based) logratio changes, the y-axis
  corresponds to \code{-log10(p)}. The median of all
  logratio changes is set to zero (this corresponds to the
  test assumptions).
}
\examples{
\dontrun{
    data(xilac_peptides)
    p <- ngq.peptides(xilac_peptides)
    q <- ngq.proteins(p)
    s <- ngq.significance(p, q, p.adjust.method="BH")
    par(mfrow=c(1,2))
    # plot CDF
    ngq.plot.proteins(s, type="cdf", median=TRUE)
    # and a volcano plot
    ngq.plot.proteins(s, type="volcano")
}
}

