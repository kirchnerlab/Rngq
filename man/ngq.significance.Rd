\name{ngq.significance}
\alias{ngq.significance}
\title{Significance calculation for protein up- or downregulation.}
\usage{
  ngq.significance(peptides, proteins,
    use.test = wilcox.test, p.adjust.method = "none")
}
\arguments{
  \item{peptides}{A \code{peptides} data frame.}

  \item{proteins}{A \code{proteins} data frame.}

  \item{use.test}{The name of the test that should be used
  to calculate p-values (e.g. \code{wilcox.test},
  \code{t.test}).  Defaults to the non-parametric Wilcoxon
  test.}

  \item{p.adjust.method}{The method that should be used to
  correct the p-values for multiple testing. See
  \code{p.adjust} for a list of valid parameters. Note that
  this defaults to \code{none}, i.e. no multiple testing
  correction!}
}
\description{
  The function compares the peptide-level abundance
  distribution of every protein against the overall
  protein-level abundance distribution and determines a
  p-value for the difference between the two distributions.
}
\section{Details}{
  The underlying assumption behind this test setup is that
  the protein-level abundance logratio distributions are
  somewhat robust estimates of the null distribution. If
  the abundance distribution of the peptides of a specific
  protein deviates from this distribution, this is an
  indicator that it is indeed differentially regulated. In
  essence, this means that we rely on the smoothing power
  of the mean/median.
}
\examples{
data(xilac_peptides)
    p <- ngq.peptides(xilac_peptides)
    q <- ngq.proteins(p)
    s <- ngq.significance(p, q, p.adjust.method="BH")
}

