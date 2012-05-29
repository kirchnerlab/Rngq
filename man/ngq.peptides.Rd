\name{ngq.peptides}
\alias{ngq.peptides}
\title{Calculate NGQ peptide-level (log-)ratios.}
\usage{
  ngq.peptides(x, ref = x$meta$labels[1])
}
\arguments{
  \item{x}{A raw ngq.peptides data frame.}

  \item{ref}{The name of the reference label whose value is
  used as the denominator in all fold change and logratio
  calculation steps. The value this parameter must be among
  x$meta$labels.}
}
\value{
  The ngq.peptides data frame with fold change and logratio
  columns appended.
}
\description{
  The function takes a raw ngq.peptides data frame and
  appends the fold change and logratio (base 2) columns.
}

