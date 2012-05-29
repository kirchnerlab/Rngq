\name{ngq.proteins}
\alias{ngq.proteins}
\title{Calculate NGQ protein-level fold changes, logratios, deviations and CVs.}
\usage{
  ngq.proteins(x, f = x$peptides$ProteinGroup)
}
\arguments{
  \item{x}{A peptide-quantified ngq.peptide data frame
  (i.e. processed with ngq.peptides().}

  \item{f}{A factor column that defines the protein group
  memberships of the peptides. This defaults to the
  ProteinGroup column in the ngq.peptide data frame.}
}
\value{
  A data frame with all protein-level quantifications
}
\description{
  This function takes an ngq.peptide data frame, derives
  protein-level quantification and variance information and
  returns an ngq.protein data frame.
}

