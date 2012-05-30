\name{ngq.write.xlsx}
\alias{ngq.write.xlsx}
\title{Store the NGQ protein and peptide quantification in an
MS Excel zipped XML file (.xlsx).}
\usage{
  ngq.write.xlsx(filename, peptides = NULL,
    proteins = NULL)
}
\arguments{
  \item{filename}{The name of the Excel file.}

  \item{peptides}{A \code{peptides} data frame with
  processed (i.e.  quantified) peptide data. See
  \code{ngq.peptides}.}

  \item{proteins}{A \code{proteins} data frame with
  protein-level quantification data.}
}
\description{
  Stores the peptide- and/or protein-level information into
  an MS Excel file. If \code{peptides} and \code{proteins}
  are both specified, this will store the data in two
  different sheets, but in the same file.
}
\section{Details}{
  Either \code{peptides} or \code{proteins} must be
  non-\code{NULL} or the function will fail.
}

\section{Excel format}{
  The function depends on the \code{XLConnect} package in
  order to store data in the Excel format. \code{XLConnect}
  is based on a Java core; should writing the file run out
  of memory, setting suitable parameters for the JVM may
  help (see examples).
}
\examples{
\dontrun{
    data(xilac_peptides)
    p <- ngq.peptides(xilac_peptides)
    # set up the virtual machine; this must be done
    # prior to loading the XLConnect library, i.e. prior
    # to the first call to ngq.write.xlsx
    options(java.parameters = "-Xmx1024m")
    ngq.write.xlsx("somefile.xlsx", peptides=p)
}
}

