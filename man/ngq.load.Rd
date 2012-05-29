\name{ngq.load}
\alias{ngq.load}
\title{Loads a raw NGQ quantification series file.}
\usage{
  ngq.load(filename, labels, min.quality = 0.95,
    na.rm = FALSE, require.single = FALSE,
    exclude = "^RRRRR")
}
\arguments{
  \item{filename}{The filename of the file to load.}

  \item{labels}{The names of the different labels.}

  \item{min.quality}{The minimum quality across all fits
  for a series to be accepted.}

  \item{na.rm}{If TRUE, removes all rows with NAs.}

  \item{require.single}{If TRUE, removes all rows that do
  not stem from a single fit (i.e. that have multiple fit
  IDs).}

  \item{exclude}{Regular expression that specifies which
  protein groups to exclude. Useful for removing reverse DB
  hits. Set to NA to disable.}
}
\value{
  A data frame with elements \item{peptides}{A data.frame
  with all peptide data.} \item{meta}{Metadata used by
  other functions in the package.}
}
\description{
  This function loads all raw peptide data from an NGQ
  result and carries out basic preprocessing steps. The
  latter may include filtering for fit quality, removing
  quantifications with missing values, requiring single-fit
  quantifications for every peptide, and removing
  quantifications that stem from reverse database hits.
}
\section{Details}{
  The function reads the column headers from the NGQ result
  file and extracts the respective columns. It depends on
  the column order used by NGQ/massquito and any change to
  the column style is likely to break the code.
}

\section{Parameters}{
  The function provides reasonable default parameters for
  most application scenarios. Note that excluding peptide
  measurements with missing values in muliplex experiments
  may influence the outcome of a later protein-level
  quantification step, in particular if simplicial
  transformations are used to mitigate NA effects.
}

