#' Loads a raw NGQ quantification series file.
#' 
#' This function loads all raw peptide data from an NGQ result and carries out
#' basic preprocessing steps. The latter may include filtering for fit
#' quality, removing quantifications with missing values, requiring single-fit
#' quantifications for every peptide, and removing quantifications that stem 
#' from reverse database hits.
#' 
#' @section Details: The function reads the column headers from
#' the NGQ result file and extracts the respective columns. It depends on the
#' column order used by NGQ/massquito and any change to the column style is
#' likely to break the code.
#' 
#' @section Parameters: The function provides reasonable default parameters for 
#' most application scenarios. Note that excluding peptide measurements with
#' missing values in muliplex experiments may influence the outcome of
#' a later protein-level quantification step, in particular if simplicial
#' transformations are used to mitigate NA effects.
#'
#' @export
#' @param filename The filename of the file to load.
#' @param labels The names of the different labels.
#' @param min.quality The minimum quality across all fits for a series to be
#'                   accepted.
#' @param na.rm If TRUE, removes all rows with NAs.
#' @param require.single If TRUE, removes all rows that do not stem from a 
#'                         single fit (i.e. that have multiple fit IDs).
#' @param exclude Regular expression that specifies which protein groups to
#'                exclude. Useful for removing reverse DB hits. Set to NA to
#'                disable.
#' @return A data frame with elements
#'         \item{peptides}{A data.frame with all peptide data.}
#'         \item{meta}{Metadata used by other functions in the package.}
ngq.load <- function(filename, labels, min.quality=0.95, na.rm=FALSE, 
        require.single=FALSE, exclude="^RRRRR") 
{
    # get the column names
    col.names <- strsplit(readLines(filename, n=1), '\t')[[1]]
    # determine the plex number
    plex <- length(grep("FitId", col.names))
    # make sure the number of labels fits with the plex
    stopifnot(plex == length(labels))
    # indexes for specific columns
    abundance.cols <- cumsum(c(10, rep(3, plex-1)))
    quality.cols <- abundance.cols + 1 
    fit.cols <- abundance.cols + 2
    # fix up the column names
    col.names[abundance.cols] <- labels
    col.names[quality.cols] <- paste(labels, "(q)", sep="")
    col.names[fit.cols] <- paste(labels, "(id)", sep="")
    # define the column classes
    col.classes <- c(rep("character",5), rep("numeric", 4 + 3*plex))
    # read the data
    x <- read.delim(filename, header=T, 
            colClasses=col.classes, sep="\t", quote="")
    names(x) <- col.names
    # remove all NAs if requested
    if (na.rm) {
        z <- apply(is.na(x[,abundance.cols]), 1, sum)
        x <- x[z==0,]
    }
    # require a minimum fit quality
    low.q <- apply(x[,quality.cols] < min.quality, 1, sum, na.rm=T)
    x <- x[low.q == 0,]
    # exclude protein groups that follow a specified pattern
    if (!is.na(exclude)) {
        x <- x[grep(exclude, x$ProteinGroup, invert=TRUE),]
    }
    # equal fit IDs
    if (require.single) {
        nFits <- apply(x[,fit.cols], 1, function(x) { length(unique(as.numeric(x))) })
        x <- x[nFits==1,]
    }
    # clean up protein group representation
    x$ProteinGroup <- gsub(";$", "", x$ProteinGroup)
    r <- list(peptides = x, 
            meta = list(
                    plex = plex,
                    labels = labels,
                    abundance.cols = abundance.cols,
                    quality.cols = quality.cols,
                    fitid.cols = fit.cols
            )
    )
    r
}