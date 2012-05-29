#' Calculate NGQ peptide-level (log-)ratios.
#' 
#' The function takes a raw ngq.peptides data frame and appends the
#' fold change and logratio (base 2) columns.
#' 
#' @export
#' @param x A raw ngq.peptides data frame.
#' @param ref The name of the reference label whose value is used as the denominator
#'            in all fold change and logratio calculation steps. The value
#'            this parameter must be among x$meta$labels.
#' @return The ngq.peptides data frame with fold change and logratio columns
#'         appended.
#' @examples
#'     library(Rngq)
#'     data(xilac_peptides)
#'     p <- ngq.peptides(xilac_peptides)
#' 
ngq.peptides <- function(x, ref=x$meta$labels[1])
{
    # make sure the refernce is among the labels
    stopifnot(ref %in% x$meta$labels)
    # store the reference info
    x$meta$ref <- ref
    # get the non-reference labels
    non.ref <- setdiff(x$meta$labels, ref)
    # get the log2 of the reference column
    log2.ref = log2(x$peptides[,ref])
    # prepare a vector of empty column names
    ratio.names <- rep("", x$meta$plex - 1)
    # prepare an empty matrix of logratios
    logratios <- matrix(NA, length(log2.ref), x$meta$plex - 1)
    # iterate over all non-reference labels and determine the logratio
    for (nr in 1:length(non.ref)) {
        logratios[,nr] <- log2(x$peptides[,non.ref[nr]]) - log2.ref
        ratio.names[nr] <- paste(non.ref[nr], ref, sep="/")
    }
    # calculate indexes
    n.cols <- ncol(x$peptides)
    x$meta$logratio.cols <- (n.cols + 1):(n.cols + x$meta$plex - 1)
    x$meta$fold.cols <- x$meta$logratio.cols + x$meta$plex - 1
    # merge into one big data.frame
    cn <- colnames(x$peptides)
    x$peptides <- cbind(x$peptides, logratios, 2^logratios)
    # adjust the column names
    colnames(x$peptides) <- c(cn, paste("log2", ratio.names), 
            paste("fold", ratio.names)
    )
    x
}