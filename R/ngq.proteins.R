#' Calculate NGQ protein-level fold changes, logratios, deviations and CVs.
#' 
#' This function takes an ngq.peptide data frame, derives protein-level
#' quantification and variance information and returns an ngq.protein data
#' frame.
#' 
#' @export
#' @param x A peptide-quantified ngq.peptide data frame (i.e. processed with
#'          ngq.peptides().
#' @param f A factor column that defines the protein group memberships of the
#'          peptides. This defaults to the ProteinGroup column in the 
#'          ngq.peptide data frame.
#' @return A data frame with all protein-level quantifications
#' 
#' @examples
#'     data(xilac_peptides)
#'     p <- ngq.peptides(xilac_peptides)
#'     q <- ngq.proteins(p)
#'  
ngq.proteins <- function(x, f=x$peptides$ProteinGroup)
{
    ratioCols <- x$meta$fold.cols
    logratioCols <- x$meta$logratio.cols
    peptides <- x$peptides
    proteins.logratios <- NULL
    proteins.mad <- NULL
    proteins.ratios <- NULL
    proteins.cv <- NULL
    nPeptides <- tapply(matrix(1, length(f), 1), f, sum)
    for (i in 1:length(logratioCols)) {
        tmp <- peptides[, logratioCols[i]]
        proteins.logratios <- cbind(proteins.logratios, tapply(tmp, f, median))
        proteins.mad <- cbind(proteins.mad, tapply(tmp, f, mad))
        tmp <- peptides[, ratioCols[i]]
        proteins.cv <- cbind(proteins.cv, tapply(tmp, f, function(x) {mad(x)/median(x)}))
        proteins.ratios <- cbind(proteins.ratios, tapply(tmp, f, median))
    }
    ix <- match(levels(as.factor(f)), f)
    desc <- peptides$Description[ix]
    ipi <- peptides$Protein[ix]
    proteins <- data.frame(
            f[ix],
            ipi,
            desc,
            as.numeric(nPeptides),
            proteins.logratios,
            proteins.mad,
            proteins.ratios,
            proteins.cv
    )
    rownames(proteins) <- NULL
    logRatioNames <- names(x$peptides[logratioCols])
    ratioNames <- names(x$peptides[ratioCols])
    names(proteins) <- c("GeneName", "AccNumber", "Description",
            "nPeptides",
            logRatioNames,
            gsub("^", "MAD ", logRatioNames),
            ratioNames,
            gsub("^", "CV ", ratioNames))
    n <- length(logRatioNames) 
    col.ix <- 4 + matrix(1:(4*n), n, 4)
    meta <- list(
            f = f,
            cols.logratio=col.ix[,1],
            cols.mad=col.ix[,2],
            cols.fold=col.ix[,3],
            cols.cv=col.ix[,4]
    )
    p <- list(
            proteins=proteins,
            meta=meta
    )
    p
}