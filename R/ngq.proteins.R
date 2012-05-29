#' Calculate NGQ protein-level fold changes, logratios, deviations and CVs.
#' 
#' This function takes an ngq.peptide data frame, derives protein-level
#' quantification and variance information and returns an ngq.protein data
#' frame.
#' 
#' @param x A peptide-quantified ngq.peptide data frame (i.e. processed with
#'          ngq.peptides().
#' @param f A factor column that defines the protein group memberships of the
#'          peptides. This defaults to the ProteinGroup column in the 
#'          ngq.peptide data frame.
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
    for (i in 1:length(logratioCols)) {
        tmp <- peptides[, logratioCols[i]]
        proteins.logratios <- cbind(proteins.logratios, tapply(tmp, f, median))
        proteins.mad <- cbind(proteins.mad, tapply(tmp, f, mad))
        tmp <- peptides[, ratioCols[i]]
        proteins.cv <- cbind(proteins.cv, tapply(tmp, f, function(x) {mad(x)/median(x)}))
        proteins.ratios <- cbind(proteins.ratios, tapply(tmp, f, median))
    }
    ix <- match(levels(as.factor(f)), f)
    #desc <- peptides$Description[ix]
    ipi <- peptides$Protein[ix]
    proteins <- data.frame(
            f,
            ipi,
    #        desc,
#            as.numeric(nPeptides),
            proteins.logratios,
            proteins.mad,
            proteins.ratios,
            proteins.cv
    )
#    rownames(proteins) <- NULL
#    names(proteins) <- c("GeneName", "AccNumber",# "Description",
#            "nPeptides",
#            gsub("ratio", "logratio", ratioNames),
#            gsub("ratio", "MAD", ratioNames),
#            ratioNames,
#            gsub("ratio", "CV", ratioNames))
    class(proteins) <- "ngq.proteins"
    proteins
}