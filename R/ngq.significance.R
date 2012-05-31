#' Significance calculation for protein up- or downregulation.
#' 
#' The function compares the peptide-level abundance distribution of every
#' protein against the overall protein-level abundance distribution and
#' determines a p-value for the difference between the two distributions.
#' 
#' @export
#' @param peptides A \code{peptides} data frame.
#' @param proteins A \code{proteins} data frame.
#' @param use.test The name of the test that should be used to calculate
#'                 p-values (e.g. \code{wilcox.test}, \code{t.test}).
#'                 Defaults to the non-parametric Wilcoxon test.
#' @param p.adjust.method The method that should be used to correct the 
#'                        p-values for multiple testing. See \code{p.adjust}
#'                        for a list of valid parameters. Note that this 
#'                        defaults to
#'                        \code{none}, i.e. no multiple testing correction!
#' 
#' @section Details: The underlying assumption behind this test setup is that
#' the protein-level abundance logratio distributions are somewhat robust
#' estimates of the null distribution. If the abundance distribution of the
#' peptides of a specific protein deviates from this distribution, this is an
#' indicator that it is indeed differentially regulated. In essence, this means
#' that we rely on the smoothing power of the mean/median.
#'
#' @examples 
#'     data(xilac_peptides)
#'     p <- ngq.peptides(xilac_peptides)
#'     q <- ngq.proteins(p)
#'     s <- ngq.significance(p, q, p.adjust.method="BH") 
#'
ngq.significance <- function(peptides, proteins, use.test=wilcox.test,
        p.adjust.method="none")
{
    # hugely convenient helper function
    pvalue.or.na <- function(x, y=NULL, test=use.test) {
        t <- try(test(na.omit(x), y), TRUE)
        if (is.null(t) || inherits(t, "try-error")) {
            p <- NA
        } else {
            p <- t$p.value
        }
        p
    }
    
    refCols <- proteins$meta$cols.logratio
    pepCols <- peptides$meta$logratio.cols
    n <- length(refCols)
    pvalues <- matrix(NA, nrow(proteins$proteins), n)
    f <- as.factor(proteins$meta$f)
    f.levels <- levels(f)
    for (i in 1:n) {
        ref <- proteins$proteins[[refCols[i]]]
        for (j in 1:length(ref)) {
            peps <- peptides$peptides[[ pepCols[i] ]]
            pvalues[j,i] <- pvalue.or.na(ref, peps[f==f.levels[j]], test=use.test)
        }
    }
    # adjust p-values in every column
    adjusted.pvalues <- apply(pvalues, 2, p.adjust, method=p.adjust.method)
    
    colnames <- c(names(proteins$proteins), paste("q-values", names(proteins$proteins)[proteins$meta$cols.logratio]))
    proteins$proteins <- cbind(proteins$proteins, adjusted.pvalues)
    names(proteins$proteins) <- colnames
    proteins$meta$pvalue.adjust.method <- "BH"
    proteins$meta$cols.qvalue <- proteins$meta$cols.cv + length(proteins$meta$cols.cv)
    proteins
}