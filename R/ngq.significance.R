#'
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