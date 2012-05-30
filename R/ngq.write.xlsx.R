#' Store the NGQ protein and peptide quantification in an
#' MS Excel zipped XML file (.xlsx).
#'
#' Stores the peptide- and/or protein-level information into
#' an MS Excel file. If \code{peptides} and \code{proteins}
#' are both specified, this will store the data in two different
#' sheets, but in the same file.
#' 
#' @export
#' @param filename The name of the Excel file.
#' @param peptides A \code{peptides} data frame with processed (i.e.
#'                 quantified) peptide data. See \code{ngq.peptides}.
#' @param proteins A \code{proteins} data frame with protein-level
#'                 quantification data.
#' 
#' @section Details: Either \code{peptides} or \code{proteins} must
#' be non-\code{NULL} or the function will fail.
#' 
#' @section Excel format: The function depends on the \code{XLConnect}
#' package in order to store data in the Excel format.
#' \code{XLConnect} is based on a Java core; should writing the file
#' run out of memory, setting suitable parameters for the JVM may help
#' (see examples).
#' 
#' @examples \dontrun{
#'     data(xilac_peptides)
#'     p <- ngq.peptides(xilac_peptides)
#'     # set up the virtual machine; this must be done
#'     # prior to loading the XLConnect library, i.e. prior
#'     # to the first call to ngq.write.xlsx
#'     options(java.parameters = "-Xmx1024m")
#'     ngq.write.xlsx("somefile.xlsx", peptides=p)
#' }
#' 
ngq.write.xlsx <- function(filename, peptides=NULL, proteins=NULL)
{
    stopifnot(!(is.null(peptides) && is.null(proteins)))
    # store into excel file; 1 peptide, 1 protein sheet
    library(XLConnect)
    unlink(filename)
    wb <- loadWorkbook(filename, create=T)
    if (!is.null(peptides)) {    
        # create a new sheet
        sheetName <- "peptides"
        wb$createSheet(name=sheetName)
        # define the destination region
        regionName <- paste(sheetName, "container", sep="_") 
        createName(wb, name=regionName,
                formula = paste(sheetName, "$A$1", sep = "!"))
        writeNamedRegion(wb, data = peptides$peptides, name=regionName, header = TRUE)
    }
    if (!is.null(proteins)) {
        # create a new sheet for the protein data
        sheetName <- "proteins"
        wb$createSheet(name=sheetName)
        # define the destination region
        regionName <- paste(sheetName, "container", sep="_") 
        createName(wb, name=regionName,
                formula = paste(sheetName, "$A$1", sep = "!"))
        writeNamedRegion(wb, data = proteins, name=regionName, header = TRUE)
    }
    # TODO: store stats
    #sheetName <- "stats"
    #wb$createSheet(name=sheetName)
    #regionName <- paste(sheetName, "basicStats", sep="_") 
    #createName(wb, name=regionName,
    #        formula = paste(sheetName, "$A$1", sep = "!"))
    #writeNamedRegion(wb, data = stats$featureExtraction, name=regionName, header = TRUE)
    #regionName <- paste(sheetName, "proteinStats", sep="_") 
    #createName(wb, name=regionName,
    #        formula = paste(sheetName, "$A$4", sep = "!"))
    #writeNamedRegion(wb, data = stats$proteinQuantification, name=regionName, header = TRUE)
    #regionName <- paste(sheetName, "peptideStats", sep="_") 
    #createName(wb, name=regionName,
    #        formula = paste(sheetName, "$A$7", sep = "!"))
    #tmp <- data.frame(rownames(stats$peptideQuantification), stats$peptideQuantification)
    #names(tmp) <- c("Label state", names(stats$peptideQuantification))
    #writeNamedRegion(wb, data = tmp, name=regionName, header = TRUE)
    saveWorkbook(wb)
}