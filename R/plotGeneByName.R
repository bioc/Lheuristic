#' plotGeneByName
#'
#' \code{plotGeneByName} plots points on a scatterplot with a 3x3 grid
#' superimposed. The name of a the gene is provided jointly with the
#' MultiAssayExperiment object and used to select the row to be plotted.
#'
#' @param geneName The name of the gene to be plotted.
#' @param  mae A MultiAssayExperiment object containing the methylation
#' and expression data for the specified gene.
#' @param filename If provided, the name of the file to save the results
#' as a PDF; defaults to NULL.
#' @param text4Title A string used as the main title for the plot. Defaults to 
#' `geneName` if set to NULL.
#' @param plotGrid A boolean parameter indicating whether to pass
#' the grid option to the `plotGeneSel` function.
#' @param figs A two-component vector defining the 2-dimensional structure of
#' the plots to be generated.
#'
#' @return a pdf with scatterplots of selected by gene name
#'
#' @keywords plot gene name
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom MultiAssayExperiment assay
#' @import ggplot2
#' @export plotGeneByName
#'
#' @examples
#' # Plot gene by name based on example data
#' # Methylation data
#' methylData <- matrix(runif(50), nrow = 10)
#' colnames(methylData) <- paste0("samp", 1:ncol(methylData))
#' rownames(methylData) <- paste0("gene", 1:nrow(methylData))
#' # Expression data
#' expresData <- matrix(rnorm(50), nrow = 10)
#' colnames(expresData) <- paste0("samp", 1:ncol(methylData))
#' rownames(expresData) <- paste0("gene", 1:nrow(methylData))
#' # ColData
#' colDat <- data.frame(
#'     sampleID = colnames(methylData),
#'     name = letters[1:ncol(methylData)]
#' )
#'
#' rownames(colDat) <- colDat$sampleID
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'     experiments = list(
#'         methylation = methylData,
#'         expression = expresData
#'     ),
#'     colData = colDat
#' )
#' plotGeneByName(gene = "gene7", mae = mae)
#'
plotGeneByName <- function(geneName, mae, filename = NULL,
    text4Title = NULL,
    plotGrid = TRUE, figs = c(2, 2)) {
    # Open PDF if filename is provided
    if (!is.null(filename)) {
    grDevices::pdf(filename)
    }

    # Set plot title   
    plotTitle <- if (!is.null(text4Title)) paste(geneName, text4Title, 
    sep = ", ") else geneName

    # Check if the gene name exists in the expression data 
    #and get the gene position
    if (geneName %in% rownames(assay(mae, "expression"))) {
    genePos <- which(rownames(assay(mae, "expression")) == geneName)
    } else {
    stop("Gene not found in expression data.")
    }

    # Plot using plotGeneSel with ggplot2
    plotGeneSel(
    mae = mae, genePos = genePos, titleText = plotTitle,
    x1 = 1/3, x2 = 2/3, plotGrid = plotGrid
    )

    # Close PDF if filename is provided
    if (!is.null(filename)) {
    grDevices::dev.off()
    }
}
