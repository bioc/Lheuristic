#' plotGeneSel
#'
#' \code{plotGeneSel} plots points on a scatterplot with a 3x3 grid overimposed.
#'
#' @param mae A MultiAssayExperiment object containing the methylation
#' and expression data for the specified gene.
#' @param genePos The index of the gene to be plotted within the
#' MultiAssayExperiment object.
#' @param titleText plot title.
#' @param x1,x2 Coordinates of vertical points in the X axis. 
#' Because it is expected to contain methylation values that vary
#' between 0 and 1. The default values are 1/3 and 2/3.
#' @param y1,y2 Coordinates of vertical points in the Y axis. 
#' Leaving them as NULL assigns them the percentiles of yVec defined by `percY1`
#' and `percY2`.
#' @param percY1,percY2 Values used to act as default for `y1`and `y2`
#' when these are set to `NULL`
#' @param plotGrid logical. Defautl to TRUE will plot gridlines over the 
#' scatterplot.
#'
#' @return a pdf with scatterplots for selected genes
#'
#' @keywords plot gene selection
#' @importFrom graphics abline
#' @importFrom MultiAssayExperiment assay
#' @import ggplot2
#' @export plotGeneSel
#'
#' @examples
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
#' 
#' plotGeneSel(mae, genePos = 7, titleText = "L-shaped gene")
#'
plotGeneSel <- function(mae, genePos, titleText, x1 = 1/3, 
    x2 = 2/3, y1 = NULL, y2 = NULL, percY1 = 1/3, 
    percY2 = 2/3, plotGrid = TRUE) {
    xVec <- as.numeric(MultiAssayExperiment::assay(mae, 
        "methylation")[genePos, ])
    yVec <- as.numeric(MultiAssayExperiment::assay(mae, 
        "expression")[genePos, ])
    minExp <- min(yVec)
    maxExp <- max(yVec)
    delta <- maxExp - minExp
    # Create a data frame for ggplot
    Methylation <- xVec
    Expression <- yVec
    plotData <- data.frame(Methylation, Expression)

    # Set y1 and y2 based on percentiles if not provided
    if (is.null(y1)) {
        y1 <- minExp + percY1 * delta
    }
    if (is.null(y2)) {
        y2 <- minExp + percY2 * delta
    }

    # Base ggplot scatterplot
    p <- ggplot2::ggplot(plotData, ggplot2::aes(x = Methylation, 
    y = Expression)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = titleText) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(minExp, maxExp)

    # Add grid lines if plotGrid is TRUE
    if (plotGrid) {
        p <- p +
        ggplot2::geom_vline(xintercept = c(x1, x2), 
        linetype = "dashed", color = "gray") +
        ggplot2::geom_hline(yintercept = c(y1, y2), 
        linetype = "dashed", color = "gray") + 
        ggplot2::theme_minimal()
    }

    # Show the plot
    p
}