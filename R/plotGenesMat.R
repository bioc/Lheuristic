#' plotGenesMat
#'
#' \code{plotGenesMat} wrapper function for plotting the scatterplots associated
#' with two matrices.
#'
#' @param mae A MultiAssayExperiment object containing the methylation
#' and expression data.
#' @param geneNames A character vector of gene names to plot. 
#' If NULL, all genes are plotted.
#' @param fileName The name of the file used to save the results as a PDF.
#' If NULL, the plot is displayed on the screen.
#' @param text4Title An optional title for the plot, incorporating the gene name
#' and L-shape score. Defaults to NULL.
#' @param x1 The x-coordinate of vertical points on the X-axis. 
#' Expected to contain methylation values ranging between 0 and 1,
#' with default values set to 1/3 and 2/3.
#' @param y1 The y-coordinate of vertical points on the Y-axis.
#' If NULL, these are set to the percentiles of yVec 
#' defined by percY1 and percY2.
#' @param percY1 Values used as defaults for y1 when it is set to NULL.
#' @param x2 The x-coordinate of vertical points on the X-axis.
#' Expected to contain methylation values ranging between 0 and 1,
#' with default values set to 1/3 and 2/3.
#' @param y2 The y-coordinate of vertical points on the Y-axis.
#' If NULL, these are set to the percentiles of yVec 
#' defined by percY1 and percY2.
#' @param percY2 Values used as defaults for y2 when it is set to NULL.
#' @param plotGrid A logical value; defaults to TRUE, indicating 
#' whether to plot gridlines on the graph.
#' @param logicSc A numeric score representing the L-shape score. 
#' Defaults to NULL.
#' @param plotsPerPage A numeric value indicating the number of plots
#' to be printed in each page. Default set to 9.
#' @param saveToPDF Logical, if TRUE, saves the plots to a PDF file 
#' specified by fileName. If FALSE (default), plots are displayed 
#' interactively without saving.
#'
#' @return a pdf with scatterplots for all genes
#'
#' @keywords scatterplot gene plot matrix
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom MultiAssayExperiment assay
#' @import ggplot2
#' @import ggpubr
#' @export plotGenesMat
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
#' sampleID = colnames(methylData),
#' name = letters[1:ncol(methylData)])
#'
#' rownames(colDat) <- colDat$sampleID
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'     experiments = list(
#'         methylation = methylData,
#'         expression = expresData
#'     ),
#'     colData = colDat
#' ) 
#' selectedGenes <- c("gene1", "gene5")

#' plotGenesMat(mae, geneNames = selectedGenes,
#' saveToPDF = FALSE)
#'

plotGenesMat <- function(mae, geneNames = NULL,
    fileName = "Scatplot_Sel_Genes.pdf", 
    text4Title = NULL, x1 = 1/3,x2 = 2/3, y1 = NULL, y2 = NULL,
    percY1 = 1/3, percY2 = 2/3, plotGrid = TRUE,logicSc = NULL,
    saveToPDF = FALSE, plotsPerPage = 9) {

    if (!is.null(geneNames)) {
    mae <- mae[geneNames, , ]}

    if (!is.null(text4Title)) {
    text4Title <- paste(rownames(assay(mae,"expression")),
    text4Title, sep = ", ")} else {
    text4Title <- if (is.null(logicSc)) {
        rownames(assay(mae, "expression"))
    } else {paste(rownames(assay(mae, "expression")),
        "\nL-shaped = ", logicSc, sep = " ")}
    }

    plotList <- list()

    for (gene in seq_len(nrow(assay(mae,"expression")))) {
    plotList[[gene]] <- plotGeneSel(mae = mae,
        genePos = gene, titleText = text4Title[gene],
        x1 = x1,x2 = x2, percY1 = percY1, percY2 = percY2,
        plotGrid = plotGrid)
    }

    if (saveToPDF) {
        grDevices::pdf(fileName)
    
    for (i in seq(1, length(plotList),
        by = plotsPerPage)) {
        ggarrange(plotlist = plotList[i:min(i + plotsPerPage - 1,
        length(plotList))], ncol = 3, 
        nrow = ceiling(plotsPerPage / 3))}

    grDevices::dev.off()} else {
    return(ggarrange(plotlist = plotList, ncol = 3,
        nrow = ceiling(length(plotList) / 3)))}
}