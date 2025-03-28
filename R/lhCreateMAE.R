#' lhCreateMAE: Create a MultiAssayExperiment 
#' with Methylation and Expression Data
#'
#' This function constructs a MultiAssayExperiment 
#' object using methylation 
#' and expression data, automatically ensuring
#' consistency in sample names, 
#' feature names, and exposing dataset names.
#' It internally uses 
#' \code{checkPairing} to validate that the provided
#' datasets have matching 
#' row and column names.
#'
#' @param xDat A numeric matrix containing methylation data 
#' (features in rows, samples in columns).
#' @param yDat A numeric matrix containing expression data 
#' (features in rows, samples in columns).
#' @param xName A string specifying the name of the methylation dataset. 
#' Default: "methylation".
#' @param yName A string specifying the name of the expression dataset. 
#' Default: "expression".
#' @param colData (Optional) A DataFrame containing sample-level metadata. 
#' Default: NULL.
#'
#' @return A MultiAssayExperiment object containing the provided datasets 
#' with sample metadata.
#'
#' @keywords MultiAssayExperiment, Data Integration, Omics Data
#' @export lhCreateMAE
#' @import MultiAssayExperiment
#' @examples
#' 
#' library(MultiAssayExperiment)
#' 
#' # Create synthetic methylation and expression data
#' methylData <- matrix(runif(50), nrow = 10)
#' colnames(methylData) <- paste0("samp", 1:ncol(methylData))
#' rownames(methylData) <- paste0("gene", 1:nrow(methylData))
#'
#' expresData <- matrix(rnorm(50), nrow = 10)
#' colnames(expresData) <- colnames(methylData)
#' rownames(expresData) <- rownames(methylData)
#'
#' # Create sample metadata
#' colDat <- data.frame(sampleID = colnames(methylData),
#' name = letters[1:ncol(methylData)])
#' rownames(colDat) <- colDat$sampleID
#'
#' # Construct MultiAssayExperiment
#' mae <- lhCreateMAE(methylData, expresData, colData = colDat)
#'
#' # Display dataset names
#' names(experiments(mae))
#'
lhCreateMAE <- function(xDat, yDat, xName = "methylation", 
    yName = "expression", colData = NULL) {

    # Validate inputs
    if (!is.matrix(xDat) || !is.matrix(yDat)) {
        stop("xDat and yDat must be numeric matrices")
    }

    # Ensure row and column names match
    if (!identical(rownames(xDat), rownames(yDat)) || 
    !identical(colnames(xDat), colnames(yDat))) {
    stop("Row and column names must match
    between xDat and yDat")
    }

    # Create sample metadata if not provided
    if (is.null(colData)) {
    colData <- DataFrame(sampleID = colnames(xDat))
    rownames(colData) <- colData$sampleID
    }

    # Explicitly create an ExperimentList with names
    experiments <- ExperimentList(setNames(list(xDat, yDat),
    c(xName, yName)))

    # Create MultiAssayExperiment
    mae <- MultiAssayExperiment(
    experiments = experiments,
    colData = colData
    )

    # Print dataset names
    message("Datasets included in the MultiAssayExperiment: ",
        paste(names(experiments(mae)), collapse = ", "))

    return(mae)
}