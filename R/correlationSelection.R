#' correlationSelection: A vector correlation function calculator
#'
#' \code{correlationSelection} Uses the function matCorrs; 
#' given two matrices X (m,n), Y (m,n) this function computes 
#' Pearson and Spearman correlation coefficients and their significance p-values
#' for every pair of row vectors.
#' @param mae a MultiAssayExperiment object containing the methylation
#' @param type specifies the correlation to choose between Spearman
#' and Pearson. Default is Spearman.
#' @param adj logical variable indicating if the p-value returned should be
#' adjusted or not. Default is TRUE, it returns an adjusted p-value.
#' @param pValCutoff the upper limit used for the  p-value. Default is 0.05.
#' @param rCutoff the upper limit used for the correlation coefficient. 
#' Default is 0, no cut off.
#' @param sortByCorrs logical; if TRUE, results are ordered by ascending 
#' p-value. Default set to FALSE.
#'
#' @return dataframe with the correlations selected
#'
#' @keywords Correlation Selection Calculator
#' @export correlationSelection
#' @import stats
#' @importFrom energy dcor
#' @examples
#'
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
#' correlationSelection(mae, pValCutoff = 0.25, rCutoff = 0.1, 
#' type = "Spearman", sortByCorrs = TRUE)

correlationSelection <- function(mae, type = "Spearman", adj = TRUE,
    pValCutoff = 0.05,
    rCutoff = 0, sortByCorrs = FALSE) {
    # Extract the methylation and expression data from the MAE
    X <- assay(mae, "methylation")
    Y <- assay(mae, "expression")

    # Calculate the correlation matrix
    corsMat <- matAllCorrs(X, Y, sortByCorrs = sortByCorrs)

    # Filter based on correlation type
    if (type == "Spearman") {
        selected <- corsMat[, c("r (Sp)", "p (Sp)", 
        "adj.Spear.Pval", "distCor")]
        lShaped <- if (adj) {
        (!is.na(selected[, "r (Sp)"]) & selected[, "r (Sp)"] < rCutoff) &
        (!is.na(selected[, "adj.Spear.Pval"]) & 
        selected[, "adj.Spear.Pval"] < pValCutoff)
    } else {
        (!is.na(selected[, "r (Sp)"]) & 
        selected[, "r (Sp)"] < rCutoff) &
        (!is.na(selected[, "p (Sp)"]) & 
        selected[, "p (Sp)"] < pValCutoff)
    }
    } else {
    selected <- corsMat[, c("r (Pear)", "p (Pear)", 
        "adj.Pear.Pval", "distCor")]
    lShaped <- if (adj) {
        (!is.na(selected[, "r (Pear)"]) & 
        selected[, "r (Pear)"] < rCutoff) &
        (!is.na(selected[, "adj.Pear.Pval"]) & 
        selected[, "adj.Pear.Pval"] < pValCutoff)
    } else {
        (!is.na(selected[, "r (Pear)"]) & 
        selected[, "r (Pear)"] < rCutoff) &
        (!is.na(selected[, "p (Pear)"]) & 
        selected[, "p (Pear)"] < pValCutoff)
    }
    }

    # Combine results and label significant correlations
    selected2 <- data.frame(selected, lShaped)
    colnames(selected2) <- c(colnames(selected), "SigNegCorr")
    return(selected2)
}