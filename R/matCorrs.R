#' A row-wise correlation function calculator
#'
#' \code{matCorrs} Given a MultiAssayExperiment object, the function
#' computes Pearson and Spearman correlation coefficients
#' and their significance p-values for every pair of row vectors.
#' @param mae A MultiAssayExperiment object containing the methylation
#' and expression data.
#'
#' @return matrix with a correlation value for each gene
#'
#' @keywords correlation
#'
#' @import Hmisc
#' @export matCorrs
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
#' matCorrs(mae)
#'
matCorrs <- function(mae) {
    corrsList <- matrix(NA, nrow = nrow(assay(mae, "methylation")), ncol = 4)
    colnames(corrsList) <- c("r (Sp)", "r (Pear)", "p (Sp)", "p (Pear)")
    rownames(corrsList) <- rownames(assay(mae, "methylation"))

    for (i in seq_len(nrow(assay(mae, "methylation")))) {
        xVec <- as.numeric(assay(mae, "methylation")[i, ])
        yVec <- as.numeric(assay(mae, "expression")[i, ])
        corrs <- vecCorrs(xVec, yVec)
        corrsList[i, ] <- corrs
    }

    return(corrsList)
}



