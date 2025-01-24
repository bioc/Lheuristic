#' numScore
#'
#' \code{numScore} A function to score scatterplot using a weight matrix.
#' The scoring does not incorporate logical conditions such as 'if xij < C ...'
#' @param aGrid A matrix of counts computed by the `calcFreqs` function.
#' @param LShaped A boolean indicating whether the scatterplot can be
#' classified as "L-shaped".
#' @param aWeightMifL A matrix of weights applied to score counts if
#' the scatterplot is classified as "L".
#' @param aWeightMifNonL A matrix of weights applied to score counts if
#' the scatterplot is classified as "non-L".
#'
#' @return a numeric score for each scatterplot
#'
#' @keywords scatterplot weights
#'
#' @export numScore
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
#' trueFreq <- calcFreqs(mae, x1 = 1 / 3, x2 = 2 / 3)
#' LShaped <- FALSE
#' weightsIfL <- matrix(c(2, -1, -99, 1, 0, -1, 1, 1, 2),
#'     nrow = 3, byrow = TRUE
#' )
#' weightsIfNonL <- matrix(c(0, -1, -99, 0, 0, -1, 0, 0, 0),
#'     nrow = 3, byrow = TRUE
#' )
#' numScore(
#'     trueFreq, LShaped,
#'     weightsIfL, weightsIfNonL
#' )
#'
numScore <- function(aGrid, LShaped, aWeightMifL, aWeightMifNonL) {
    if (!LShaped) {
        scoresM <- aGrid * aWeightMifNonL
    } else {
        scoresM <- aGrid * aWeightMifL
    }
    return(sum(scoresM))
}
