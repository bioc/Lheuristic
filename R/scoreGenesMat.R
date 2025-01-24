#' scoreGenesMat
#'
#' \code{scoreGenesMat} scores scatterplots using a binary and a numeric
#' schemes row-wise.
#' @param mae MultiAssayExperiment object containing methylation and expression 
#' matrices.
#' @param aReqPercentsMat A matrix specifying the minimum and maximum percentage
#' counts required in each cell.
#' @param aWeightMifL A matrix of weights applied to score counts when the
#' scatterplot is classified as "L".
#' @param aWeightMifNonL A matrix of weights applied to score counts when the
#' scatterplot is classified as "non-L".
#' @param x1,x2 Coordinates of vertical points on the X-axis. 
#' Expected to contain methylation values ranging between 0 and 1,
#' with default values set to 1/3 and 2/3.
#' @param y1,y2 Coordinates of vertical points on the Y-axis. If set to NULL,
#' they default to the percentiles of yVec as defined by `percY1` and `percY2`.
#' @param percY1,percY2 Values used to act as default for `y1`and `y2` when
#' these are set to `NULL`.
#' @export scoreGenesMat
#' @return A data frame with two columns: 'logicSc' (logical score indicating if
#' the gene is 'active') and 'numericSc' (a numerical score).
#'
#' @examples
#' # Score genes based on example data
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
#' sampleSize <- ncol(MultiAssayExperiment::experiments(mae)[[1]])
#' reqPercentages <- matrix(c(3, 20, 5, 5, 40, 20, 4, 1, 2),
#'     nrow = 3, byrow = TRUE
#' )
#' (theWeightMifL <- matrix(c(2, -2, -sampleSize / 5, 1, 0, -2, 1, 1, 2),
#'     nrow = 3, byrow = TRUE
#' ))
#' (theWeightMifNonL <- matrix(c(0, -2, -sampleSize / 5, 0, 0, -2, 0, 0, 0),
#'     nrow = 3, byrow = TRUE
#' ))
#' scoreGenesMat(mae,
#'     x1 = 1 / 3, x2 = 2 / 3,
#'     y1 = NULL, y2 = NULL, percY1 = 1 / 3, percY2 = 2 / 3,
#'     aReqPercentsMat = reqPercentages,
#'     aWeightMifL = theWeightMifL,
#'     aWeightMifNonL = theWeightMifNonL
#' )
scoreGenesMat <- function(mae, x1 = 1 / 3, x2 = 2 / 3, y1 = NULL, y2 = NULL,
    percY1 = 1 / 3,
    percY2 = 2 / 3, aReqPercentsMat, aWeightMifL = 0.5, aWeightMifNonL = 0.25) {
    if (sum(aReqPercentsMat) != 100) {
        stopifnot(`Percentages must add up to 100` 
            = sum(aReqPercentsMat) == 100)
    }
    stopifnot(prod(names(MultiAssayExperiment::experiments(mae)) == c(
        "methylation",
        "expression"
    )) == 1)
    mets <- MultiAssayExperiment::assay(mae, "methylation")
    expres <- MultiAssayExperiment::assay(mae, "expression")
    N <- ncol(mets)
    Ngenes <- nrow(mets)
    scores <- data.frame(logicSc = rep(FALSE, Ngenes),
        numericSc = rep(0, Ngenes))
    rownames(scores) <- rownames(mets)
    minmaxCounts <- toReqMat(N, aReqPercentMat = aReqPercentsMat)
    indexes <- seq(from = 1, to = Ngenes, by = 1)
    for (gene in indexes) {
        geneNum <- rownames(expres)[gene]
        xMet <- mets[geneNum, ]
        yExp <- expres[geneNum, ]
        geneGrid <- calcFreqs(mae, geneNum, x1, x2,
            y1 = NULL, y2 = NULL, percY1 = 1 / 3,
            percY2 = 2 / 3
        )
        binSc <- binScore(geneGrid, minmaxCounts)
        scores[gene, "logicSc"] <- binSc
        numSc <- numScore(geneGrid, LShaped = binSc, aWeightMifL = aWeightMifL,
            aWeightMifNonL = aWeightMifNonL)
        scores[gene, "numericSc"] <- numSc
    }
    return(scores)
}
