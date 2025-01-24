#' calcFreqs
#'
#' \code{calcFreqs} Given a MultiAssayExperiment with methylation and
#' expression data, the function overlays a grid on a YMet ~ Xmet scatterplot
#' and returns a 3x3 matrix with point counts per grid cell based on  the
#' vertical and horizontal lines.
#'
#' @param mae MultiAssayExperiment object containing methylation and expression
#' matrices.
#' @param geneNum row of expression/methylation matrix for which the frequencies
#' will be computed.
#' @param x1,x2 Coordinates of vertical points in the X axis. Because it is
#' expected to contain.
#' methylation values that vary between 0 and 1, default values are 1/3 and 2/3.
#' @param y1,y2 Coordinates of vertical points on the Y-axis. If set to NULL,
#' they will be automatically assigned as the percentiles of `y` defined 
#' by `percY1` and `percY2`.
#' @param percY1,percY2 Default values for `y1` and `y2` when set to `NULL`
#'
#' @return a matrix with calculated frequencies
#'
#' @keywords calculation frequencies
#' @export calcFreqs
#'
#' @examples
#'
#' # Methylation data
#' methylData <- matrix(runif(50), nrow = 10)
#' colnames(methylData) <- paste0(
#'     "samp",
#'     1:ncol(methylData)
#' )
#' rownames(methylData) <- paste0(
#'     "gene",
#'     1:nrow(methylData)
#' )
#' # Expression data
#' expresData <- matrix(rnorm(50), nrow = 10)
#' colnames(expresData) <- paste0(
#'     "samp",
#'     1:ncol(methylData)
#' )
#' rownames(expresData) <- paste0(
#'     "gene",
#'     1:nrow(methylData)
#' )
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
#' geneRow <- 1
#' x1 <- 1 / 3
#' x2 <- 2 / 3
#' y1 <- NULL
#' y2 <- NULL
#' percY1 <- 1 / 3
#' percY2 <- 2 / 3
#'
#' calcFreqs(
#'     mae = mae, geneNum = geneRow,
#'     x1, x2, y1, y2, percY1, percY2
#' )
calcFreqs <- function(mae, geneNum, x1, x2, y1 = NULL, y2 = NULL,
    percY1 = 1 / 3, percY2 = 2 / 3) {
    # if (prod(names(mae@ExperimentList) == c('methylation' ,'expression'
    # ))!=1) stop('Error: Names of layers must be 'methylation' and
    # 'expression'')
    mets <- MultiAssayExperiment::assay(mae, "methylation")
    expres <- MultiAssayExperiment::assay(mae, "expression")
    xMet <- mets[geneNum, ]
    yExp <- expres[geneNum, ]
    freqsMat <- matrix(0, nrow = 3, ncol = 3)
    xVals <- c(x1, x2)
    minExp <- min(yExp)
    maxExp <- max(yExp)
    delta <- maxExp - minExp
    if (is.null(y1)) {
        y1 <- minExp + percY1 * delta
    }
    if (is.null(y2)) {
        y2 <- minExp + percY2 * delta
    }
    yVals <- c(y1, y2)
    condX <- c("(xMet<=x1)", "((xMet>x1) &
                (xMet<=x2))", "(xMet>x2)")
    condY <- c("(yExp>y2)", "((yExp<=y2) &
                (yExp>y1))", "(yExp<=y1)")
    for (i in seq_len(3)) {
        for (j in seq_len(3)) {
            condij <- paste(condX[j], condY[i], sep = "&")
            freqsMat[i, j] <- sum(eval(parse(text = condij)))
        }
    }
    return(freqsMat)
}
