#' toReqMat
#'
#' \code{toReqMat} can be used to turn a matrix of required percentages into
#' a matrix of required counts to facilitate its scoring.
#'
#' @param numPoints Number of points in a scatterplot. Used to turn the required percentages into required counts.
#' @param aReqPercentMat Matrix of required percentages
#' @export toReqMat
#'
#' @examples
#' reqPercentages <- matrix (c(3, 20, 5, 5, 40, 20, 4, 1, 2), nrow=3, byrow=TRUE)
#' numberOfPoints <- 100
#' reqMat <- toReqMat (numPoints=numberOfPoints, aReqPercentMat=reqPercentages)
#'
toReqMat <- function (numPoints, aReqPercentMat){
  return(round(aReqPercentMat*numPoints/100,0)) 
}

