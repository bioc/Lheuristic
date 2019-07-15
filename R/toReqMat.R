#' toReqMat
#' \code{toReqMat} can be used to turn a matrix of required percentages into
#' a matrix of required counts to facilitate its scoring
#' @param numPoints Number of points in a scatterplot. Used to turn the required percentages into required counts.
#' @param aReqPercentMat Matrix of required percentages
#' @examples
#' \dontrun{
#'
#' }
#'
toReqMat <- function (numPoints, aReqPercentMat){
  return(round(aReqPercentMat*numPoints/100,0)) #  return(round(aReqPercentMat*sum(numPoints)/100,0))
}

