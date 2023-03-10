#' numScore
#'
#' \code{numScore} A function to score scatterplot using a weight matrix.
#'
#' The scoring does not incorporate logical conditions such as "if xij < C ..."
#' @param aGrid A matrix of counts as computed by `calcFreqs` function
#' @param LShaped A boolean value indicating if this scatterplot can be seen as LShaped.
#' @param aWeightMifL A matrix of weights to score the previous counts if the scatterplot has been classified as L.
#' @param aWeightMifNonL A matrix of weights to score the previous counts if the scatterplot has been classified as non-L
#' @keywords scatterplot weights
#' @export numScore
#' @examples
#' xVecT <- matrix(runif(1000), nrow=100)[1,]
#' yVecT<- matrix(rnorm(1000), nrow=100)[1,]
#' trueFreq <- calcFreqs(xMet=xVecT, yExp=yVecT, x1=1/3, x2=2/3)
#' LShaped <- FALSE
#' weightsIfL    <- matrix (c(2,-1,-99,1,0,-1,1,1,2), nrow=3, byrow=TRUE)
#' weightsIfNonL <- matrix (c(0,-1,-99,0,0,-1,0,0,0), nrow=3, byrow=TRUE)
#' numScore (trueFreq, LShaped, weightsIfL, weightsIfNonL)
#'
numScore <- function(aGrid, LShaped, aWeightMifL, aWeightMifNonL){
  if (!LShaped) {
    scoresM <- aGrid * aWeightMifNonL
  }else{
    scoresM <- aGrid * aWeightMifL
  }
  return(sum(scoresM))
}
