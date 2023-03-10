#' binScore
#'
#' \code{binScore} can be used to score scatterplots by directly comparing
#' the sample counts with a matrix of minimal or maximal percentages/counts
#' to be found in each cell. It implements the three bands rule implicitly
#' by setting threshold values.
#'
#' @param aGrid A matrix of counts as computed by `calcFreqs` function
#' @param aReq A matrix of minimum or maximum counts to be found in each cell
#' if L-shape is TRUE
#'
#' @keywords binary scoring
#' @export binScore
#' @examples
#' # Generate some example data
#' aGrid <- matrix(c(20, 3, 0, 10, 2, 2, 20, 10, 20), nrow = 3, ncol = 3, byrow = TRUE)
#' aReq <- matrix(c(15, 5, 0, 0, 5, 5, 10, 10, 15), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Calculate the bin score
#' binScore(aGrid, aReq)
#'


binScore <- function(aGrid, aReq){
  comp <- (aGrid[1,1]>= aReq[1,1])&&(aGrid[1,2]<= aReq[1,2])&&(aGrid[1,3]<= aReq[1,3]) &&
    (aGrid[2,1]>=aReq[2,1]) &&(aGrid[2,2]<=aReq[2,2]) &&(aGrid[2,3]<=aReq[2,3]) &&
    (aGrid[3,1]>=aReq[3,1])&&(aGrid[3,2]>=aReq[3,2])&&(aGrid[3,3]>=aReq[3,3])
  return(comp)
}

