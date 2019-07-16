#' checkPairing
#'
#' \code{checkPairing} is a  function to check if two matrices have the same dimensions (same rows and columns) and same row and column names.
#'
#' @param X matrix
#' @param Y matrix
#' @keywords matrix
#'
#' @export checkPairing
#'
#' @examples
#' \dontrun{
#' (X <- round(matrix (rnorm(30)*10, ncol=6),1)) + 1:10
#' (Y <- round(X + matrix (rnorm(30)*10, ncol=6),1)) - 10:1
#' (rownames(X)=rownames(Y)=letters[1:nrow(X)])
#' (m1<-checkPairing(X,Y))
#'}
#'
checkPairing <- function (X, Y){
  check <- (sum(rownames(X)!=rownames(Y)) == 0) && (sum(colnames(X)!=colnames(Y)) == 0)
}

