#' checkPairing
#'
#' \code{checkPairing} function to check if two matrices have the same 
#' dimensions row and column names.
#'
#' @param X matrix
#' @param Y matrix
#'
#' @return a logical value indicating if the two matrices have the 
#' same row and column names.
#'
#' @keywords matrix
#'
#' @export checkPairing
#'
#' @examples
#' (X <- round(matrix(rnorm(30) * 10, ncol = 6), 1)) + 1:10
#' (Y <- round(X + matrix(rnorm(30) * 10, ncol = 6), 1)) - 10:1
#' (rownames(X) <- rownames(Y) <- letters[1:nrow(X)])
#' (m1 <- checkPairing(X, Y))
#'
checkPairing <- function(X, Y) {
    # Check if row names and column names match
    row_check <- sum(rownames(X) != rownames(Y)) == 0
    col_check <- sum(colnames(X) != colnames(Y)) == 0

    # Return TRUE if both rows and columns match, FALSE otherwise
    return(row_check && col_check)
}

