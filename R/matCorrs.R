#' A row-wise correlation function calculator
#'
#' \code{matCorrs} Given two matrices X (m,n), Y(m,n) this function computes (m) Pearson and Spearman correlation coefficients
#' and their significance p-values for every pair of row vectors.
#'
#' @param X First matrix
#' @param Y Second matrix. Must have the same dimensions as X
#' @keywords correlation
#' @import Hmisc
#' @export matCorrs
#' @examples
#' (X <- round(matrix (rnorm(30)*10, ncol=6),1)) + 1:10
#' (Y <- round(X + matrix (rnorm(30)*10, ncol=6),1)) - 10:1
#' (rownames(X)=rownames(Y)=letters[1:nrow(X)])
#' (m1<-matCorrs(X,Y))
#' matCorrs
#'
#'

matCorrs <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  corrsList<- matrix(NA, nrow=nrow(X), ncol=4)
  colnames(corrsList) <- c("r (Sp)", "r (Pear)",  "p (Sp)","p (Pear)")
  for (i in seq_along(1:nrow(X))){
    corrs<- vecCorrs(X[i,],Y[i,])
    corrsList[i,] <- corrs
  }
  rownames(corrsList) <- rownames(X)
  return(corrsList)
}

