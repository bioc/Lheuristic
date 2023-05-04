#' correlationSelection: A vector correlation function calculator
#'
#' \code{correlationSelection} Uses the function matCorrs; given two matrices X (m,n) , Y (m,n) this function computes Pearson and Spearman correlation coefficients
#' and their significance p-values  for every pair of row vectors.
#'
#' @param X First matrix
#' @param Y Second matrix. Must have the same dimensions as X.
#' @param type specifies the correlation to choose between Spearman and Pearson. Default is Spearman.
#' @param adj logical variable indicating if the p-value returned should be adjusted or not. Default set to TRUE, which will return an adjusted p-value.
#' @param pValCutoff the upper limit to be used for the  p-value. Default is 0.05.
#' @param rCutoff the upper limit to be used for the correlation coefficient. Default is 0, no cut off.
#' @param sortByCorrs logical; if TRUE, results will be ordered in ascending order by p-value. Default set to FALSE.
#'
#' @return dataframe with correlations selected
#'
#' @keywords Correlation Selection Calculator
#' @export correlationSelection
#' @import stats
#' @importFrom energy dcor
#' @examples
#' 
#' (X <- round(matrix (rnorm(30)*10, ncol=6),1))
#' (Y <- round(X + matrix (rnorm(30)*10, ncol=6),1))
#' (rownames(X)=rownames(Y)=letters[1:nrow(X)])
#' correlationSelection(X, Y, pValCutoff=0.25, rCutoff=0.1, type="Spearman", sortByCorrs=TRUE)
#' 
#' 
correlationSelection <- function (X, Y, type="Spearman",
                            adj=TRUE, pValCutoff=0.05, rCutoff=0,
                            sortByCorrs=FALSE){
  corsMat <- matAllCorrs (X, Y, sortByCorrs=sortByCorrs)
  if (type=="Spearman"){
    selected <- corsMat[,c("r (Sp)", "p (Sp)", "adj.Spear.Pval", "distCor")]
    if(adj){
      lShaped <-(selected[,"r (Sp)"]< rCutoff) & (selected[,"adj.Spear.Pval"] < pValCutoff)
    }else{
      lShaped <-(selected[,"r (Sp)"]< rCutoff) & (selected[,"p (Sp)"] < pValCutoff)
    }
  }else{
    selected <- corsMat[,c("r (Pear)", "p (Pear)", "adj.Pear.Pval", "distCor")]
    if(adj){
      lShaped <-(selected[,"r (Pear)"]< rCutoff) & (selected[,"adj.Pear.Pval"] < pValCutoff)
    }else{
      lShaped <-(selected[,"r (Pear)"]< rCutoff) & (selected[,"p (Pe)"] < pValCutoff)
    }
  }
  selected2 <- data.frame (selected, lShaped)
  colnames(selected2) <- c(colnames(selected), "SigNegCorr")
  return(selected2)
}




