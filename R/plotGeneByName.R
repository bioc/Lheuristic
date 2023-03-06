#' plotGeneByName
#'
#' \code{plotGeneByName} plots points on a scatterplot with a 3x3 grid superimposed.
#' The name of a the gene is provided jointly with the matrix and used to select the row to be plotted.
#'
#' @param geneName name of the gene to be plotted.
#' @param mets matrix containing the methylation data of the named gene.
#' @param expresmatrix containing the expression data of the named gene.
#' @param filename NULL, name of the file to store the results as pdf if a name is passed.
#' @param text4Title String to be used as main title for the plot. If NULL it is set to geneName
#' @param plotGrid boolean parameter to be passed to `plotGeneSel` function
#' @param figs 2-components vector defining the 2-dimensional structure of the plots that will be drawn
#' @keywords plot gene name
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export plotGeneByName
#'
#' @examples
#' \dontest{
#' plotGeneByName (gene="HOOK1", mets=falseLMet, expresmatrix=falseLExpr, fileName=NULL)
#'}
#'
plotGeneByName <- function(geneName, mets, expresmatrix, filename, text4Title=NULL,
                           plotGrid=TRUE, figs=c(2,2)){
  if (!is.null(filename))
    grDevices::pdf(filename)
  if (!is.null(text4Title)){
    text4Title <- paste(geneName, text4Title, sep=", ")
  }else{
    text4Title<- geneName
  }
  if (geneName %in% rownames(expresmatrix)){
    genePos <- which(rownames(expresmatrix)==geneName)
  }else{
    genePos <- NULL
  }
  if(!(is.null(genePos))){
    xVec<- as.numeric(mets[genePos,])
    yVec<- as.numeric(expresmatrix[genePos,])
    plotGeneSel(xMet=xVec, yExp=yVec, titleText=text4Title,
               x1=1/3, x2=2/3, plotGrid=plotGrid)
  }
  if (!is.null(filename))
    grDevices::dev.off()
}



