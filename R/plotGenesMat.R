#' plotGenesMat
#'
#' \code{plotGenesMat} wrapper function for plotting the scatterplots associated with two matrices.
#'
#' @param mets matrix with methylation data.
#' @param expres matrix with expression data.
#' @param fileName name of the file used to save the results as pdf. If NULL, plot goes to screen
#' @param text4Title NULL, name for the plot title containing the gene name and the L-shape score.
#' @param x1, x-Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1 the default values are 1/3 and 2/3.
#' @param y1, y-Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY1, Values used to act as default for `y1` when it is set to `NULL`.
#' @param x2, x-Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1 the default values are 1/3 and 2/3.
#' @param y2, y-Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY2 Values used to act as default for `y2` when these are set to `NULL`.

#' @param plotGrid logical; default set to TRUE will plot gridlines on the graph.
#' @param logicSc NULL, numeric score representing the L-shape score.
#'
#'@return a pdf with scatterplots for all genes
#'
#' @keywords scatterplot gene plot matrix
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export plotGenesMat
#' @examples
#' mets <- matrix(runif(1000), nrow=100)
#' expres <- matrix(rnorm(1000), nrow=100)
#' rownames(mets) <- paste0("Gene", 1:nrow(mets))
#' rownames(expres) <- paste0("Gene", 1:nrow(expres))
#' # plotGenesMat (mets=mets, expres=expres, fileName = "PlotAllGenes.pdf")
#' 

plotGenesMat <- function(mets, expres, fileName = NULL, text4Title=NULL,
                         x1=1/3, x2=2/3,
                         y1=NULL, y2=NULL,
                         percY1=1/3, percY2=2/3,
                         plotGrid=TRUE, logicSc = NULL){
  if (!is.null(fileName))
    grDevices::pdf(fileName)
  if (!is.null(text4Title)){
    text4Title <- paste(rownames(expres),text4Title, sep=", ")
  }else{
    if (is.null(logicSc)){
      text4Title<- rownames(expres)
    }else{
      text4Title<- paste(rownames(expres), "\n L-shaped = ", logicSc, sep = " ") #text4Title<- rownames(expres)
    }
  }
  #  opt<-par(mfrow=c(2,2))
  for (gene in seq_len(nrow(expres))){
    xVec<- as.numeric(mets[gene,])
    yVec<- as.numeric(expres[gene,])
    plotGeneSel(xMet=xVec, yExp=yVec, titleText=text4Title[gene],
               x1=x1, x2=x2, percY1=percY1, percY2=percY2, plotGrid=plotGrid) #x1=1/3, x2=2/3, plotGrid=plotGrid)
  }
  #  par(opt)
  if (!is.null(fileName))
    grDevices::dev.off()
}

