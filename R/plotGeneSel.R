#' plotGeneSel
#'
#' \code{plotGeneSel} plots points on a scatterplot with a 3x3 grid overimposed.
#'
#' @param xMet vector with methylation data.
#' @param yExp vector for expression data.
#' @param titleText plot title.
#' @param x1,x2 Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1 the default values are 1/3 and 2/3.
#' @param y1,y2 Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY1,percY2 Values used to act as default for `y1`and `y2` when these are set to `NULL`
#' @param plotGrid logical. Defautl to TRUE will plot gridlines over the scatterplot.
#' 
#' @keywords plot gene selection
#' @importFrom graphics abline
#' @export plotGeneSel
#'
#' @examples
#' \dontrun{
#' myGene1 <-rownames(myGenes)[1]
#' xVec<- as.numeric(myMet[myGene1,])
#' yVec<-as.numeric(myExpr[myGene1,])
#' titleT <- paste (myGene1, "(May be GRM)")
#' plotGeneSel(xMet=xVec, yExp=yVec, titleText=titleT, x1=1/3, x2=2/3)
#'}
#'
plotGeneSel <- function(xMet, yExp, titleText,
                       x1=1/3, x2=2/3, y1=NULL, y2=NULL,
                       percY1=1/3, percY2=2/3, plotGrid=TRUE)
{
  minExp<-min(yExp); maxExp <- max(yExp); delta<- maxExp-minExp
  plot(xMet,yExp,  xlim=c(0,1), ylim=c(minExp, maxExp), main=titleText)
  if (plotGrid){
    if (is.null(y1)) y1<- minExp + percY1*delta
    if (is.null(y2)) y2<- minExp + percY2*delta
    graphics::abline(v=x1);  graphics::abline(v=x2)
    graphics::abline(h=y1);  graphics::abline(h=y2)
  }
}

