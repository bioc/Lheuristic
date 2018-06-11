#' cmiSelection: Wrapper function to compute CMI and select L-shaped genes
#'
#' \code{cmiSelection}cMI function computes cMI values for different t values, from 0 to 1 and a step of 0.01. Output is stored in a data frame called cMI_values. For each gene, the optimal threshold is the t-value that results the minimum CMI.
#'
#' @param methData matrix containing methylation data
#' @param exprData matrix containing expression data
#' @param h number used for tuning  of the kernel width, and empirically set at a default of \code{$h = 0.2$}.
#' @param smallR numeric value representing the ratio mincM \code{$I(t)/cM I(0)$}. Default is 0.25.
#' @param minCMI Minimum value of unconditioned CMI. default set at 0.1.
#' @keywords conditional mutual information CMI.
#' @export cmiSelection
#' @details   L-shaped genes (regulated by methylation) are selected according to three conditions
#' (parameters were chosen according to a random permutation test as in Liu 2012).
#' 1. Ratio mincM \code{$I(t)/cM I(0)$} must be small enough, \code{$r<0.25$}.
#' 2. Minimum  value of unconditioned CMI must be large enough, \code{$cMI(0)>0.1$}.
#' 3. Expression values must be higher on the left side of the plot than on the right side.
#' @examples
#' \dontrun {
#' methData<- metthylarray; exprData<- exprarray; h=0.2; smallR=0.25; minCMI=0.1
#' cmiSelection
#' }
cmiSelection <- function (methData, exprData, h=0.2, smallR=0.25, minCMI=0.1){

  t <- seq(0,1,by=0.01)
  cMI_values <- matrix(NA,nrow = nrow(exprData),
                       ncol = length(t)+4,
                       byrow = TRUE)
  rownames(cMI_values) <- rownames(exprData)
  colnames(cMI_values) <- c(t,"cMI_min", "t_opt", "ratio", "meth_regulated")
  cMI_values <- as.data.frame(cMI_values)
  for (i in 1:nrow(methData)){
    dataMeth <- methData[i,]
    dataExp <- exprData[i,]
    for (j in 1:length(t)){
      cMI_values[i,j] <- cMI(dataMeth,dataExp,t[j], h=h)
    }
  }
  # columns cMI_min and t_opt store minimum cMI value computed above and t value used for it, respectively.
  for (i in 1:nrow(cMI_values)){
    cMI_values$cMI_min[i] <- min(cMI_values[i,1:length(t)])
  }
  for (i in 1:nrow(cMI_values)){
    cMI_values$t_opt[i] <- names(which.min(cMI_values[i,1:length(t)]))
  }

  for (i in 1:nrow(cMI_values)){
    cMI_values$ratio [i] <- cMI_values$cMI_min[i]/cMI_values[i,1]
  }
  mean_exp_left <- c()
  mean_exp_right <- c()
  for (i in 1:nrow(cMI_values)){
    dataMeth <- methData[i,]
    dataExp <- exprData[i,]
    filt <- dataMeth <= cMI_values$t_opt[i]
    right <- dataExp[!filt]
    left <- dataExp[filt]
    if (sum(left) == 0){
      mean_exp_left[i] <- 0
      mean_exp_right[i] <- mean(right)
    }
    if (sum(right) == 0){
      mean_exp_right[i] <- 0
      mean_exp_left[i] <- mean(left)
    }
    else{
      mean_exp_left[i] <- mean(left)
      mean_exp_right[i] <- mean(right)
    }
  }
  # Genes are labeled according three conditions above:
  for (i in 1:nrow(methData)){
    if (cMI_values$ratio[i] < smallR &
        cMI_values[i,1] > minCMI &
        mean_exp_left[i] > mean_exp_right[i]){
      cMI_values$meth_regulated[i] <- TRUE
    }
    else{
      cMI_values$meth_regulated[i] <- FALSE
    }
  }
  return(cMI_values[,c("cMI_min", "t_opt", "ratio", "meth_regulated")])
}
