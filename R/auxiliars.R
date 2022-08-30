vecCorrs <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  return(unlist(list(rhoS=corrS$r[1,2],rhoP=corrP$r[1,2],
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

distCorrs <- function (x, y){
  distCorr <- energy::dcor(x, y, index=1.0)
  return(unlist(list(distCorr)))
}

matDistCorr <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  distCorrsList<- matrix(NA, nrow=nrow(X), ncol=1)
  colnames(distCorrsList) <- c("DistCor")
  for (i in 1:nrow(X)){
    DistCorr<- distCorrs(X[i,],Y[i,])
    distCorrsList[i,] <- DistCorr
  }
  rownames(distCorrsList) <- rownames(X)
  return(distCorrsList)
}

allCorrs  <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  distCorr <- energy::dcor(x, y, index=1.0)
  return(unlist(list(rhoS=corrS$r[1,2],  rhoP=corrP$r[1,2], distCorr,
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

sort1 <- function (X, col,DEC=TRUE, ...){
  return(X[sort.list(X[,col], decreasing=DEC), ])
}

matAllCorrs  <- function (X, Y, sortByCorrs = FALSE){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  corrsList<- matrix(NA, nrow=nrow(X), ncol=5)
  colnames(corrsList) <- c("r (Sp)", "r (Pear)", "distCor", "p (Sp)",  "p (Pear)")
  for (i in 1:nrow(X)){
    corrs<- allCorrs(X[i,],Y[i,])
    corrsList[i,] <- corrs
  }
  corrsList<- cbind(corrsList, adj.Spear.Pval= p.adjust(corrsList[,"p (Sp)"],"fdr"))
  corrsList<- cbind(corrsList, adj.Pear.Pval = p.adjust(corrsList[,"p (Pear)"],"fdr"))
  rownames(corrsList) <- rownames(X)
  if (sortByCorrs) corrsList <- sort1(corrsList,4, DEC=FALSE)
  return(corrsList)
}


zeros <-function(x){which (x==0)}
discard <- function(x, where, howMany){
  return(length(zeros(x[where])) > howMany)
}

discardA <- function(x, percentage){
  maxZeros <- ceiling (length(x)*percentage)
  return (discard(x,1:length(x),maxZeros))
}


countNAs <- function (X){
  sumaX <- apply(X, 1,sum);
  # sum(is.na(sumaX))
  xAmbNAs <- X[is.na(sumaX),]
  NAs <- apply(xAmbNAs, 1, function(unVec){sum(is.na(unVec))})
}


