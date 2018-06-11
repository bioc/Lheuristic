vecCorrs <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  return(unlist(list(rhoS=corrS$r[1,2],rhoP=corrP$r[1,2],
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

distCorrs <- function (x, y){
  distCorr <- dcor(x, y, index=1.0)
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
  distCorr <- dcor(x, y, index=1.0)
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


multivCorr <- function(X,Y){
  # stopifnot(require(energy))
  # stopifnot(require(FactoMineR))
  # RV1 <-coeffRV(X, Y)
  # cat("p-value : ", RV1$p.value,"\n")
  # dcor1<-dcor(X, Y)
  # cat("DistCorr: ", dcor1,"\n")
  coin1 <- cia(X, Y)
  cat("RV coeff: ", coin1$coinertia$RV,"\n")
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

DimMat <- function (y, title="", cols2show=5){
  show(dim(y))
  cat("\nFirst rows and columns\n")
  show(head(y[,1:cols2show]))
  cat("\nLast rows first columns\n")
  show(tail(y[,1:cols2show]))
  #  cat("\nLast rows last columns\n")
  #  show(tail(y[,(dim(y)[2]-cols2show+1):(dim(y)[2])]))
}

MI <- function(x,y,h){
  # x : vector of methylations
  # y : vector of expressions
  # h : the std of the Gaussian kernel for density estimation

  if(length(x) != length(y))  stop("Different number of samples!")

  M <- length(x) # samples
  aux <- 0
  two.h.square <- 2*h^2

  for(i in 1:M){
    # kernel distance between the i.th data point and all other points j for each gene k
    tmpx <- x - x[i]
    tmpx <- exp(-tmpx^2/(two.h.square))
    tmpy <- y - y[i]
    tmpy <- exp(-tmpy^2/(two.h.square))
    tmp <- tmpx*tmpy
    aux <- aux + log(M*sum(tmp)/(sum(tmpx)*sum(tmpy)))
  }
  aux/M
}

cMI <- function(dataMeth, dataExp, t, h=0.3){
  # dataMeth : input data methylation (a vector, not a matrix!)
  # dataExp  : input data expression (a vector, not a matrix!)
  # t : moves from 0 to 1
  n <- length(dataMeth)
  if(length(dataExp) != n)  stop("Different number of samples!")
  if(t < 0 | t > 1)  stop("t value is out of range")
  filter <- dataMeth < t
  ss <- sum(filter)
  if(ss != 0){
    x <- dataMeth[filter]
    y <- dataExp[filter]
    aux <- MI(x,y,h)*ss/n
  }
  else{
    aux <- 0
  }
  ss <- sum(!filter)
  if(ss != 0){
    x <- dataMeth[!filter]
    y <- dataExp[!filter]
    aux + MI(x,y,h)*ss/n
  }
  else{
    aux
  }
}

computeCMI <- function (methData, exprData){
  # check consistency
  # Provar un try/catch
  if((nrow(exprData)!=nrow(methData))||(ncol(exprData)!=ncol(methData)))
    stop("Expression amd methylation data must have same dimensions")
  tt <- seq(0,1,by=0.01)
  nt <- length(tt)
  ngenes <- nrow(exprData)
  nsamples <- ncol(exprData)
  cmi <- matrix(rep(0,nt*ngenes), ncol=nt)
  for (j in 1:nrow(methData)){
    # cat (j, " ")
    methVec <- methData[j,]
    exprVec <- exprData[j,]
    # df.xy <- na.omit(data.frame(methVec=methData[j,], exprVec=exprData[j,]))
    for(i in 1:nt){
      cmi[j ,i] <- cMI(methVec, exprVec, t=tt[i])
    }
  }
  return(cmi)
}
