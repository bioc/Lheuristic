vecCorrs <- function(x, y) {
    corrS <- rcorr(x, y, type = "spearman")
    corrP <- rcorr(x, y, type = "pearson")
    return(unlist(list(
        rhoS = corrS$r[1, 2],
        rhoP = corrP$r[1, 2], pvalS = corrS$P[1, 2],
        pvalP = corrP$P[1, 2]
    )))
}

distCorrs <- function(x, y) {
    distCorr <- energy::dcor(x, y, index = 1)
    return(unlist(list(distCorr)))
}

matDistCorr <- function(X, Y) {
    if ((nrow(X) != nrow(Y)) || (ncol(X) != ncol(Y))) {
        stop("matrices dimensions do not match")
    }
    distcorL <- matrix(NA, nrow = nrow(X), ncol = 1)
    colnames(distcorL) <- c("DistCor")
    for (i in seq_len(nrow(X))) {
        DistCorr <- distCorrs(X[i, ], Y[i, ])
        distcorL[i, ] <- DistCorr
    }
    rownames(distcorL) <- rownames(X)
    return(distcorL)
}

allCorrs <- function(x, y) {
    corrS <- rcorr(x, y, type = "spearman")
    corrP <- rcorr(x, y, type = "pearson")
    distCorr <- energy::dcor(x, y, index = 1)
    return(unlist(list(
        rhoS = corrS$r[1, 2],
        rhoP = corrP$r[1, 2],
        distCorr, pvalS = corrS$P[1, 2], pvalP = corrP$P[1, 2]
    )))
}

sort1 <- function(X, col, DEC = TRUE, ...) {
    return(X[sort.list(X[, col], 
    decreasing = DEC), ])
}

matAllCorrs <- function(X, Y, sortByCorrs = FALSE) {
    if ((nrow(X) != nrow(Y)) || (ncol(X) != ncol(Y))) {
        stop("matrices dimensions do not match")
    }
    corL <- matrix(NA, nrow = nrow(X), ncol = 5)
    colnames(corL) <- c("r (Sp)", "r (Pear)", 
        "distCor", "p (Sp)", "p (Pear)")
    for (i in seq_len(nrow(X))) {
        corrs <- allCorrs(X[i, ], Y[i, ])
        corL[i, ] <- corrs
    }
    corL <- cbind(corL, 
        adj.Spear.Pval = p.adjust(corL[, "p (Sp)"], "fdr"))
    corL <- cbind(corL, 
        adj.Pear.Pval = p.adjust(corL[, "p (Pear)"], "fdr"))
    rownames(corL) <- rownames(X)
    if (sortByCorrs) {
        corL <- sort1(corL, 4, DEC = FALSE)
    }
    return(corL)
}


zeros <- function(x) {
    which(x == 0)
}
discard <- function(x, where, howMany) {
    return(length(zeros(x[where])) > howMany)
}

discardA <- function(x, percentage) {
    maxZeros <- ceiling(length(x) * percentage)
    return(discard(x, seq_len(length(x)), maxZeros))
}


countNAs <- function(X) {
    sumaX <- apply(X, 1, sum)
    # sum(is.na(sumaX))
    xAmbNAs <- X[is.na(sumaX), ]
    NAs <- apply(xAmbNAs, 1, function(unVec) {
        sum(is.na(unVec))
    })
}


read2 <- function(expresFName, metFName, dataDirectory = ".", 
    sepChar = ";", decChar = ".") {
    expres <- utils::read.table(
        file = file.path(dataDirectory, expresFName), header = TRUE,
        sep = sepChar, dec = decChar, row.names = 1
    )
    mets <- utils::read.table(
        file = file.path(dataDirectory, metFName), header = TRUE,
        sep = sepChar, dec = decChar, row.names = 1
    )
    if (checkPairing(expres, mets)) {
        result <- list(expres, mets)
    } else {
        result <- NULL
    }
    return(result)
}
