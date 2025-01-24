## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## -----------------------------------------------------------------------------
library(Lheuristic)

# library(kableExtra)
# library(VennDiagram)

data("TCGAexp")
data("TCGAmet")

## -----------------------------------------------------------------------------
dim(TCGAexp)
ifelse(checkPairing(TCGAexp, TCGAmet),
    "Match OK", "Check matching"
)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(MultiAssayExperiment)
doubleExp <- list(
    "methylation" = TCGAmet,
    "expression" = TCGAexp
)
mae1 <- MultiAssayExperiment(experiments = doubleExp)

## -----------------------------------------------------------------------------
cl <- correlationSelection(mae1, type = "Spearman",
    pValCutoff = 0.25, rCutoff = -0.5, adj = TRUE)

correlationL <- cl[!is.na(cl$SigNegCorr) & 
    cl$SigNegCorr, ]
correlationNoL <- cl[!is.na(cl$SigNegCorr)& 
    !cl$SigNegCorr, ]

message(
    "The number of genes selected 
    with the correlation method is: ", sum(correlationL$SigNegCorr),
    "\n"
)

## ----fig.width = 6, fig.height=4----------------------------------------------
d <- density(correlationL[, 1])
x2 <- data.frame(x = d$x, y = d$y)

library(ggplot2)
ggplot(x2, aes(x,y)) + geom_line() +
    labs(
    title = "Significant correlations in the TCGA dataset",
    x = "Correlation",
    y = "Density"
    )+
    theme_minimal()


## -----------------------------------------------------------------------------
head(correlationL)

## ----fig.height=6, fig.width=6------------------------------------------------
genes2plot <- rownames(correlationL)[1:4]

plotGenesMat(mae1, geneNames = genes2plot,
    fileName = NULL, text4Title = correlationL[rownames(correlationL), ""]
)

## -----------------------------------------------------------------------------
sampleSize <- dim(mae1[[2]])[2]
numGenes <- dim(mae1[[2]])[1]

reqPercentages <- matrix(c(2, 20, 5, 5, 40, 20, 3, 3, 2), 
    nrow = 3, byrow = TRUE)
sum(reqPercentages)
(maxminCounts <- toReqMat(sampleSize, reqPercentages))

(theWeightMifL <- matrix(c(2, -2, -sampleSize / 5, 1, 0, -2, 1, 1, 2), 
    nrow = 3, byrow = TRUE))
(theWeightMifNonL <- matrix(c(0, -2, -sampleSize / 5, 0, 0, -2, 0, 0, 0),
    nrow = 3,
    byrow = TRUE
))

heur <- scoreGenesMat(mae1,
    aReqPercentsMat = reqPercentages, aWeightMifL = theWeightMifL,
    aWeightMifNonL = theWeightMifNonL
)

message("Number of scatterplots scored  : ", dim(heur)[1], "\n")
message("Number of L-shape scatterplots : ", sum(heur[, 1]), "\n")

heurL <- heur[heur$logicSc, ]
heurNoL <- heur[!heur$logicSc, ]

## -----------------------------------------------------------------------------
knitr::kable(heurL)

## ----fig.height=6, fig.width=6------------------------------------------------
genes2plot2 <- rownames(mae1[[2]]) %in% rownames(heurL)[1:4]

plotGenesMat(mae1,
    fileName = NULL, text4Title = heurL[genes2plot2, "numeriSc"]
)


## -----------------------------------------------------------------------------
inCommonL <- intersect(rownames(correlationL), 
    rownames(heurL))
inCorrelationLOnly <- setdiff(rownames(correlationL),
    inCommonL)
inheurLLOnly <- setdiff(rownames(heurL), inCommonL)

## ----fig.height=6, fig.width=6------------------------------------------------
par(mfrow = c(2, 2))
myGene1 <- inCommonL[1]
titleT <- paste(myGene1, "(May be GRM)")
plotGeneSel(mae1,myGene1,
    titleText = titleT,
    x1 = 1/3, x2 = 2/3)

myGene2 <- inCommonL[2]
titleT <- paste(myGene2, "(May be GRM)")
plotGeneSel(mae1, myGene2,
    titleText = titleT,
    x1 = 1/3, x2 = 2/3)


myGene3 <- inCommonL[3]
titleT <- paste(myGene3, "(May be GRM)")
plotGeneSel(mae1,myGene3,
    titleText = titleT,
    x1 = 1/3, x2 = 2/3)

myGene5 <- inCommonL[5]
titleT <- paste(myGene5, "(May be GRM)")
plotGeneSel(mae1, myGene5,
    titleText = titleT,
    x1 = 1/3, x2 = 2/3)

## -----------------------------------------------------------------------------
sessionInfo()

