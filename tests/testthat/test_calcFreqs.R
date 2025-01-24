library(MultiAssayExperiment)

# Methylation data
methylData <- matrix(runif(50), nrow = 10)
colnames(methylData) <- paste0("samp", 1:ncol(methylData))
rownames(methylData) <- paste0("gene", 1:nrow(methylData))
# Expression data
expresData <- matrix(rnorm(50), nrow = 10)
colnames(expresData) <- paste0("samp", 1:ncol(methylData))
rownames(expresData) <- paste0("gene", 1:nrow(methylData))
# ColData
colDat <- data.frame(sampleID = colnames(methylData), name = letters[1:ncol(methylData)])
rownames(colDat) <- colDat$sampleID


mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(
    methylation = methylData,
    expression = expresData
), colData = colDat)
geneRow <- 1
x1 <- 1 / 3
x2 <- 2 / 3
y1 <- NULL
y2 <- NULL
percY1 <- 1 / 3
percY2 <- 2 / 3

calcFreqs(mae = mae, geneNum = geneRow, x1, x2, y1, y2, percY1, percY2)
colData(mae)
experiments(mae)
sampleMap(mae)
prod(names(experiments(mae)) == c("methylation", "expression")) == 1
