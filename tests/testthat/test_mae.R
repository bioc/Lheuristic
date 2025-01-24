library(MultiAssayExperiment)

# Methylation data

methylData <- matrix(runif(50), nrow = 10)
colnames(methylData) <- paste0("samp", 1:ncol(methylData))
rownames(methylData) <- paste0("gene", 1:nrow(methylData))


# Expression data

expresData <- matrix(rnorm(50), nrow = 10)
colnames(expresData) <- paste0("samp", 1:ncol(methylData))
rownames(expresData) <- paste0("gene", 1:nrow(methylData))

# Information about samples goes to 'ColData'

colDat <- data.frame(sampleID = colnames(methylData), name = letters[1:ncol(methylData)])
rownames(colDat) <- colDat$sampleID

# Information about features goes to 'RowData' It seems this property is
# available in SummarizedExperiments, but not in MultiAssayExperiments

# The main reason seems to be because `MultiAssayExperiment` requires the data
# slots to have the same column numbers & names but allows for distinct rows
# numbers/names Check QFeaturesClass

# rowDat <- data.frame(featureID = rownames(methylData), isOncogene
# =c(rep(1,3), rep(0,7)), known2Methylate =sample(x= c(0,1), size= 10,
# rep=TRUE, prob=c(.8,.2))) rownames(rowDat)<- rowDat$featureID


mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(
    methylation = methylData,
    expression = expresData
), colData = colDat)

mae@colData
mae@ExperimentList
mae@sampleMap

MultiAssayExperiment::experiments(mae)[[1]]
mae@ExperimentList[[1]]
MultiAssayExperiment::experiments(mae)[[2]]
