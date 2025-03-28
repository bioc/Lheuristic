library(testthat)
library(MultiAssayExperiment)

# Create example data
set.seed(123)
methylData <- matrix(runif(50), nrow = 10)
colnames(methylData) <- paste0("samp", 1:ncol(methylData))
rownames(methylData) <- paste0("gene", 1:nrow(methylData))

expresData <- matrix(rnorm(50), nrow = 10)
colnames(expresData) <- paste0("samp", 1:ncol(methylData))
rownames(expresData) <- paste0("gene", 1:nrow(methylData))

colDat <- data.frame(sampleID = colnames(methylData), name = letters[1:ncol(methylData)])
rownames(colDat) <- colDat$sampleID

mae <- lhCreateMAE(methylData, expresData, colData = colDat)

test_that("calcFreqs returns a 3x3 matrix with expected properties", {
  # Set parameters
    geneRow <- 1
    x1 <- 1 / 3
    x2 <- 2 / 3
    y1 <- NULL
    y2 <- NULL
    percY1 <- 1 / 3
    percY2 <- 2 / 3
    
    # Run function
    freqs <- calcFreqs(mae, geneRow, x1, x2, y1, y2, percY1, percY2)
  
  # Tests
  expect_type(freqs, "double")  # Expect numeric matrix
  expect_equal(dim(freqs), c(3, 3))  # Expect 3x3 matrix
  expect_true(all(freqs >= 0))  # Expect non-negative values
  expect_equal(sum(freqs), ncol(methylData))  # Total count should match the number of samples
})

test_that("calcFreqs handles out-of-bounds geneNum", {

  numGenes <- nrow(MultiAssayExperiment::assay(mae, 1))  # Get number of genes
  
  # Test valid geneNum (should not throw errors)
  expect_silent(calcFreqs(mae, 1, 1/3, 2/3))  # First row
  expect_silent(calcFreqs(mae, numGenes, 1/3, 2/3))  # Last row
  
  # Test out-of-bounds geneNum
  warnings <- testthat::capture_warnings(calcFreqs(mae, 0, 1/3, 2/3))
  
  # Check that the specific warnings appear
  expect_true(any(grepl("no non-missing arguments to min", warnings)))
  expect_true(any(grepl("no non-missing arguments to max", warnings)))

  expect_error(calcFreqs(mae, numGenes + 1, 1/3, 2/3), regexp = "subscript out of bounds")
  
  # Alternative: Catch any generic warnings without expecting a specific message
  # expect_warning(calcFreqs(mae, 0, 1/3, 2/3))
  # expect_error(calcFreqs(mae, numGenes + 1, 1/3, 2/3))
})

