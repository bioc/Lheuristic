library(testthat)
library(MultiAssayExperiment)

test_that("scoreGenesMat generates valid scores", {
  # Generate test data
  set.seed(42)
  methylData <- matrix(runif(50), nrow = 10)
  colnames(methylData) <- paste0("samp", 1:ncol(methylData))
  rownames(methylData) <- paste0("gene", 1:nrow(methylData))
  
  expresData <- matrix(rnorm(50), nrow = 10)
  colnames(expresData) <- paste0("samp", 1:ncol(methylData))
  rownames(expresData) <- paste0("gene", 1:nrow(methylData))
  
  colDat <- data.frame(sampleID = colnames(methylData), name = letters[1:ncol(methylData)])
  rownames(colDat) <- colDat$sampleID
  
  mae <- lhCreateMAE(methylData,
      expresData, colData = colDat
  )
  
  sampleSize <- ncol(MultiAssayExperiment::experiments(mae)[[1]])
  reqPercentages <- matrix(c(3, 20, 5, 5, 40, 20, 4, 1, 2),
                           nrow = 3, byrow = TRUE
  )
  
  theWeightMifL <- matrix(c(2, -2, -sampleSize / 5, 1, 0, -2, 1, 1, 2),
                          nrow = 3, byrow = TRUE
  )
  
  theWeightMifNonL <- matrix(c(0, -2, -sampleSize / 5, 0, 0, -2, 0, 0, 0),
                             nrow = 3, byrow = TRUE
  )
  
  # Run the function
  scores <- scoreGenesMat(mae,
                          x1 = 1 / 3, x2 = 2 / 3,
                          y1 = NULL, y2 = NULL, percY1 = 1 / 3, percY2 = 2 / 3,
                          aReqPercentsMat = reqPercentages,
                          aWeightMifL = theWeightMifL,
                          aWeightMifNonL = theWeightMifNonL
  )
  
  # Check that the output is a data frame
  expect_s3_class(scores, "data.frame")
  
  # Check that it has the correct column names
  expect_true(all(c("logicSc", "numericSc") %in% colnames(scores)))
  
  # Check that logicSc is logical (TRUE/FALSE)
  expect_type(scores$logicSc, "logical")
  
  # Check that numericSc is numeric
  expect_type(scores$numericSc, "double")
  
  # Check that the number of rows matches the number of genes
  expect_equal(nrow(scores), nrow(methylData))
  
  # Ensure numeric scores are finite values
  expect_true(all(is.finite(scores$numericSc)))
})
