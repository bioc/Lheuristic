library(MultiAssayExperiment)
library(testthat)

test_that("lhCreateMAE creates a valid MultiAssayExperiment", {
  
  # Create synthetic methylation and expression data
  methylData <- matrix(runif(50), nrow = 10)
  colnames(methylData) <- paste0("samp", 1:ncol(methylData))
  rownames(methylData) <- paste0("gene", 1:nrow(methylData))
  
  expresData <- matrix(rnorm(50), nrow = 10)
  colnames(expresData) <- colnames(methylData)
  rownames(expresData) <- rownames(methylData)
  
  # Create sample metadata
  colDat <- data.frame(sampleID = colnames(methylData),
                       name = letters[1:ncol(methylData)])
  rownames(colDat) <- colDat$sampleID
  
  # Call function
  mae <- lhCreateMAE(methylData, expresData, colData = colDat)
  
  # Check that the output is a MultiAssayExperiment object
  expect_s4_class(mae, "MultiAssayExperiment")
  
  # Check that dataset names match the expected ones
  expect_equal(names(experiments(mae)), c("methylation", "expression"))
  
  # Check if colData was correctly assigned
  expect_equal(rownames(colData(mae)), colDat$sampleID)
  
  # Check if experiments were correctly stored
  expect_equal(dim(experiments(mae)[["methylation"]]), dim(methylData))
  expect_equal(dim(experiments(mae)[["expression"]]), dim(expresData))
})

test_that("lhCreateMAE works when row/column names match", {
    methylData <- matrix(runif(50), nrow = 10)
    colnames(methylData) <- paste0("samp", 1:ncol(methylData))
    rownames(methylData) <- paste0("gene", 1:nrow(methylData))
    expresData <- matrix(rnorm(50), nrow = 10)
    colnames(expresData) <- colnames(methylData)  # Match column names
    rownames(expresData) <- rownames(methylData)

    # Suppress messages and ensure no error occurs
    expect_silent(suppressMessages(lhCreateMAE(methylData, expresData)))
    })

test_that("lhCreateMAE assigns default colData if missing", {
  methylData <- matrix(runif(50), nrow = 10)
  colnames(methylData) <- paste0("samp", 1:ncol(methylData))
  rownames(methylData) <- paste0("gene", 1:nrow(methylData))
  
  expresData <- matrix(rnorm(50), nrow = 10)
  colnames(expresData) <- colnames(methylData)
  rownames(expresData) <- rownames(methylData)
  
  mae <- lhCreateMAE(methylData, expresData)  # No colData provided
  
  expect_s4_class(mae, "MultiAssayExperiment")
  expect_equal(rownames(colData(mae)), colnames(methylData))
})
