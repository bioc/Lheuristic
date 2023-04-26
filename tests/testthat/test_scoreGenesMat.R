library(testthat)# Define the test data
set.seed(123)
mets <- matrix(runif(100), nrow = 10)
expres <- matrix(rnorm(100), nrow = 10)
reqPercentages <- matrix(c(3, 20, 5, 5, 40, 20, 4, 1, 2), nrow = 3, byrow = TRUE)
theWeightMifL <- matrix(c(2, -2, -5, 1, 0, -2, 1, 1, 2), nrow = 3, byrow = TRUE)
theWeightMifNonL <- matrix(c(0, -2, -5, 0, 0, -2, 0, 0, 0), nrow = 3, byrow = TRUE)

# Define the test
test_that("scoreGenesMat returns a data frame with two columns", {
  # Call the function
  result <- scoreGenesMat(mets = mets, expres = expres,
                          aReqPercentsMat = reqPercentages,
                          aWeightMifL = theWeightMifL,
                          aWeightMifNonL = theWeightMifNonL)
  
  # Check that the result is a data frame with two columns
  expect_is(result, "data.frame")
  expect_equal(ncol(result), 2)
  
  # Check that the column names are "logicSc" and "numericSc"
  expect_identical(colnames(result), c("logicSc", "numericSc"))
})


library(testthat)

context("scoreGenesMat")

test_that("sum of percentages in aReqPercentsMat should be 100", {
  mets <- matrix(runif(1000), nrow=100)
  expres <- matrix(rnorm(1000), nrow=100)
  reqPercentages <- matrix(c(3, 20, 5, 5, 40, 20, 4, 1, 2), nrow = 3, byrow = TRUE)
  expect_error(scoreGenesMat(mets, expres, aReqPercentsMat = reqPercentages + 1))
  expect_error(scoreGenesMat(mets, expres, aReqPercentsMat = reqPercentages - 1))
  expect_no_error(scoreGenesMat(mets, expres, aReqPercentsMat = reqPercentages))
})
