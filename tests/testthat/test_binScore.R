# Test that binScore returns the expected output for a given input
library(testthat)
test_that("binScore returns the expected output for a given input", {
    # Generate example data
    aGrid <- matrix(c(20, 3, 0, 10, 2, 2, 20, 10, 20), nrow = 3, ncol = 3, byrow = TRUE)
    aReq <- matrix(c(15, 5, 0, 0, 5, 5, 10, 10, 15), nrow = 3, ncol = 3, byrow = TRUE)

    # Test binScore with the example data
    expected_output <- TRUE
    output <- binScore(aGrid, aReq)
    expect_equal(output, expected_output, info = "binScore should return TRUE for the example input data")
})
