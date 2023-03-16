aGrid <- matrix(c(20, 3, 0, 10, 2, 2, 20, 10, 20), nrow = 3, ncol = 3, byrow = TRUE)
aReq <- matrix(c(15, 5, 0, 0, 5, 5, 10, 10, 15), nrow = 3, ncol = 3, byrow = TRUE)
binSc <- binScore(aGrid, aReq)

test_that("L shape score", {
  expect_that(binSc , is.logical)
})