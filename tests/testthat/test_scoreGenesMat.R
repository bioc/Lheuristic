aReqPercentsMat1 <- matrix (c(3, 20, 5, 5, 40, 20, 4, 1, 2), nrow=3, byrow=TRUE)

test_that("Sum percentage matrix" ,{
  expect_identical(100, sum(aReqPercentsMat1))})

aReqPercentsMat2 <- matrix (c(3, 15, 5, 5, 40, 20, 4, 1, 2), nrow=3, byrow=TRUE)
test_that("Sum percentage matrix", {
  expect_error(sum(aReqPercentsMat2) != 100, "Percentages must add up to 100")})