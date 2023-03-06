
test_that("Sum percentage matrix" ,{
  expect_identical(100, sum(aReqPercentsMat))})

test_that("Sum percentage matrix", {
  expect_error(sum(aReqPercentsMat) != 100, "Percentages must add up to 100")})

# test_that("L shape score", {
#   expect_equal(binSc , "FALSE")
#   expect_equal(binSc, "TRUE")
# })