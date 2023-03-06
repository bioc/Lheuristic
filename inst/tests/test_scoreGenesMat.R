expect_identical(100, sum(aReqPercentsMat))
expect_error(sum(aReqPercentsMat) != 100, "Percentages must add up to 100")

