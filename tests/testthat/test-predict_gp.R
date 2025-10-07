test_that("predict_gp() returns expected output", {
  x.test <- seq(1, 10, length=6)
  result <- mle_gp(y=1:5, x=matrix(c(1,3,4,6,7), ncol=1))
  pred <- predict_gp(result, x.test)
  expect_true("mup" %in% names(pred))
})
