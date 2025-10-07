test_that("multiplication works", {
  y.test <- seq(11, 20, length=10)
  x.test <- seq(1, 10, length=10)
  result <- mle_gp(y=1:5, x=matrix(c(1,3,4,6,7), ncol=1))
  pred <- predict_gp(result, x.test)
  score_value <- score(y.test, pred$mup, pred$Sigmap)
  expect_type(score_value, "double")
  expect_length(score_value, 1)
})
