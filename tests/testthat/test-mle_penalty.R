test_that("mle_penalty() returns expected output", {
  skip_on_cran()
  result <- gp_cv(y=1:5, x=matrix(c(1,3,4,6,7), ncol=1), lambda=0.5)
  out <- mle_penalty(result)
  expect_true("theta" %in% names(out))
})
