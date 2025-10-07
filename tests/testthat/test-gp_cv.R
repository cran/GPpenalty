test_that("gp_cv() returns expected output ", {
  skip_on_cran()
  result <- gp_cv(y=1:5, x=matrix(1:10, ncol=2))
  expect_true("lambda.min" %in% names(result))
})
