test_that("cv returns expected output", {
  y <- rnorm(3)
  x <- matrix(1:3, ncol=1)
  lambda <- 0.1
  initialvals <- matrix(c(0.3, 0.5), ncol=1)
  n <- length(y)
  k <- length(y)
  p <- 1
  d <- 1
  results <- cv(y=y, x=x, lambda=lambda, initialvals=initialvals, n=n, k=k, p=p,
                d=d, fixed_g=0.0000001)
  expect_true("mse_min" %in% names(results))
})
