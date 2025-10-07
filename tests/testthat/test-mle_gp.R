test_that("mle_gp() returns expected output", {
  result <- mle_gp(y=1:5, x=matrix(1:10, ncol=2))
  expect_type(result, "list")
  expect_true("theta" %in% names(result))
})
