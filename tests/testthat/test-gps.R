test_that("Error is thrown when lengths of u and v differ", {
  expect_error(gps(c(0.1, 0.1), c(0.1)), "Lengths of u and v differ")
})

test_that("Error is thrown when indices of max u and v elements coincide", {
  expect_error(gps(c(0.1, 0.2, 0.3), c(0.1, 0.15, 0.3)), "Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.")
})
