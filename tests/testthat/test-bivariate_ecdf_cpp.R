set.seed(42)

test_that("Test bivariate_ecdf_cpp on some trivial cases", {
  testSample <- seq(0.1, 1.0, length.out=10)
  expect_equal(bivariate_ecdf_cpp(testSample, testSample), testSample)

  testSample <- c(1.0, 0.9, 0.5, 0.2, 0.7, 0.1, 0.6, 0.3, 0.4, 0.8)
  expect_equal(bivariate_ecdf_cpp(testSample, testSample), testSample)
})

test_that('Test bivariate_ecdf_cpp with duplicates in the sample', {
  testSample <- c(0.1,0.2,0.2,0.3,0.4)
  expect_equal(bivariate_ecdf_cpp(testSample, testSample), c(0.2,0.6,0.6,0.8,1.0))
})
