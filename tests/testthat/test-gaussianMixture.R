test_that("Sample is drawn without error in simple use case", {
  set.seed(42)
  sam<-rmixture(1e2, 0.1, 1, 2)
})
